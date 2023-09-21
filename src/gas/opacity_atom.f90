!No dissolution and no mag yet!
module Opacity_atom

   use atom_type
   use grid
   use parametres
   use broad, Only               : line_damping
   use voigts, only              : voigt
   use occupation_probability, only : f_dissolve
   use gas_contopac, only        : H_bf_Xsection, alloc_gas_contopac, background_continua_lambda, &
                                     dealloc_gas_contopac, hnu_k
   use wavelengths, only         :  n_lambda
   use wavelengths_gas, only     : Nlambda_max_line, Nlambda_max_cont, n_lambda_cont, tab_lambda_cont, tab_lambda_nm, &
                                    peak_gauss_limit, Nlambda_line_gauss_log, Nlambda_line_gauss_lin
   use constantes, only          : c_light
   use molecular_emission, only  : v_proj, ds
   use utils, only               : linear_1D_sorted
   !$ use omp_lib

   implicit none

   real(kind=dp), allocatable, dimension(:,:) :: chi_cont, eta_cont
   !local profile for cell id in direction iray for all atoms and b-b trans
   real(kind=dp), allocatable :: Itot(:,:,:), psi(:,:,:), phi_loc(:,:,:,:,:), vlabs(:,:)
   real(kind=dp), allocatable :: eta_atoms(:,:,:), Uji_down(:,:,:,:), chi_up(:,:,:,:), chi_down(:,:,:,:), chi_tot(:), eta_tot(:)
   integer, parameter 		   :: NvspaceMax = 151
   integer, dimension(:), allocatable :: Nvel_points
   real(kind=dp), dimension(:,:), allocatable :: omegav
   logical 		               :: lnon_lte_loop
   integer                    :: N_gauss

   contains

   function gmax_line(line)
   !compute the damping max of a given line line.
   !assumes that eletronic densities, populations and thermodynamics
   !quantities are set
      type (AtomicLine), intent(in) :: line
      integer :: i
      real(kind=dp) :: gmax_line
   !could be para
      gmax_line = 0.0_dp ! damping * vth
      do i=1,n_cells
         if (icompute_atomRT(i)>0.0) then
            gmax_line = max(gmax_line,line_damping(i,line)*&
                        vbroad(T(i),line%atom%weight, vturb(i)))
         endif
      enddo
      return
   end function gmax_line
   function gmin_line(line)
   !compute the damping min of a given line line.
   !assumes that eletronic densities, populations and thermodynamics
   !quantities are set
      type (AtomicLine), intent(in) :: line
      integer :: i
      real(kind=dp) :: gmin_line
   !could be para
      gmin_line = 1d50 ! damping * vth
      do i=1,n_cells
         if (icompute_atomRT(i)>0.0) then
            gmin_line = min(gmin_line,line_damping(i,line)*&
                        vbroad(T(i),line%atom%weight, vturb(i)))
         endif
      enddo
      return
   end function gmin_line

   subroutine set_max_damping()
   !sets also the minimum
   !result in m/s
      integer :: nat, kr
      do nat=1, n_atoms
         do kr=1,atoms(nat)%p%nline
            if (.not.atoms(nat)%p%lines(kr)%lcontrib) cycle
            if (atoms(nat)%p%lines(kr)%Voigt) then
               atoms(nat)%p%lines(kr)%damp_max = gmax_line(atoms(nat)%p%lines(kr))
               atoms(nat)%p%lines(kr)%damp_min = gmin_line(atoms(nat)%p%lines(kr))
            endif
         enddo
      enddo
      return
   end subroutine set_max_damping

   subroutine deactivate_lines()
      integer :: nact, kr
      do nact=1, N_Atoms
         do kr=1, Atoms(nact)%p%Nline
            atoms(nact)%p%lines(kr)%lcontrib = .false.
         enddo
      enddo
      return 
   end subroutine deactivate_lines
   subroutine activate_lines()
      integer :: nact, kr
      do nact=1, N_Atoms
         do kr=1, Atoms(nact)%p%Nline
            atoms(nact)%p%lines(kr)%lcontrib = .true.
         enddo
      enddo
      return 
   end subroutine activate_lines
   subroutine deactivate_continua()
      integer :: nact, kr
      do nact=1, N_Atoms
         do kr=1, Atoms(nact)%p%Ncont
            atoms(nact)%p%continua(kr)%lcontrib = .false.
         enddo
      enddo
      return 
   end subroutine deactivate_continua
   subroutine activate_continua()
      integer :: nact, kr
      do nact=1, N_Atoms
         do kr=1, Atoms(nact)%p%ncont
            atoms(nact)%p%continua(kr)%lcontrib = .true.
         enddo
      enddo
      return 
   end subroutine activate_continua

   subroutine setup_activelines_list()
   !TO DO: smarter way
   !TO DO: add continuum
      integer :: l, nat, kr, nb, nr, ntr
      real(kind=dp) :: dl
      type (AtomType), pointer :: at

      ! allocate(art_lines(n_lambda))
      ! allocate(art_lines(1)%art(1))
      ! art_lines(1)%art(1)%p => hydrogen%lines(1)

      ! write(*,*) art_lines(1)%art(1)%p%lambda0, hydrogen%lines(1)%lambda0
      ! stop

      allocate(art_lines(n_lambda))
      !first, count number of transitions active at a given frequency
      !could be better
      do l=1, n_lambda
         art_lines(l)%Ntrans = 0
         do nat=1, n_atoms
            at => atoms(nat)%p
            tr_loop : do kr = 1,at%Nline
               if (.not.at%lines(kr)%lcontrib) cycle tr_loop
               nb = at%lines(kr)%Nb; nr = at%lines(kr)%Nr
               dl = max(at%lines(kr)%lambda0-tab_lambda_nm(Nb),tab_lambda_nm(Nr)-at%lines(kr)%lambda0)
               if (abs(tab_lambda_nm(l) - at%lines(kr)%lambda0) <= dl) then
                     art_lines(l)%Ntrans = art_lines(l)%Ntrans + 1
               endif
            enddo tr_loop
         enddo
         ! write(*,*) " There are", art_lines(l)%Ntrans," active lines at that lambda", tab_lambda_nm(l)
      enddo

      !then setup the pointers
      do l=1, n_lambda
         allocate(art_lines(l)%art(art_lines(l)%Ntrans))
         ntr = 0
         do nat=1, n_atoms
            at => atoms(nat)%p
            tr_loopb : do kr = 1,at%Nline
               if (.not.at%lines(kr)%lcontrib) cycle tr_loopb
               nb = at%lines(kr)%Nb; nr = at%lines(kr)%Nr
               dl = max(at%lines(kr)%lambda0-tab_lambda_nm(Nb),tab_lambda_nm(Nr)-at%lines(kr)%lambda0)
               if (abs(tab_lambda_nm(l) - at%lines(kr)%lambda0) <= dl) then
                     ntr = ntr + 1
                     art_lines(l)%art(ntr)%p => at%lines(kr)
               endif
            enddo tr_loopb
         enddo
      enddo

      return
   endsubroutine setup_activelines_list

   !could be parralel
   subroutine alloc_atom_opac(N,x,limage)
      integer, intent(in) :: N
      real(kind=dp), dimension(N) :: x
      logical, intent(in) :: limage
      integer :: nat, kr, icell, nb, nr
      type(AtomType), pointer :: atm
      real(kind=dp) :: vth
      integer(kind=8) :: mem_loc, mem_contopac

      mem_loc = 0

      ! call alloc_gas_contopac(N,x)
      !-> on a small grid and interpolated later
      call alloc_gas_contopac(n_lambda_cont,tab_lambda_cont)
      if (limage) then
         N_gauss = 0
         do nat=1, n_Atoms
            N_gauss = max(N_gauss, atoms(nat)%p%n_speed_rt)
         enddo
      else
         N_gauss = (Nlambda_line_gauss_log + Nlambda_line_gauss_lin) + mod(Nlambda_line_gauss_log + Nlambda_line_gauss_lin+1,2)
      endif

      allocate(omegav(NvspaceMax,nb_proc),Nvel_points(nb_proc)); omegav = 0.0_dp
      mem_loc = mem_loc + sizeof(omegav) + sizeof(Nvel_points)

      call setup_activelines_list()
      mem_loc = mem_loc + sizeof(art_lines)
      ! write(*,*) "size art_lines:", sizeof(art_lines) / 1024.**3, 'GB'

      do nat=1, n_atoms
         atm => atoms(nat)%p
         mem_loc = mem_loc + sizeof(atm%vth)
         do kr=1,atm%nline
               if (.not.atm%lines(kr)%lcontrib) cycle
               if (atm%lines(kr)%Voigt) then
               !-> do not allocate if using thomson and humlicek profiles
                  ! allocate(atm%lines(kr)%v(atm%lines(kr)%Nlambda),atm%lines(kr)%phi(atm%lines(kr)%Nlambda,n_cells))
                  allocate(atm%lines(kr)%a(n_cells)); atm%lines(kr)%a(:) = 0.0_dp
                  allocate(atm%lines(kr)%b(n_cells),atm%lines(kr)%c(n_cells))
                  mem_loc = mem_loc + 2 * sizeof(atm%lines(kr)%b)
                  mem_loc = mem_loc + sizeof(atm%lines(kr)%a)!+sizeof(atm%lines(kr)%phi)+sizeof(atm%lines(kr)%v)
               endif
         enddo

         do icell=1, n_cells
            if (icompute_atomRT(icell) <= 0) cycle
            !tmp because of vbroad, recomputed after in m/s
            vth = vbroad(T(icell),atm%weight, vturb(icell))
            !Voigt
            do kr=1,atm%nline
               if (.not.atm%lines(kr)%lcontrib) cycle
               if (atm%lines(kr)%Voigt) then
                  nb = atm%lines(kr)%nb; nr = atm%lines(kr)%nr
                  atm%lines(kr)%a(icell) = line_damping(icell,atm%lines(kr))
               !-> do not allocate if using thomson and humlicek profiles
                  ! atm%lines(kr)%v(:) = c_light * (x(nb:nr)-atm%lines(kr)%lambda0)/atm%lines(kr)%lambda0 / vth
                  ! atm%lines(kr)%phi(:,icell) = Voigt(atm%lines(kr)%Nlambda, &
                  !                                  atm%lines(kr)%a(icell), &
                  !                                  atm%lines(kr)%v(:)) / (vth * sqrtpi)
                  !f Thomson
                  atm%lines(kr)%b(icell) = (vth**5. + 2.69269*vth**4. * (atm%lines(kr)%a(icell)*vth) + 2.42843*vth**3. * (atm%lines(kr)%a(icell)*vth)**2. + &
                     4.47163*vth**2. *(atm%lines(kr)%a(icell)*vth)**3. + 0.07842*vth*(atm%lines(kr)%a(icell)*vth)**4. + (atm%lines(kr)%a(icell)*vth)**5.)**(0.2)

                  !eta Thomson
                  atm%lines(kr)%c(icell) = 1.36603*(vth*atm%lines(kr)%a(icell)/atm%lines(kr)%b(icell)) - &
                     0.47719*(vth*atm%lines(kr)%a(icell)/atm%lines(kr)%b(icell))**2 + &
                     0.11116*(vth*atm%lines(kr)%a(icell)/atm%lines(kr)%b(icell))**3
               endif
            enddo
         enddo

         !-> do not allocate if using thomson and humlicek profiles
         ! do kr=1,atm%nline
         !    if (.not.atm%lines(kr)%lcontrib) cycle
         !    if (atm%lines(kr)%Voigt) then
         !       nb = atm%lines(kr)%nb; nr = atm%lines(kr)%nr
         !       !-> will not change during the non-LTE loop.
         !       atm%lines(kr)%v(:) = c_light * (x(nb:nr)-atm%lines(kr)%lambda0)/atm%lines(kr)%lambda0 !m/s
         !    endif
         ! enddo

         do kr = 1, atm%Ncont
            if (.not.atm%continua(kr)%lcontrib) cycle
            ! allocate(atm%continua(kr)%twohnu3_c2(atm%continua(kr)%Nlambda))
            ! atm%continua(kr)%twohnu3_c2(:) = twohc/x(atm%continua(kr)%Nb:atm%continua(kr)%Nr)**3
            ! allocate(atm%continua(kr)%alpha(atm%continua(kr)%Nlambda))
            ! if (atm%continua(kr)%hydrogenic) then
            !    atm%continua(kr)%alpha(:) = H_bf_Xsection(atm%continua(kr), x(atm%continua(kr)%Nb:atm%continua(kr)%Nr))
            ! else
            !    atm%continua(kr)%alpha(:) = linear_1D_sorted(size(atm%continua(kr)%alpha_file),&
            !    atm%continua(kr)%lambda_file,atm%continua(kr)%alpha_file,atm%continua(kr)%Nlambda,x(atm%continua(kr)%Nb:atm%continua(kr)%Nr))
            ! endif
            !-> on a small grid and interpolated later
            allocate(atm%continua(kr)%twohnu3_c2(atm%continua(kr)%Nlambdac))
            atm%continua(kr)%twohnu3_c2(:) = twohc / &
                 tab_lambda_cont(atm%continua(kr)%Nbc:atm%continua(kr)%Nrc)**3
            allocate(atm%continua(kr)%alpha(atm%continua(kr)%Nlambdac))
            if (atm%continua(kr)%hydrogenic) then
               atm%continua(kr)%alpha(:) = H_bf_Xsection(atm%continua(kr), &
                    tab_lambda_cont(atm%continua(kr)%Nbc:atm%continua(kr)%Nrc))
            else
               atm%continua(kr)%alpha(:) = linear_1D_sorted(size(atm%continua(kr)%alpha_file),&
                    atm%continua(kr)%lambda_file,atm%continua(kr)%alpha_file,atm%continua(kr)%Nlambdac,&
                    tab_lambda_cont(atm%continua(kr)%Nbc:atm%continua(kr)%Nrc))
               atm%continua(kr)%alpha(atm%continua(kr)%Nlambdac) = atm%continua(kr)%alpha_file(size(atm%continua(kr)%alpha_file))
               !to chheck the edge
            endif

         enddo

         atm => null()
      enddo

      !now if not limit_mem, store the initial value of the continuum opacities
      !for image mode or for 1st iteration of the non-LTE loop;
      mem_contopac = 0
      if (.not.llimit_mem) then
         allocate(chi_cont(n_lambda_cont,n_cells), eta_cont(n_lambda_cont,n_cells))
         mem_contopac = 2 * sizeof(chi_cont)
         mem_loc = mem_loc + mem_contopac
         call calc_contopac()!needs continua(:)%alpha etc
      endif

      write(*,'("allocate "(1F6.2)" GB for background opacities")') real(mem_loc) / 1024.0/1024./1024.0
      if (mem_contopac>0) write(*,'("  -> with "(1F6.2)" GB for contopac.")') real(mem_contopac) / 1024.0/1024./1024.0

      return
   end subroutine alloc_atom_opac

   subroutine dealloc_atom_opac()
      integer :: nat, kr
      type(AtomType), pointer :: atm

      call dealloc_gas_contopac()
      if (allocated(chi_cont)) deallocate(chi_cont,eta_cont)

      if (allocated(omegav)) deallocate(omegav,Nvel_points)
      if (allocated(art_lines)) deallocate(art_lines)

      do nat=1, n_atoms
         atm => atoms(nat)%p

         do kr=1,atm%nline
            if (allocated( atm%lines(kr)%a)) deallocate(atm%lines(kr)%a )
            if (allocated( atm%lines(kr)%v)) deallocate(atm%lines(kr)%v,atm%lines(kr)%phi)
            if (allocated( atm%lines(kr)%b)) deallocate(atm%lines(kr)%b,atm%lines(kr)%c)
         enddo

         do kr = 1, atm%Ncont
            if (allocated(atm%continua(kr)%twohnu3_c2)) then
               deallocate(atm%continua(kr)%twohnu3_c2)
               deallocate(atm%continua(kr)%alpha)
            endif
         enddo

         atm => null()
      enddo


      return
   end subroutine dealloc_atom_opac

   subroutine calc_contopac_loc(icell)
      integer, intent(in) :: icell

      call background_continua_lambda(icell, n_lambda_cont, tab_lambda_cont, chi_cont(:,icell), eta_cont(:,icell))
      call opacity_atom_bf_loc(icell, n_lambda_cont, tab_lambda_cont, chi_cont(:,icell), eta_cont(:,icell))

      return
   end subroutine calc_contopac_loc

   subroutine calc_contopac()
      integer :: icell, ibar, n_cells_done

      write(*,*) " ** Init continuous background opacities..."
      ibar = 0
      n_cells_done = 0

      !$omp parallel &
      !$omp default(none) &
      !$omp private(icell) &
      !$omp shared(ibar,n_cells_done,n_cells,icompute_atomrt)
      !$omp do schedule(dynamic,1)
      do icell=1, n_Cells
         if (icompute_atomRT(icell)>0) then
            call calc_contopac_loc(icell)
            !$omp atomic
            n_cells_done = n_cells_done + 1
            if (real(n_cells_done) > 0.02*ibar*n_cells) then
               call progress_bar(ibar)
               !$omp atomic
               ibar = ibar+1
            endif
         endif
      enddo
      !$omp end do
      !$omp end parallel
      call progress_bar(50)
      !write(*,*) " " !for progress bar

      return
   end subroutine calc_contopac

   subroutine contopac_atom_loc(icell,N,lambda,chi,snu)
      integer, intent(in) :: icell, N
      real(kind=dp), intent(in), dimension(N) :: lambda
      real(kind=dp), intent(inout), dimension(N) :: chi, Snu
      real(kind=dp), dimension(N_lambda_cont) :: chic, snuc
      integer :: la, lac, i0
      real(kind=dp) :: w

      if (llimit_mem) then
         !init continuous opacity with background gas continuum.
         call background_continua_lambda(icell, n_lambda_cont, tab_lambda_cont, chic, Snuc)
         !Snu = Snu + scat(lambda, icell) * Jnu(:,icell)
         !accumulate b-f
         call opacity_atom_bf_loc(icell, n_lambda_cont, tab_lambda_cont, chic, Snuc)
      else
         chic(:) = chi_cont(:,icell)
         snuc(:) = eta_cont(:,icell)
      endif

      !linear interpolation
      i0 = 2
      do la=1, N
         loop_i : do lac=i0, n_lambda_cont
            if (tab_lambda_cont(lac) > lambda(la)) then
               w = (lambda(la) - tab_lambda_cont(lac-1)) / (tab_lambda_cont(lac) - tab_lambda_cont(lac-1))
               chi(la) = (1.0_dp - w) * chic(lac-1)  + w * chic(lac)
               snu(la) = (1.0_dp - w) * snuc(lac-1)  + w * snuc(lac)
               i0 = lac
               exit loop_i
            endif
         enddo loop_i
      enddo
      chi(N) = chic(n_lambda_cont)
      snu(N) = snuc(n_lambda_cont)

      return
   end subroutine contopac_atom_loc
   subroutine contopac_atom(icell,i0,ilam,chic,snuc,chi,snu)
      integer, intent(in) :: icell, ilam
      integer, intent(inout)  :: i0
      real(kind=dp), intent(inout) :: chi, Snu
      real(kind=dp), dimension(N_lambda_cont), intent(in) :: chic, snuc
      integer :: la, lac, N
      real(kind=dp) :: w, lam

      !linear interpolation
      lam = tab_lambda_nm(ilam)
      loop_i : do lac=i0, n_lambda_cont
         if (tab_lambda_cont(lac) > lam) then
            w = (lam - tab_lambda_cont(lac-1)) / (tab_lambda_cont(lac) - tab_lambda_cont(lac-1))
            chi = (1.0_dp - w) * chic(lac-1)  + w * chic(lac)
            snu = (1.0_dp - w) * snuc(lac-1)  + w * snuc(lac)
            i0 = lac
            return
            exit loop_i
         endif
      enddo loop_i

      if (ilam==n_lambda) then
         chi = chic(n_lambda_cont)
         snu = snuc(n_lambda_cont)
      endif

      return
   end subroutine contopac_atom

   subroutine opacity_atom_bf_loc(icell,N,lambda,chi,Snu)
   !to do: remove lambda dep since it must be consistent with Nr, Nb
      integer, intent(in) :: icell, N
      real(kind=dp), intent(in), dimension(N) :: lambda
      real(kind=dp), intent(inout), dimension(N) :: chi, Snu
      integer :: nat, Nred, Nblue, kr, i, j
      type(AtomType), pointer :: atom
      real(kind=dp) :: gij, ni_on_nj_star
      real(kind=dp), dimension(N) :: ehnukt
      ! real(kind=dp), dimension(Nlambda_max_cont) :: dissolve

      ehnukt(:) = exp(hnu_k/T(icell))!exphckT(icell)**(lambda_base/lambda(:))

      atom_loop : do nat = 1, N_atoms
         atom => Atoms(nat)%p

         tr_loop : do kr = 1,atom%Ncont

            if (.not.atom%continua(kr)%lcontrib) cycle

            !beware Nc here, assumed tab_lambda_cont
            Nred = atom%continua(kr)%Nrc; Nblue = atom%continua(kr)%Nbc
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            !ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
            ni_on_nj_star = atom%nstar(i,icell)/(atom%nstar(j,icell) + 1d-100)

            gij = ni_on_nj_star * exp(-hc_k/T(icell)/atom%continua(kr)%lambda0)
            if ((atom%n(i,icell) - atom%n(j,icell) * gij) <= 0.0_dp) then
               cycle tr_loop
            endif

            ! 1 if .not. ldissolve (lambdamax==lambda0)
            ! dissolve(1:Nred-Nblue+1) = f_dissolve(T(icell), ne(icell), hydrogen%n(1,icell), atom%continua(kr), Nred-Nblue+1, lambda(Nblue:Nred))

            !should be the same as directly exp(-hc_kT/lambda)
            chi(Nblue:Nred) = chi(Nblue:Nred) + atom%continua(kr)%alpha(:) * (atom%n(i,icell) - &
               ni_on_nj_star * ehnukt(Nblue:Nred) * atom%n(j,icell))! * dissolve(:)

            Snu(Nblue:Nred) = Snu(Nblue:Nred) + atom%n(j,icell) * atom%continua(kr)%alpha(:) * atom%continua(kr)%twohnu3_c2(:) *&
               ni_on_nj_star * ehnukt(Nblue:Nred)! * dissolve(:)

         end do tr_loop

         atom => NULL()
      end do atom_loop


      return
   end subroutine opacity_atom_bf_loc

   function calc_vloc(icell,u,v,w,x,y,z,x1,y1,z1)
      !computes the local mean velocity in direction (u,v,w) of the cell icell
      !by averaging between points (x,y,z) and (x1,y1,z1).
      real(kind=dp) :: calc_vloc
      integer, intent(in) :: icell
      real(kind=dp), intent(in) :: x, y, z, x1, y1, z1
      real(kind=dp), intent(in) :: u, v, w

      calc_vloc = 0.5 * (v_proj(icell,x1,y1,z1,u,v,w) + v_proj(icell,x,y,z,u,v,w))
      return
   end function calc_vloc

   subroutine opacity_atom_bb_loc(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l_void_before,l_contrib, &
         iterate,N,lambda,chi,Snu)
   !to do: remove lambda dep since it must be consistent with Nr, Nb
      integer, intent(in) :: id, icell, iray,N
      logical, intent(in) :: iterate
      real(kind=dp), intent(in), dimension(N) :: lambda
      real(kind=dp), intent(inout), dimension(N) :: chi, Snu
      real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l_void_before,l_contrib
      integer :: nat, Nred, Nblue, kr, i, j, Nlam
      real(kind=dp) :: dv, vg_max, vth, peak_g
      type(AtomType), pointer :: atom
      real(kind=dp), dimension(Nlambda_max_line) :: phi0, vline
      integer :: nvel_rt
      real(kind=dp), dimension(N_gauss) :: v_gauss, phi_gauss


      dv = 0.0_dp
      if (lnon_lte_loop.and..not.iterate) then !not iterate but non-LTE
         dv = calc_vloc(icell,u,v,w,x,y,z,x1,y1,z1) - vlabs(iray,id)
      endif

      atom_loop : do nat = 1, N_Atoms
         atom => Atoms(nat)%p

         vth = vbroad(T(icell),Atom%weight, vturb(icell))

         !any line of that atom with a Gaussian profile ? 
         !use the same profile for all
         !TO DO: before each mode (nlte, image/flux)
         if (atom%lany_gauss_prof) then
            if (lnon_lte_loop) then               
            !completely linear for Gaussian with the sum of Ngauss_log and Ngauss_lin (see line_lambda_grid in gas/wavelengths_gas.f90)
               nvel_rt = N_gauss
               peak_g = peak_gauss_limit / sqrt(pi) / vth
               vg_max = vth * sqrt(-log(peak_g * sqrt(pi) * vth))
            else 
            !image mode or flux, linear grid in hv
               nvel_rt = atom%n_speed_rt
               vg_max = 1d3 * atom%vmax_rt
            endif
            v_gauss(1:nvel_rt) = span_dp(-vg_max,vg_max,nvel_rt,1)
            phi_gauss(1:nvel_rt) = gauss_profile(id,icell,iray,iterate,nvel_rt,v_gauss,vth,&
                                    x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
         endif

         tr_loop : do kr = 1,atom%Nline

            if (.not.atom%lines(kr)%lcontrib) cycle

            !if .not.labs (image or LTE), Nred and Nblue includes the maxium extension of line
            !due to velocity fields (and not only the natural width + damping)
            Nred = atom%lines(kr)%Nr; Nblue = atom%lines(kr)%Nb
            Nlam = atom%lines(kr)%Nlambda
            i = atom%lines(kr)%i; j = atom%lines(kr)%j

            if ((atom%n(i,icell) - atom%n(j,icell)*atom%lines(kr)%gij) <= 0.0_dp) cycle tr_loop

            !Expand the edge of a profile during the non-LTE loop if necessary.
            !the condition is only possibly reached if dv > 0 (lnon_lte_loop = .true.)
            if (lnon_lte_loop) then
               if (abs(dv)>1.0*vbroad(T(icell),Atom%weight, vturb(icell))) then
                  Nred = atom%lines(kr)%Nover_sup
                  Nblue = atom%lines(kr)%Nover_inf
                  Nlam = Nred - Nblue + 1
               endif
            endif

            if (atom%lines(kr)%voigt) then
               phi0(1:Nlam) = profile_art(atom%lines(kr),id,icell,iray,iterate,Nlam,lambda(Nblue:Nred),&
                                 x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
            else
               if (Nlam==nvel_rt) then
                  phi0(1:Nlam) = phi_gauss(1:Nlam)
               else
                  !TO DO: compare the local interp with master and clean. The profile is already shifted
                  ! we just want to interpolate it on another grid. Check that we don't apply the shift twice.
                  vline(1:Nlam) = c_light * (lambda(Nblue:Nred) - atom%lines(kr)%lambda0)/atom%lines(kr)%lambda0
                  phi0(1:Nlam) = linear_1D_sorted(nvel_rt,v_gauss,phi_gauss(1:nvel_rt),Nlam,vline(1:Nlam))
               endif
            endif

            chi(Nblue:Nred) = chi(Nblue:Nred) + &
               hc_fourPI * atom%lines(kr)%Bij * phi0(1:Nlam) * (atom%n(i,icell) - atom%lines(kr)%gij*atom%n(j,icell))

            Snu(Nblue:Nred) = Snu(Nblue:Nred) + &
               hc_fourPI * atom%lines(kr)%Aji * phi0(1:Nlam) * atom%n(j,icell)

            if ((iterate.and.atom%active)) then
               phi_loc(1:Nlam,atom%ij_to_trans(i,j),atom%activeindex,iray,id) = phi0(1:Nlam)
            endif


            !if (lmagnetized) then
            !endif


         end do tr_loop

         atom => null()

      end do atom_loop


      return
   end subroutine opacity_atom_bb_loc
   subroutine lineopac_atom(id, icell, iray, iterate, ilam, chi, Snu)

      integer, intent(in) :: id, icell, iray, ilam
      logical, intent(in) :: iterate
      real(kind=dp), intent(inout) :: chi, Snu
      integer :: nat, Nred, Nblue, kr, i, j, Nlam, la
      real(kind=dp) :: dv, lam, phi0(1)
      type(AtomicLine), pointer :: line

      lam = tab_lambda_nm(ilam)

   !-> this should be outside + need common gauss prof outside
      ! dv = 0.0_dp
      ! if (lnon_lte_loop.and..not.iterate) then !not iterate but non-LTE
      !    !the mean of omegav is essentially the same as vloc (2points).
      !    ! dv = sum(omegav(1:Nvel_points(id),id))/real(Nvel_points(id))
      !    dv = 0.5*(omegav(1,id)+omegav(Nvel_points(id),id))
      !   ! the - vlabs is already accounted for here.
      ! endif

      !loop over lines at that frequency
      do kr = 1, art_lines(ilam)%Ntrans
         line => art_lines(ilam)%art(kr)%p

         Nred = line%Nr; Nblue = line%Nb
         Nlam = line%Nlambda
         i = line%i; j = line%j

         !move outside ?
         ! if ((line%atom%n(i,icell) - line%atom%n(j,icell)*line%gij) <= 0.0_dp) cycle


         ! if (lnon_lte_loop) then
         !    if (abs(dv)>1.0*vbroad(T(icell),line%Atom%weight, vturb(icell))) then
         !       Nred = line%Nover_sup
         !       Nblue = line%Nover_inf
         !       Nlam = Nred - Nblue + 1
         !    endif
         ! endif

         phi0 = 1.0/sqrtpi/line%atom%vth(icell)!profile_art_mono(line,id,icell,lam)
! TO DO: store hc_fourPI * line%Bij * (line%atom%n(i,icell) - line%gij*line%atom%n(j,icell))
! and hc_fourPI * line%Aji * line%atom%n(j,icell) before the loop over frequency
         ! chi = chi + &
         !    hc_fourPI * line%Bij * phi0(1) * (line%atom%n(i,icell) - line%gij*line%atom%n(j,icell))

         ! Snu = Snu + &
         !    hc_fourPI * line%Aji * phi0(1) * line%atom%n(j,icell)
         chi = chi + phi0(1) * line%atom%n(i,icell)

         Snu = Snu + phi0(1) * line%atom%n(j,icell)


         if ((iterate.and.line%atom%active)) then
            !index on Nlam
            la = ilam - Nblue + 1 !1 if Nblue == la
                                  !when la == Nred -> Nlam
            phi_loc(la,line%atom%ij_to_trans(i,j),line%atom%activeindex,iray,id) = phi0(1)
         endif


         line => null()

      end do


      return
   end subroutine lineopac_atom


   subroutine xcoupling(id, icell, iray)
      integer, intent(in) :: id, icell, iray
      integer :: nact, j, i, kr, Nb, Nr, la, Nl
      integer :: i0, j0, la0, N1, N2
      type (AtomType), pointer :: aatom
      real(kind=dp) :: gij, chicc, wl, ni_on_nj_star, wphi, twohnu3_c2, wt
      real(kind=dp) :: term1(Nlambda_max_cont), term2(Nlambda_max_cont), term3(Nlambda_max_cont)
      real(kind=dp), dimension(Nlambda_max_line) :: phi0, wei_line

      Uji_down(:,:,:,id) = 0.0_dp
      chi_down(:,:,:,id) = 0.0_dp
      chi_up(:,:,:,id)   = 0.0_dp

      eta_atoms(:,:,id) = 0.0_dp

      aatom_loop : do nact=1, Nactiveatoms

         aatom => ActiveAtoms(nact)%p

         cont_loop : do kr = 1, aatom%Ncont

            if (.not.aatom%continua(kr)%lcontrib) cycle cont_loop

            j = aatom%continua(kr)%j
            i = aatom%continua(kr)%i
            Nb = aatom%continua(kr)%Nbc; Nr = aatom%continua(kr)%Nrc
            Nl = Nr-Nb+1
            N1 = aatom%continua(kr)%Nb; N2 = aatom%continua(kr)%Nr

            !ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
            ni_on_nj_star = aatom%nstar(i,icell)/(aatom%nstar(j,icell) + 1d-100)


            gij = ni_on_nj_star * exp(-hc_k/T(icell)/aatom%continua(kr)%lambda0)


            if (aatom%n(i,icell) - gij*aatom%n(j,icell) <= 0.0_dp) cycle cont_loop

            !The wavelength integration weight cannot be used in that loop because we interpolate after.
            freq_loop : do la=1, Nl
               gij = ni_on_nj_star * exp(hnu_k(Nb-1+la)/T(icell))
               term1(la) = fourpi_h * aatom%continua(kr)%alpha(la) * (aatom%n(i,icell) - gij*aatom%n(j,icell))
               term2(la) = aatom%continua(kr)%alpha(la) * aatom%continua(kr)%twohnu3_c2(la) * gij
               term3(la) = aatom%continua(kr)%alpha(la) * aatom%continua(kr)%twohnu3_c2(la) * gij * aatom%n(j,icell)
            enddo freq_loop
            !linear interpolation + adding wavelength integration weight
            !wl = 0.5*(tab_lambda_nm(N1+1)-tab_lambda_nm(N1)) / tab_lambda_nm(N1)
            ! chi_down(N1,j,nact,id) = chi_down(N1,j,nact,id) + term1(1)! * wl
            ! chi_up(n1,i,nact,id) = chi_up(n1,i,nact,id) + term1(1)! * wl
            ! Uji_down(n1,j,nact,id) = Uji_down(N1,j,nact,id) + term2(1)
            ! eta_atoms(N1,nact,id) = eta_atoms(N1,nact,id) + term3(1)
            i0 = 2
            do la=N1,N2
               ! if (la>1) then
               !    wl = 0.5*(tab_lambda_nm(la+1)-tab_lambda_nm(la-1)) / tab_lambda_nm(la)
               ! !otherwise, wl is the previous one, first point of the grid
               ! endif
               loop_i : do la0=i0, Nl
                  if (tab_lambda_cont(Nb+la0-1) > tab_lambda_nm(la)) then
                     wt = (tab_lambda_nm(la) - tab_lambda_cont(Nb+la0-2)) / (tab_lambda_cont(Nb+la0-1) - tab_lambda_cont(Nb+la0-2))

                     chi_down(la,j,nact,id) = chi_down(la,j,nact,id) + ( (1.0_dp - wt) *term1(la0-1)  + wt * term1(la0) )!*wl
                     chi_up(la,i,nact,id) = chi_up(la,i,nact,id) + ( (1.0_dp - wt) *term1(la0-1)  + wt * term1(la0) )!*wl
                     Uji_down(la,j,nact,id) = Uji_down(la,j,nact,id) + (1.0_dp - wt) *term2(la0-1)  + wt * term2(la0)
                     eta_atoms(la,nact,id) = eta_atoms(la,nact,id) + (1.0_dp - wt) *term3(la0-1)  + wt * term3(la0)

                     i0 = la0
                     exit loop_i
                  endif
               enddo loop_i
            enddo
            ! wl = 0.5*(tab_lambda_nm(N2)-tab_lambda_nm(N2-1)) / tab_lambda_nm(N2-1)
            chi_down(n2,j,nact,id) = chi_down(n2,j,nact,id) + term1(Nl)! * wl
            chi_up(n2,i,nact,id) = chi_up(n2,i,nact,id) + term1(Nl)! * wl
            Uji_down(n2,j,nact,id) = Uji_down(n2,j,nact,id) + term2(Nl)
            eta_atoms(n2,nact,id) = eta_atoms(n2,nact,id) + term3(Nl)

         enddo cont_loop

         do la=1, n_lambda
             if (la==1) then
                wl = 0.5*(tab_lambda_nm(la+1)-tab_lambda_nm(la)) / tab_lambda_nm(la)
             elseif (la==n_lambda) then
                wl = 0.5*(tab_lambda_nm(la)-tab_lambda_nm(la-1)) / tab_lambda_nm(la-1)
             else
                wl = 0.5*(tab_lambda_nm(la+1)-tab_lambda_nm(la-1)) / tab_lambda_nm(la)
             endif
             chi_down(la,:,nact,id) =  chi_down(la,:,nact,id)  * wl
             chi_up(la,:,nact,id) = chi_up(la,:,nact,id) * wl
         enddo

         line_loop : do kr = 1, aatom%Nline

            if (.not.aatom%lines(kr)%lcontrib) cycle line_loop

            j = aatom%lines(kr)%j
            i = aatom%lines(kr)%i

            if (aatom%n(i,icell) - aatom%lines(kr)%gij*aatom%n(j,icell) <= 0.0_dp) cycle line_loop

            Nb = aatom%lines(kr)%Nb; Nr = aatom%lines(kr)%Nr
            Nl = Nr - Nb + 1

            wphi = 0.0
            do la=1,Nl
               if (la==1) then
                  wl = 0.5*(tab_lambda_nm(Nb+1)-tab_lambda_nm(Nb)) * c_light / aatom%lines(kr)%lambda0
               elseif (la==Nl) then
                  wl = 0.5*(tab_lambda_nm(Nr)-tab_lambda_nm(Nr-1)) * c_light / aatom%lines(kr)%lambda0
               else
                  wl = 0.5*(tab_lambda_nm(Nb+la)-tab_lambda_nm(Nb+la-2)) * c_light / aatom%lines(kr)%lambda0
               endif
               wei_line(la) = wl
               phi0(la) = phi_loc(la,kr,nact,iray,id)
               wphi = wphi + wl * phi0(la)
            enddo

            freq2_loop : do la=1, Nl
               wl = wei_line(la)

               Uji_down(Nb-1+la,j,nact,id) = Uji_down(Nb-1+la,j,nact,id) + hc_fourPI * aatom%lines(kr)%Aji * phi0(la)!/wphi

               !there is a hc_fourpi factor that simplifies here, because integral is over dnu/hnu dOmega = dv/hc dOmega * hc/4pi
               !dv dOmega/4pi which is whtat is contained in wl (for dv) and what the angular integration provides (dOmega/4pi)
               chicc = wl * aatom%lines(kr)%Bij * (aatom%n(i,icell) - aatom%lines(kr)%gij*aatom%n(j,icell)) * phi0(la)/wphi

               chi_down(Nb-1+la,j,nact,id) = chi_down(Nb-1+la,j,nact,id) + chicc
               chi_up(Nb-1+la,i,nact,id) = chi_up(Nb-1+la,i,nact,id) + chicc

               eta_atoms(Nb-1+la,nact,id) = eta_atoms(Nb-1+la,nact,id) + &
                  hc_fourPI * aatom%lines(kr)%Aji * aatom%n(j,icell) * phi0(la)!/wphi


            enddo freq2_loop

         enddo line_loop

         aatom => null()
      enddo aatom_loop

    return
   end subroutine xcoupling

   function profile_art(line,id,icell,iray,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
   use voigts, only : VoigtThomson
      ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: id, icell, N, iray
      type (AtomicLine), intent(in)          :: line
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer 											:: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm, vth
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: u0, profile_art, u1, u0sq

      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nr; Nblue = line%Nb
      vth = vbroad(T(icell),line%Atom%weight, vturb(icell))

      u0(:) = (lambda - line%lambda0)/line%lambda0  * ( c_light/vth )

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
         Nvspace = 1
         omegav_mean = v0
      else

         Omegav = 0.0
         omegav(1) = v0
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/vth*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
         omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)
      endif
      !in non-LTE:
      !the actual cell icell_nlte must be centered on 0 (moving at vmean).
      !the other cells icell crossed must be centered in v(icell) - vmean(icell_nlte)
      if (lsubstract_avg) then!labs == .true.
         omegav(1:Nvspace) = omegav(1:Nvspace) - omegav_mean
         vlabs(iray,id) = omegav_mean
      else
         if (lnon_lte_loop) omegav(1:Nvspace) = omegav(1:Nvspace) - vlabs(iray,id)
      endif


      if (line%voigt) then
         u1(:) = u0(:) - omegav(1)/vth
         ! profile_art(:) = Voigt(N, line%a(icell), u1(:))
         ! do nv=2, Nvspace
         !    u1(:) = u0(:) - omegav(nv)/vth
         !    profile_art(:) = profile_art(:) + Voigt(N, line%a(icell), u1(:))
         ! enddo
         if (lnon_lte_loop) then!approximate for non-LTE
            profile_art(:) = VoigtThomson(N,line%a(icell), u1(:),vth)
            do nv=2, Nvspace
               u1(:) = u0(:) - omegav(nv)/vth
               profile_art(:) = profile_art + VoigtThomson(N,line%a(icell), u1(:),vth)
            enddo
         else!accurate for images
            profile_art(:) = Voigt(N, line%a(icell), u1(:))
            do nv=2, Nvspace
               u1(:) = u0(:) - omegav(nv)/vth
               profile_art(:) = profile_art(:) + Voigt(N, line%a(icell), u1(:))
            enddo
         endif

      else
         u0sq(:) = u0(:)*u0(:)
         !Note: u1 = (u0 - omegav(nv)/vth)**2
         u1(:) = u0sq(:) + (omegav(1)/vth)*(omegav(1)/vth) - 2*u0(:) * omegav(1)/vth
         profile_art(:) = exp(-u1(:))
         do nv=2, Nvspace
            u1(:) = u0sq(:) + (omegav(nv)/vth)*(omegav(nv)/vth) - 2*u0(:) * omegav(nv)/vth
            profile_art(:) = profile_art(:) + exp(-u1(:))
         enddo
      endif

      norm = Nvspace * vth * sqrtpi
      profile_art(:) = profile_art(:) / norm

      return
   end function profile_art
   function profile_art_mono(line,id,icell,lam)
   use voigts, only : VoigtThomson_b
      ! phi = Voigt / sqrt(pi) / vbroad(icell)
!TO DO: move omegav outside for each atom
      integer, intent(in)                    :: id, icell
      type (AtomicLine), intent(in)          :: line
      real(kind=dp), intent(in)              :: lam
      integer 											:: Nvspace
      real(kind=dp)                          :: norm, vth
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp)                          :: u0, profile_art_mono(1), u1(1), u0sq

   
      i = line%i; j = line%j
      vth = line%atom%vth(icell)
      Nred = line%Nr; Nblue = line%Nb

      Nvspace = Nvel_points(id)

      u0 = (lam - line%lambda0)/line%lambda0  * ( c_light/vth )

      if (line%voigt) then
         u1(1) = u0 - omegav(1,id)/vth
         if (lnon_lte_loop) then!approximate for non-LTE
            profile_art_mono(1) = VoigtThomson_b(line%a(icell),line%b(icell),line%c(icell), u1(1),vth)
            do nv=2, Nvspace
               u1(1) = u0 - omegav(nv,id)/vth
               profile_art_mono(1) = profile_art_mono(1) + VoigtThomson_b(line%a(icell),line%b(icell),line%c(icell), u1(1),vth)
            enddo
         else!accurate for images
            profile_art_mono = Voigt(1, line%a(icell), u1(1))
            do nv=2, Nvspace
               u1(1) = u0 - omegav(nv,id)/vth
               profile_art_mono = profile_art_mono + Voigt(1, line%a(icell), u1(1))
            enddo
         endif

      else
         u0sq = u0*u0
         !Note: u1 = (u0 - omegav(nv)/vth)**2
         u1(1) = u0sq + (omegav(1,id)/vth)*(omegav(1,id)/vth) - 2*u0 * omegav(1,id)/vth
         profile_art_mono(1) = exp(-u1(1))
         do nv=2, Nvspace
            u1(1) = u0sq + (omegav(nv,id)/vth)*(omegav(nv,id)/vth) - 2*u0 * omegav(nv,id)/vth
            profile_art_mono(1) = profile_art_mono(1) + exp(-u1(1))
         enddo
      endif

      norm = Nvspace * vth * sqrtpi
      profile_art_mono = profile_art_mono / norm

      return
   end function profile_art_mono
   subroutine calc_omegav (id,icell,iray,lsubstract_avg,x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      integer, intent(in)                    :: id, icell,iray
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer 											:: Nvspace
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean, vth
      integer                                :: nv, nat
      type(AtomType), pointer                :: at

      Nvspace = NvspaceMax
      v0 = v_proj(icell,x,y,z,u,v,w)
      Nvel_points(id) = Nvspace !hmmm tmp

      if (lvoronoi) then
         omegav(1,id) = v0
         Nvspace = 1
         omegav_mean = v0
         if (lsubstract_avg) then!labs == .true.
            omegav(1:Nvspace,id) = omegav(1:Nvspace,id) - omegav_mean
            vlabs(iray,id) = omegav_mean
         else
            if (lnon_lte_loop) omegav(1:Nvspace,id) = omegav(1:Nvspace,id) - vlabs(iray,id)
         endif
         return
      endif

!TO DO
      !common Nvspace for all atoms at the moment
      vth = 1d50
      do nat=1, n_atoms
         at => atoms(nat)%p
         vth = min(vth,vbroad(T(icell),at%weight, vturb(icell)))
      enddo
      at => null()


      Omegav(:,id) = 0.0
      omegav(1,id) = v0
      v1 = v_proj(icell,x1,y1,z1,u,v,w)
      dv = abs(v1-v0)
      Nvspace = min(max(2,nint(dv/vth*20.)),NvspaceMax)
      Nvel_points(id) = Nvspace !hmmmm tmp

      do nv=2, Nvspace-1
         delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
         xphi=x+delta_vol_phi*u
         yphi=y+delta_vol_phi*v
         zphi=z+delta_vol_phi*w
         omegav(nv,id) = v_proj(icell,xphi,yphi,zphi,u,v,w)
      enddo
      omegav(Nvspace,id) = v1
      omegav_mean = sum(omegav(1:Nvspace,id))/real(Nvspace,kind=dp)

      if (lsubstract_avg) then!labs == .true.
         omegav(1:Nvspace,id) = omegav(1:Nvspace,id) - omegav_mean
         vlabs(iray,id) = omegav_mean
      else
         if (lnon_lte_loop) omegav(1:Nvspace,id) = omegav(1:Nvspace,id) - vlabs(iray,id)
      endif
   end subroutine calc_omegav

   function profile_art_i(line,id,icell,iray,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      integer, intent(in)                    :: icell, N, id, iray
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer 											:: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm, vth
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      type (AtomicLine), intent(in)          :: line
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: uloc, u0, profile_art_i, u1, u0sq

      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nr; Nblue = line%Nb
      vth = vbroad(T(icell),line%Atom%weight, vturb(icell))

      u0(:) = (lambda - line%lambda0)/line%lambda0  * ( c_light/vth )

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
         Nvspace = 1
         omegav_mean = v0
      else

         Omegav = 0.0
         omegav(1) = v0
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/vth*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
         omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)
      endif
      if (lsubstract_avg) then
         omegav(1:Nvspace) = omegav(1:Nvspace) - omegav_mean
         vlabs(iray,id) = omegav_mean
      else
         if (lnon_lte_loop) omegav(1:Nvspace) = omegav(1:Nvspace) - vlabs(iray,id)
      endif

      if (line%voigt) then
         u1(:) = u0(:) - omegav(1)/vth
         uloc(:) = line%v(:) / vth
         profile_art_i(:) = linear_1D_sorted(line%nlambda,uloc(:),line%phi(:,icell),N,u1(:))
         do nv=2, Nvspace
            u1(:) = u0(:) - omegav(nv)/vth
            profile_art_i(:) = profile_art_i(:) + linear_1D_sorted(line%nlambda,uloc(:),line%phi(:,icell),N,u1(:))
         enddo

      else
         u0sq(:) = u0(:)*u0(:)
         u1(:) = u0sq(:) + (omegav(1)/vth)*(omegav(1)/vth) - 2*u0(:) * omegav(1)/vth
         profile_art_i(:) = exp(-u1(:)) / sqrtpi / vth
         do nv=2, Nvspace
            u1(:) = u0sq(:) + (omegav(nv)/vth)*(omegav(nv)/vth) - 2*u0(:) * omegav(nv)/vth
            profile_art_i(:) = profile_art_i(:) + exp(-u1(:)) / sqrtpi / vth
         enddo
      endif

      profile_art_i(:) = profile_art_i(:) / Nvspace

      return
   end function profile_art_i

   function gauss_profile(id,icell,iray,lsubstract_avg,N,vel,vth,x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
   use voigts, only : VoigtThomson
      ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: id,icell, iray,N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: vel
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib, & !physical length of the cell
                                             vth
      integer 											:: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      integer                                :: nv
      real(kind=dp), dimension(N)            :: u0, gauss_profile, u1, u0sq

      Nvspace = NvspaceMax

      u0(:) = vel/vth

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
         Nvspace = 1
         omegav_mean = v0
      else

         Omegav = 0.0
         omegav(1) = v0
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/vth*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
         omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)
      endif
      !in non-LTE:
      !the actual cell icell_nlte must be centered on 0 (moving at vmean).
      !the other cells icell crossed must be centered in v(icell) - vmean(icell_nlte)
      if (lsubstract_avg) then!labs == .true.
         omegav(1:Nvspace) = omegav(1:Nvspace) - omegav_mean
         vlabs(iray,id) = omegav_mean
      else
         if (lnon_lte_loop) omegav(1:Nvspace) = omegav(1:Nvspace) - vlabs(iray,id)
      endif

   
      u0sq(:) = u0(:)*u0(:)
      !Note: u1 = (u0 - omegav(nv)/vth)**2
      u1(:) = u0sq(:) + (omegav(1)/vth)*(omegav(1)/vth) - 2*u0(:) * omegav(1)/vth
      gauss_profile(:) = exp(-u1(:))
      do nv=2, Nvspace
         u1(:) = u0sq(:) + (omegav(nv)/vth)*(omegav(nv)/vth) - 2*u0(:) * omegav(nv)/vth
         gauss_profile(:) = gauss_profile(:) + exp(-u1(:))
      enddo

      norm = Nvspace * vth * sqrtpi
      gauss_profile(:) = gauss_profile(:) / norm

      return
   end function gauss_profile

   subroutine write_opacity_emissivity_bin(Nlambda,lambda)
   !To do: store opacity in 3d arrays in space instead of icell
      integer, intent(in) :: Nlambda
      real(kind=dp), intent(in) :: lambda(Nlambda)
      integer :: unit, unit2, status = 0
      integer :: alloc_status, id, icell, m, Nrec
      real(kind=dp), allocatable, dimension(:,:,:) :: chi_tmp, eta_tmp, rho_tmp
      real(kind=dp), allocatable, dimension(:,:,:) :: chic_tmp, etac_tmp
      character(len=11) :: filename_chi="chi.b"
      character(len=50) :: filename_eta="eta.b"
      character(len=18) :: filename_rho="magnetoopt.b"
      real(kind=dp) :: zero_dp, u, v, w

      zero_dp = 0.0_dp
      u = zero_dp; v = zero_dp; w = zero_dp
      write(*,*) " ** Writing emissivity and opacity (zero-velocity)..."
      ! if (lmagnetized) then
      !    Nrec = 4
      !    allocate(rho_tmp(Nlambda,n_cells,Nrec-1),stat=alloc_status)
      !    if (alloc_status /= 0) call error("Cannot allocate rho_tmp !")
      !    call warning("  Polarized case with an observer // to z!")
      ! else
         Nrec = 1
      ! endif
      allocate(chi_tmp(Nlambda, n_cells, Nrec),stat=alloc_status)
      if (alloc_status /= 0) call error("(write_opac_bin) Cannot allocate chi_tmp !")
      allocate(eta_tmp(Nlambda, n_cells, Nrec),stat=alloc_status)
      if (alloc_status /= 0) call error("(write_opac_bin) Cannot allocate eta_tmp !")
      allocate(chic_tmp(Nlambda, n_cells, 1),stat=alloc_status)
      if (alloc_status /= 0) call error("(write_opac_bin) Cannot allocate chic_tmp !")
      allocate(etac_tmp(Nlambda, n_cells, 1),stat=alloc_status)
      if (alloc_status /= 0) call error("(write_opac_bin) Cannot allocate etac_tmp !")

      chi_tmp = 0.0; eta_tmp = 0.0; chic_tmp = 0.0; etac_tmp = 0.0

      call ftgiou(unit,status)
      call ftgiou(unit2,status)
      open(unit, file=trim(filename_chi),form="unformatted",status='unknown',access="stream",iostat=status)
      open(unit2, file=trim(filename_eta),form="unformatted",status='unknown',access="stream",iostat=status)
      !write wavelength first
      write(unit,iostat=status) shape(chi_tmp)
      write(unit2,iostat=status) shape(eta_tmp)
      write(unit,iostat=status) lambda
      write(unit2,iostat=status) lambda
      id = 1
      do icell=1, n_cells
         !$ id = omp_get_thread_num() + 1
         if (icompute_atomRT(icell) > 0) then
            call contopac_atom_loc(icell,Nlambda,lambda,chi_tmp(:,icell,1),eta_tmp(:,icell,1))
            chic_tmp(:,icell,1) = chi_tmp(:,icell,1)
            etac_tmp(:,icell,1) = eta_tmp(:,icell,1)
            call opacity_atom_bb_loc(id,icell,1,1d0,zero_dp,zero_dp,1d0,zero_dp,zero_dp,u,v,w,&
               zero_dp,zero_dp,.false.,Nlambda,lambda,chi_tmp(:,icell,1), eta_tmp(:,icell,1))
            ! do m=2,Nrec
            !    !m=1, unpolarized already filled.
            !    !etaQUV, chiQUV only for Q(=1), U, V
            !    eta_tmp(:,m,icell) = etaQUV_p(:,m-1,id)
            !    chi_tmp(:,m,icell) = chiQUV_p(:,m-1,id)
            !    !only QUV for rho_p
            !    rho_tmp(:,m-1,icell) = rho_p(:,m-1,id)
            ! enddo
         endif
      enddo

      write(unit,iostat=status) chi_tmp
      write(unit,iostat=status) chic_tmp
      write(unit2,iostat=status) eta_tmp
      write(unit2,iostat=status) etac_tmp
      ! if (lmagnetized) then
      !   write(unit, iostat=status) rho_tmp
      !   deallocate(rho_tmp)
      ! endif
      deallocate(chi_tmp, eta_tmp, chic_tmp, etac_tmp)
      close(unit)
      close(unit2)

      !free unit
      call ftfiou(unit, status)
      call ftfiou(unit2, status)

      return
   end subroutine write_opacity_emissivity_bin

end module Opacity_atom
