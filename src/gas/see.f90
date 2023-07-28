module see

    use atom_type
    use elements_type
    use grid
    use parametres
    use gas_contopac, only        : H_bf_Xsection
    use wavelengths, only         :  n_lambda
    use wavelengths_gas, only     : Nlambda_max_line, Nlambda_max_trans, Nlambda_max_cont, n_lambda_cont, &
         tab_lambda_cont, tab_lambda_nm, line_weights, cont_weights
    use utils, only               : gaussslv, solve_lin, is_nan_infinity_vector, linear_1D_sorted, is_nan_infinity_matrix, jacobi
    use opacity_atom, only : phi_loc, psi, chi_up, chi_down, uji_down, Itot, eta_atoms, psi_odiag, Ujdown_odiag, cells_id, lambda_ij!, eta_odiag
    use messages, only : warning, error
    use collision_atom, only : collision_rates_atom_loc, collision_rates_hydrogen_loc
    use fits_utils, only : print_error

    implicit none

    !populations below that threshold not taken into account in convergence.
    real(kind=dp), parameter :: frac_limit_pops = 1d-15!1d-50
    !Variables for Non-LTE loop and MALI method
    real(kind=dp), allocatable ::  ngpop(:,:,:,:)
    real(kind=dp), allocatable :: tab_Aji_cont(:,:,:), tab_Vij_cont(:,:,:)


    !variables for non-LTE H (+He) ionisation
    integer :: Neq_ne, Neq_ng
    integer, allocatable, dimension(:,:) :: gs_ion !index of the ground states of each ion
    real(kind=dp), allocatable :: gtot(:,:,:), gr_save(:,:,:,:), dgrdne(:,:,:), dgcdne(:,:,:)
    real(kind=dp), allocatable :: pops_ion(:,:,:), npop_dag(:,:)
    real(kind=dp), allocatable :: fvar(:,:), xvar(:,:), dfvar(:,:,:)!Jacobian

    contains

    function innz(nl,i,l,j,lp)
        ! For a matrix A of (nl*n,(1+m)*nl) with nl levels, n non empty cells and m neighbours
        ! returns the index innz on a (nl**2 * n * (m+1)) vector of element A^{i,j}_{lp,l}
        !    /                                                                                   \
        !   | A^{i,j}_{lp,l} A^{i,j}_{lp+1,l} A^{i,j+1}_{lp,l} A^{i,j+1}_{lp+1,l}                 |
        !   | A^{i,j}_{lp,l+1} A^{i,j}_{lp+1,l+1} A^{i,j+1}_{lp,l+1} A^{i,j+1}_{lp+1,l+1}         |
        ! A | A^{i+1,j}_{lp,l} A^{i+1,j}_{lp+1,l} A^{i+1,j+1}_{lp,l} A^{i+1,j+1}_{lp+1,l}         |
        !   | A^{i+1,j}_{lp,l+1} A^{i+1,j}_{lp+1,l+1} A^{i+1,j+1}_{lp,l+1} A^{i+1,j+1}_{lp+1,l+1} |
        !   \                                                                                    /
        !
        integer :: innz
        integer, intent(in) :: nl
        integer, intent(in) :: l, i, j, lp
        integer :: n, m

        n = n_non_empty_cells
        m = 1 + n_neighbours_max

        innz = nl*n*m*(i-1) + m*nl*(l-1) + m*(j-1) + lp
        return
    end function innz

    subroutine alloc_nlte_var(n_rayons_max,mem)
    !allocate space for non-lte loop (only)
        integer, intent(in) :: n_rayons_max
        integer(kind=8), intent(in), optional :: mem
        integer :: Nmaxlevel,Nmaxline,Nmaxcont,NlevelTotal,Nmaxstage
        integer :: alloc_status
        integer(kind=8) :: mem_alloc_local, mem_alloc_non_local_alo
        integer :: kr, n, icell, l, i, j
        real(kind=dp) :: wl, anu1, e1, b1
        real(kind=dp) :: anu2, e2, b2!higher prec integ
        type(AtomType), pointer :: atom

        mem_alloc_local = 0
        if (present(mem)) then
            mem_alloc_local = mem_alloc_local + mem
        endif
        NmaxLevel = 0 !maximum number of levels
        NlevelTotal = 0 !sum of all levels among all active atoms
        NmaxLine = 0 !maximum number of lines
        Nmaxcont = 0 !maximum number of cont
        NmaxStage = 0 !maximum number of ionisation stages
        !neglect size of rates
        do n=1, NactiveAtoms
            atom => ActiveAtoms(n)%p
            Nmaxlevel = max(Nmaxlevel,atom%Nlevel)
            Nmaxcont = max(NmaxCont,atom%Ncont)
            NmaxLine = max(NmaxLine,atom%Nline)
            NlevelTotal = NlevelTotal + atom%Nlevel
            Nmaxstage = max(NmaxStage,atom%Nstage)
            !the local rate matrix, with diagonal operator
            allocate(atom%Gamma(atom%Nlevel,atom%Nlevel,nb_proc))
            mem_alloc_local = mem_alloc_local + sizeof(atom%Gamma)
            do kr=1, atom%Nline
                allocate(atom%lines(kr)%Rij(nb_proc),atom%lines(kr)%Rji(nb_proc))
                allocate(atom%lines(kr)%Cij(nb_proc),atom%lines(kr)%Cji(nb_proc))
                call line_weights(atom%lines(kr))
                mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kr)%wei)
            enddo
            do kr=1, atom%Ncont
                allocate(atom%continua(kr)%Rij(nb_proc),atom%continua(kr)%Rji(nb_proc))
                allocate(atom%continua(kr)%Cij(nb_proc),atom%continua(kr)%Cji(nb_proc))
                call cont_weights(atom%continua(kr),atom%continua(kr)%Nb,atom%continua(kr)%Nr,n_lambda,tab_lambda_nm)
                mem_alloc_local = mem_alloc_local + sizeof(atom%continua(kr)%wei)
            enddo
            atom => null()
        enddo

        !allocate space for the non-local ALO.
        !Currently, the local variabls, %Gamma(:,:,id), etc are kept because
        !electronic density is solved locally with the local populations only.
        !
        call count_neighbours(n_neighbours_max,alo_order) !2, 8, 26 for spherical grid 1D, 2D, 3D
		n_non_empty_cells = size(pack(icompute_atomRT,mask=icompute_atomrt>0))
		allocate(tab_index_cell(n_non_empty_cells),tab_index_i(n_cells))
		allocate(tab_index_neighb(n_non_empty_cells,n_neighbours_max))
        write(*,*) " --- The number of maximum neighbours for each cell is ", n_neighbours_max
        write(*,*) " --- There are ", n_non_empty_cells, " non-empty cells"
        i = 0
        do icell=1,n_cells
            tab_index_i(icell) = 0
            if (icompute_atomrt(icell)>0) then
                i = i + 1
                tab_index_cell(i) = icell
                tab_index_i(icell) = i
                ! 0 means there are less neighbours for that cell than n_neighbours_max.
                tab_index_neighb(i,:) = index_neighbours(n_neighbours_max,icell,alo_order)
            endif
        enddo
        mem_alloc_non_local_alo = sizeof(tab_index_cell)+sizeof(tab_index_neighb)+sizeof(tab_index_i)
        allocate(lambda_ij(n_cells,n_cells),stat=alloc_status)
        write(*,*) " *** Allocating Lambda*_ij to check the sampling of neighbours in the non-local approximate op."
        Lambda_ij = 0.0; mem_alloc_non_local_alo = mem_alloc_non_local_alo  + sizeof(Lambda_ij)
        do n=1, NactiveAtoms
            atom => ActiveAtoms(n)%p
            !Built locally, frequency-angle integrated, for all neighbours of a cell.
            !non-zero elements of the (Nlevel*n_cells,Nlevel*(N_neighbour + 1)) rate matrix
            !NOTE: there are some 0 elements if the number of neighbour of a given cell is < n_neighbours_max.
            allocate(atom%g_nnz(atom%Nlevel,atom%Nlevel, n_non_empty_cells, (1 + n_neighbours_max)))
            if (alloc_Status > 0) then
                write(*,*) 8 * (n_non_empty_cells*(1+n_neighbours_max))*atom%Nlevel**2 / 1024.**3, 'GB'
                call error("Allocation error atom%w")
            endif
            mem_alloc_non_local_alo = mem_alloc_non_local_alo + sizeof(atom%G_nnz)
        enddo
        atom => null()
        if (alo_order > 0) Then
            allocate(psi_odiag(n_lambda, n_neighbours_max, nb_proc), stat=alloc_status); psi_odiag(:,:,:) = 0.0_dp
            write(*,*) " size psi_odiag:", sizeof(psi_odiag) / 1024./1024./1024.," GB"
            mem_alloc_non_local_alo = mem_alloc_non_local_alo + sizeof(psi_odiag)
            if (alloc_Status > 0) call error("Allocation error psi_odiag in nlte loop")
            ! allocate(eta_odiag(n_lambda,NactiveAtoms,n_neighbours_max, nb_proc),stat=alloc_status)
            ! mem_alloc_non_local_alo = mem_alloc_non_local_alo + sizeof(eta_odiag)
            ! if (alloc_Status > 0) call error("Allocation error eta_odiag")
            ! write(*,*) " size eta_atoms:", sizeof(eta_odiag) / 1024./1024./1024.," GB"
            allocate(Ujdown_odiag(n_lambda,Nmaxlevel,NactiveAtoms,n_neighbours_max,nb_proc),stat=alloc_status)
            mem_alloc_non_local_alo = mem_alloc_non_local_alo + sizeof(Ujdown_odiag)
            if (alloc_Status > 0) call error("Allocation error Ujdown_odiag")
            write(*,*) " size Ujdown_odiag:", sizeof(Ujdown_odiag) / 1024./1024./1024.," GB"
            !eta_odiag = sum_i Ujdown_odiag(dim=2) * n
        endif

        !-> allocate space for non-LTE only quantities !
        ! NOTE: index 1 is always the current solution
        !   previous solution is 2, previous previous is 3 etc..
        !   size is NactiveAtoms + 1 for current iterate of electronic density (ne_new)
        Neq_ng = max(Ng_Norder+2,3) !storing 3 by default
        if (lNg_acceleration) then
            !Automatically set if not provided (Ng_Nperiod = -1)
            if (Ng_Nperiod == -1) Ng_Nperiod = Ng_Norder + 2
        endif
        !even if no Ng acceleration allows for storing Neq_ng previous solutions
        !   for the non-LTE populations: TO DO improve convergence criterion.
        allocate(ngpop(NmaxLevel,NactiveAtoms+1,n_cells,Neq_ng),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error ngpop")
        
        !initialize electronic density
        ngpop(1,NactiveAtoms+1,:,1) = ne(:)
        write(*,*) " size ngpop:", sizeof(ngpop)/1024./1024., " MB"
        allocate(psi(n_lambda, n_rayons_max, nb_proc), stat=alloc_status); psi(:,:,:) = 0.0_dp
        write(*,*) " size psi:", sizeof(psi) / 1024./1024.," MB"
        if (alloc_Status > 0) call error("Allocation error psi in nlte loop")
        allocate(eta_atoms(n_lambda,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error eta_atoms")
        write(*,*) " size eta_atoms:", sizeof(eta_atoms) / 1024./1024.," MB"
        call warning("(SEE) eta_atoms used in local rate matrix future deprec (sum Uji * nj)")
        allocate(Uji_down(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error Uji_down")
        write(*,*) " size cross-coupling:", 3 * sizeof(Uji_down) / 1024./1024./1024.," GB"
        allocate(chi_down(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error chi_down")
        allocate(chi_up(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error chi_up")

        allocate(Itot(n_lambda,n_rayons_max,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error Itot")
        allocate(phi_loc(Nlambda_max_line,Nmaxline,NactiveAtoms,1,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error phi_loc")

        write(*,*) " size Itot:", sizeof(Itot) / 1024./1024./1024.," GB"
        write(*,*) " size phi_loc:", sizeof(phi_loc) / 1024./1024./1024.," GB"

        mem_alloc_local = mem_alloc_local + sizeof(psi) + sizeof(eta_atoms) + sizeof(Itot) + sizeof(phi_loc) + &
            sizeof(uji_down)+sizeof(chi_down)+sizeof(chi_up) + sizeof(ngpop)

        allocate(tab_Aji_cont(NmaxCont,NactiveAtoms,n_cells))
        allocate(tab_Vij_cont(Nlambda_max_cont,NmaxCont,NactiveAtoms))
        tab_Aji_cont = 0.0_dp; tab_Vij_cont = 0.0_dp
        write(*,*) " size tab_Aji_cont:", sizeof(tab_Aji_cont) / 1024./1024./1024.," GB"
        write(*,*) " size tab_Vij_cont:", sizeof(tab_Vij_cont) / 1024./1024./1024.," GB"
        !integrate in frequency Uji which is fourpi/h * anu * 2hnu3/c2 * exp(-hnu/kT)
        !for each cont get the cross-sections
        mem_alloc_local = mem_alloc_local + sizeof(tab_Aji_cont) + sizeof(tab_Vij_cont)
        do n=1, NactiveAtoms
            atom => ActiveAtoms(n)%p
            ct_loop : do kr=1,atom%Ncont
                if (.not.atom%continua(kr)%lcontrib) cycle ct_loop
                i = atom%Continua(kr)%i
                j = atom%continua(kr)%j
                if (atom%continua(kr)%hydrogenic) then
                   tab_Vij_cont(1:atom%continua(kr)%Nlambda,kr,n) = H_bf_Xsection(atom%continua(kr), &
                        tab_lambda_nm(atom%continua(kr)%Nb:atom%continua(kr)%Nr))
                else
                    tab_Vij_cont(1:atom%continua(kr)%Nlambda,kr,n) = linear_1D_sorted(size(atom%continua(kr)%alpha_file),&
                         atom%continua(kr)%lambda_file,atom%continua(kr)%alpha_file,atom%continua(kr)%Nlambda,&
                         tab_lambda_nm(atom%continua(kr)%Nb:atom%continua(kr)%Nr))
                    tab_Vij_cont(atom%continua(kr)%Nlambda,kr,n) = atom%continua(kr)%alpha_file(size(atom%continua(kr)%alpha_file))
                endif
                !loop for all cells here can be long
                !could be para
                do icell=1, n_cells
                    if (icompute_atomRT(icell) > 0) then
                        tab_Aji_cont(kr,n,icell) = 0.0_dp !frequency and angle integrated!
                        ! do l=2,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1
                        !     anu1 = tab_Vij_cont(l,kr,n)
                        !     anu2 = tab_Vij_cont(l-1,kr,n)

                        !     e1 = exp(-hc_k/T(icell)/tab_lambda_nm(atom%continua(kr)%Nb+l-1))
                        !     e2 = exp(-hc_k/T(icell)/tab_lambda_nm(atom%continua(kr)%Nb+l-2))

                        !     b1 = twohc/tab_lambda_nm(atom%continua(kr)%Nb+l-1)**4
                        !     b2 = twohc/tab_lambda_nm(atom%continua(kr)%Nb+l-2)**4

                        !     wl = (tab_lambda_nm(atom%continua(kr)%Nb+l-1) - tab_lambda_nm(atom%continua(kr)%Nb+l-2))

                        !     tab_Aji_cont(kr,n,icell) = tab_Aji_cont(kr,n,icell) + wl * fourpi_h * 0.5*(anu1*e1*b1 + anu2*e2*b2)

                        ! enddo
                        !use cont%weights here
                        do l=1,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1
                            if (l==1) then
                               wl = 0.5*(tab_lambda_nm(atom%continua(kr)%Nb+1)-tab_lambda_nm(atom%continua(kr)%Nb)) / &
                                    tab_lambda_nm(atom%continua(kr)%Nb)
                            elseif (l==n_lambda) then
                               wl = 0.5*(tab_lambda_nm(atom%continua(kr)%Nr)-tab_lambda_nm(atom%continua(kr)%Nr-1)) / &
                                    tab_lambda_nm(atom%continua(kr)%Nr-1)
                            else
                               wl = 0.5*(tab_lambda_nm(atom%continua(kr)%Nb+l)-tab_lambda_nm(atom%continua(kr)%Nb+l-2)) / &
                                    tab_lambda_nm(atom%continua(kr)%Nb+l-1)
                            endif
                            anu1 = tab_Vij_cont(l,kr,n)

                            e1 = exp(-hc_k/T(icell)/tab_lambda_nm(atom%continua(kr)%Nb+l-1))

                            b1 = twohc/tab_lambda_nm(atom%continua(kr)%Nb+l-1)**3

                            tab_Aji_cont(kr,n,icell) = tab_Aji_cont(kr,n,icell) + fourpi_h * anu1*e1*b1 * wl

                        enddo
                    endif
                enddo
            enddo ct_loop
        enddo
        atom => null()

        !see + ionisation for all atoms ??
        if (n_iterate_ne > 0) then
            Neq_ne = NlevelTotal + 1 !adding charge conservation
            write(*,*) " There are ", Neq_ne, " equations for non-LTE SEE + ionisation"
            allocate(gtot(NlevelTotal,NlevelTotal,nb_proc),gr_save(NmaxLevel,NlevelTotal,NactiveAtoms,nb_proc))
            allocate(dgrdne(NlevelTotal,NlevelTotal,nb_proc),dgcdne(NlevelTotal,NlevelTotal,nb_proc))
            allocate(gs_ion(NmaxStage,NactiveAtoms),pops_ion(NmaxStage-1,NactiveAtoms,nb_proc))
            !pops_ion does not include neutral stage.
            do n=1,NactiveAtoms
                allocate(activeAtoms(n)%p%dgdne(activeAtoms(n)%p%Nlevel,activeAtoms(n)%p%Nlevel,nb_proc))
                mem_alloc_local = mem_alloc_local + sizeof(activeAtoms(n)%p%dgdne)
                gs_ion(1,n) = 1
                do j=2, activeatoms(n)%p%Nstage
                    gs_ion(j,n) = find_continuum(activeatoms(n)%p, gs_ion(j-1,n))
                enddo
                !index of the ground state of each ion.
                !if H: gs_ion(1) = 1; gs_ion(2) = nlevel
            enddo
            allocate(xvar(Neq_ne,nb_proc),fvar(Neq_ne,nb_proc),dfvar(Neq_ne,Neq_ne,nb_proc),npop_dag(Neq_ne,nb_proc))
            mem_alloc_local = mem_alloc_local + sizeof(xvar)*3 + sizeof(dfvar)+3*sizeof(gtot)+sizeof(gr_save)
            mem_alloc_local = mem_alloc_local + sizeof(gs_ion) + sizeof(pops_ion)
        endif


        mem_alloc_local = mem_alloc_local +  mem_alloc_non_local_alo
        write(*,'("  Total memory allocated in NLTEloop:"(1F14.3)" GB")') mem_alloc_local / 1024./1024./1024.
        write(*,'("  -- with "(1F14.3)" GB for the non-local rate matrix.")') mem_alloc_non_local_alo / 1024.**3

        return
    end subroutine alloc_nlte_var

    subroutine dealloc_nlte_var()
    !free space after non-lte loop.
        integer :: n, kr
        type (AtomType), pointer :: atom

        deallocate(psi,eta_atoms,chi_up,chi_down,uji_down)
        deallocate(itot,phi_loc,ngpop)

        deallocate(tab_Aji_cont, tab_Vij_cont)

        deallocate(tab_index_cell,tab_index_neighb,tab_index_i)
        if (allocated(lambda_ij)) deallocate(lambda_ij)
        if (allocated(psi_odiag)) deallocate(psi_odiag,Ujdown_odiag)!,eta_odiag)

        do n=1, NactiveAtoms
            atom => ActiveAtoms(n)%p
            deallocate(atom%Gamma,atom%g_nnz)
            do kr=1, atom%Nline
                deallocate(atom%lines(kr)%Rij,atom%lines(kr)%Rji)
                deallocate(atom%lines(kr)%Cij,atom%lines(kr)%Cji)
            enddo
            do kr=1, atom%Ncont
                deallocate(atom%continua(kr)%Rij,atom%continua(kr)%Rji)
                deallocate(atom%continua(kr)%Cij,atom%continua(kr)%Cji)
            enddo
            atom => null()
        enddo

        if (n_iterate_ne > 0) then
            deallocate(gtot,gr_save,dgrdne,dgcdne)
            deallocate(xvar,fvar,dfvar,npop_dag,gs_ion,pops_ion)
            do n=1,NactiveAtoms
                deallocate(activeAtoms(n)%p%dgdne)
            enddo
        endif

        return
    end subroutine dealloc_nlte_var


    subroutine see_atom(id,icell,atom,dM)
    ! --------------------------------------------------------------------!
    ! For atom atom solves for the Statistical Equilibrium Equations (SEE)
    ! with a digonal approximate operator.
    !
    !
    ! For numerical stability, the row with the largest populations is
    ! eliminated (i.e., it is replaced by the mass conservation law).
    !
    ! see Hubeny & Mihalas Stellar atmospheres, p. 448-450
    !	eqs 14.6 to 14.8c
    !
    !
    ! --------------------------------------------------------------------!

        integer, intent(in) :: icell, id
        type(AtomType), intent(inout) :: atom
        integer :: lp, imaxpop, l, Nsmall_pops, Nneg_or_null_pops
        integer :: level_index(atom%Nlevel)
        real(kind=dp), intent(out) :: dM
        real(kind=dp), dimension(atom%Nlevel) :: ndag, delta
        real(kind=dp) :: ntotal, Gamma_dag(atom%Nlevel,atom%Nlevel) !debug
        real(kind=dp) :: dn_n

        !if molecules (or H-), ntotal must be different from A*nHtot and found through chemical equilibrium.
        ntotal = atom%Abund*nHtot(icell) !in atomic form here. Could be stored on mem.
        ndag(:) = atom%n(:,icell) / ntotal
        !index of the equation to replace with mass conservation
        imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))

        !see Hubeny & Mihalas 2014 eq. 14.8a
        !calc the term in delta(l,l') (diagonal elements of rate matrix)
        !delta(l,l') = Gamma(l,l) = sum_l' R_ll' + C_ll' = sum_l' -Gamma(l',l)
        ! do l = 1, atom%Nlevel
        !     atom%Gamma(l,l,id) = 0.0_dp
        !     !Gamma(j,i) = -(Rij + Cij); Gamma(i,j) = -(Rji + Cji)
        !     !Gamma(i,i) = Rij + Cij = -sum(Gamma(j,i))
        !     atom%Gamma(l,l,id) = sum(-atom%Gamma(:,l,id)) !positive
        ! end do
        !check that diagonal is > 0 and off-diagonals < 0
        !call check_elements(id, icell, atom%Nlevel, atom%Gamma(:,:,id), "in see_atom: "//atom%ID)

        atom%n(:,icell) = 0d0
        atom%n(imaxpop,icell) = 1.0_dp

        Gamma_dag = atom%Gamma(:,:,id) !with diagonal
        atom%Gamma(imaxpop,:,id) = 1.0_dp

        ! solve for residual
        delta = atom%n(:,icell) - matmul(atom%Gamma(:,:,id), ndag)

        call GaussSlv(atom%Gamma(:,:,id), delta(:), atom%Nlevel)
   	    atom%n(:,icell) = ndag(:) + delta(:) ! * jacobi_sor
        ! call gaussslv(atom%Gamma(:,:,id), atom%n(:,icell), atom%Nlevel)


        if ((maxval(atom%n(:,icell)) < 0.0)) then
            write(*,*) ""
            write(*,*) atom%ID, " id=",id, " icell=",icell
            write(*,'("nstar: "*(ES14.5E3))') (ndag(l),l=1,atom%Nlevel) !*ntotal
            write(*,'("n: "*(ES14.5E3))') (atom%n(l,icell),l=1,atom%Nlevel) !*ntotal
            write(*,'("b: "*(ES14.5E3))') (atom%n(l,icell)/ndag(l),l=1,atom%Nlevel) !*ntotal
            do l=1, atom%Nlevel
                write(*, '(1I3, *(ES14.5E3))') l, (atom%Gamma(l,lp,id), lp=1, atom%Nlevel)
            enddo
            call warning("All pops are negative after SEE!")
            !or error ?
        endif

        if ((is_nan_infinity_vector(atom%n(:,icell))>0)) then
            write(*,*) atom%ID
            write(*,*) "(SEE) BUG pops", " id=",id, " icell=",icell
            write(*,*) "T=", t(icell), ' ne=', ne(icell)," nH=", nHtot(icell)
            write(*,'("ilevel: "*(1I4))') (l, l=1, atom%Nlevel)
            write(*,'("n: "*(ES14.5E3))') (atom%n(l,icell),l=1,atom%Nlevel)
            write(*,'("ndag: "*(ES14.5E3))') (ndag(l),l=1,atom%Nlevel)
            write(*,*) "Gamma:"
            write(*,'(*(I14))') (l, l=1, atom%Nlevel)
            do l=1, atom%Nlevel
                write(*, '(1I3, *(ES14.5E3))') l, (Gamma_dag(l,lp), lp=1, atom%Nlevel)
                write(*, '(1I3, "out:", *(ES14.5E3))') l, (atom%Gamma(l,lp,id), lp=1, atom%Nlevel)
            enddo
            write(*,*) "Rates"
            do l=1, atom%Nline
                write(*,*) " line ", atom%lines(l)%i, atom%lines(l)%j
                write(*,*) "-> Rij=",atom%lines(l)%Rij(id)," Rji=",atom%lines(l)%Rji(id)
                write(*,*) "-> Cij=",atom%lines(l)%Cij(id)," Cji=",atom%lines(l)%Cji(id)
            enddo
            do l=1, atom%Ncont
                write(*,*) " cont ", atom%continua(l)%i, atom%continua(l)%j
                write(*,*) "-> Rij=",atom%continua(l)%Rij(id)," Rji=",atom%continua(l)%Rji(id)
                write(*,*) "-> Cij=",atom%continua(l)%Cij(id)," Cji=",atom%continua(l)%Cji(id)
            enddo
            call error("nan or infinity found in n after SEE!")
        end if

        dM = 0.0_dp
        ndag = ndag * ntotal

        !Handle negative pops and very small populations
        Nsmall_pops = 0
        Nneg_or_null_pops = 0
        lp = 1
        do l=1,atom%Nlevel
            atom%n(l,icell) = atom%n(l,icell) * ntotal

            if (atom%n(l,icell) < frac_limit_pops * ntotal) then
                Nsmall_pops = Nsmall_pops + 1
                level_index(lp) = l
                lp = lp + 1
                if (atom%n(l,icell) <= 0.0_dp) then
                    atom%n(l,icell) = abs(atom%n(l,icell))
                    Nneg_or_null_pops = Nneg_or_null_pops + 1
                endif
            else
                dn_n = (1.0_dp - ndag(l))/atom%n(l,icell)
                dM = max(dM, abs(1.0_dp - ndag(l)/atom%n(l,icell)))
                ! dM = max(dM, abs(atom%n(l,icell)-ndag(l))/ndag(l))
            endif

        enddo
    ! 		if (Nsmall_pops >= 0.75*atom%Nlevel) then
    ! 			write(*,*) "*************** Warning SEE ***************"
    ! 			write(*,'("id #"(1I2)", cell #"(1I6)," ne "(ES14.5E3)" m^-3, nTot "(ES14.5E3)" m^-3")') id, icell, ne(icell), ntotal
    ! 			write(*,'(" for atom "(1A2)" with "(1I3)" levels")') atom%ID, atom%Nlevel
    ! 			write(*,'("--> Found "(1I3)" impurity levels "(1F12.5)" %")') Nsmall_pops, 100.0*real(Nsmall_pops)/real(atom%Nlevel)
    ! 			write(*,'(" --> with "(1I3)" below "(ES14.5E3)" m^-3")') Nneg_or_null_pops, prec_pops * ntotal
    ! 			write(*,'("ilevel: "*(1I4))') (level_index(lp), lp=1, Nsmall_pops)
    ! 			write(*,'("n: "*(ES14.5E3)" m^-3")') (atom%n(level_index(lp),icell), lp=1, Nsmall_pops)
    ! 			write(*,'("n/nTot: "*(ES14.5E3)" m^-3")') (atom%n(level_index(lp),icell)/ntotal, lp=1, Nsmall_pops)
    ! 			write(*,'("n/ne: "*(ES14.5E3)" m^-3")') (atom%n(level_index(lp),icell)/ne(icell), lp=1, Nsmall_pops)
    ! 			write(*,*) "*************** *********** ***************"
    ! 		endif

        ngpop(1:atom%Nlevel,atom%activeindex,icell,1) = atom%n(:,icell)
        atom%n(:,icell) = ndag

        return
    end subroutine see_atom

    subroutine nlocal_SEE_atom(at)
    ! --------------------------------------------------------------------!
    ! For atom at solves for the Statistical Equilibrium Equations (SEE)
    !
    !
    ! For numerical stability, the row with the largest populations is
    ! eliminated (i.e., it is replaced by the mass conservation law).
    !
    ! see Hubeny & Mihalas Stellar atmospheres, p. 448-450
    !	eqs 14.6 to 14.8c
    !
    !
    ! --------------------------------------------------------------------!
        type(AtomType), intent(inout) :: at
        real(kind=dp) :: dM
        integer :: niter, i, imax, j, l
        real(kind=dp) :: local_sol(at%Nlevel)
        real(kind=dp), allocatable :: b(:,:), ndag(:,:)

        !debug
        if (allocated(lambda_ij)) then
            open(unit=23,file='lambda_ij.b',form='unformatted',access='stream',status='unknown')
            write(23) n_cells
            write(23) Lambda_ij
            close(23)
        endif

        allocate(b(at%Nlevel,n_non_empty_cells),ndag(at%Nlevel,n_cells))
        b(:,:) = 0.0
        ndag(:,:) = at%n !to replace because the convergence is computed after inversion
        !using the new iterate computed locally with at%gamma (diagonal of rate matrix or at%g_nnz(:,:,:,1)).
        !as a starting solution.
        at%n = ngpop(1:at%Nlevel,at%activeindex,:,1)

        !for each cell, replace an equation with the mass conservation.
        do i=1,n_non_empty_cells
            !we replace the equation for the level with the highest population at the previous iteration, for each cell. 
            imax = locate(at%n(:,tab_index_cell(i)),maxval(at%n(:,tab_index_cell(i))))
            b(imax,i) = 1.0_dp
            at%n(:,tab_index_cell(i)) = at%n(:,tab_index_cell(i)) / (at%Abund * nHtot(tab_index_cell(i)))
            at%g_nnz(imax,:,i,1) = 1.0_dp
            ! do l = 1, at%Nlevel
            !     at%g_nnz(l,l,i,2:n_neighbours_max+1) = 0.0_dp
            !     at%g_nnz(l,l,i,2:n_neighbours_max+1) = sum(-at%g_nnz(:,l,i,2:n_neighbours_max+1),dim=1) !positive
            ! end do
            at%g_nnz(imax,:,i,2:n_neighbours_max+1) = 0.0
        enddo

        call jacobi_sparse_nlocal_see(at%Nlevel,n_non_empty_cells,n_neighbours_max,at%g_nnz,b,at%n,niter,dM)
        deallocate(b)
        do l=1, at%Nlevel
            at%n(l,:) = at%n(l,:) * at%Abund * nHtot
        enddo

        ngpop(1:at%Nlevel,at%activeindex,:,1) = at%n(:,:)
        at%n(:,:) = ndag !to change that in the future with non-local inversion.
        deallocate(ndag)

        return
    end subroutine nlocal_SEE_atom

    subroutine jacobi_sparse_nlocal_see(nl,n,m,nnz,b,x,niter,diff)
    ! Note: In principle could go to utils.f90. However, this routine is so specific
    ! to how we solve the non-local multi-level populations, that it cannot be used
    ! for general Ax = b systems.
    !
    ! Solve a set of linear equations Ax = b with the Jacobi  method. Where A is
    ! is a sparse matrix subset of B (nl*n, nl*n), of size (nl*n,nl*(m+1)) or (nl,nl,n,m+1).
    ! Here, x is the population vectors (nlevel,n_cells) and A, is a large sparse, 
    ! block matrix (nlevel,nlevel,n_non_empty_cells,n_neighbours_max+1) that contain
    ! a set of SEE for each cell taking into account the coupling with the neighbours.
    ! matrix vector multiplication of the type Ax, are therefore fine tuned for that problem.
    ! The diagonal is the concatenated diagonal of each SEE for each cell. that is:
    ! do i=1,n
    !     diag(1+(i-1)*nl:i*nl) = matdiag(nnz(:,:,i,1),nl)
    ! enddo
    !for a given cell, the diagonal is matdiag(nnz(:,:,i,1),nl). The off-diagonal elements (nnz(:,:,i,j>1))
    !have matdiag(nnz(:,:,i,j>1),nl) = 0.0 by definition.

        real(kind=dp), parameter :: tol = 1d-9 ! Convergence tolerance
        !eventually, omega would be defined as a function of the eigen values of A
        real(kind=dp), parameter :: omega = 1.0 ! Damping factor (0 < omega < 1)
        integer, parameter :: nIterMax = 500000

        integer, intent(in) :: nl, n, m
        real(kind=dp), intent(in) :: nnz(nl,nl,n,m+1)
        real(kind=dp), intent(inout) :: x(nl,*), b(nl,n)
        !to monitor convergence or not
        integer, intent(out) :: niter
        real(kind=dp), intent(out) :: diff

        logical :: lconverged
        integer :: i, j, icell, icell_n, l, imax, imax_l
        real(kind=dp) :: xtot(n), suma
        real(kind=dp) :: diag(nl), x_new(nl,n), y(nl), ap(nl,nl)


        diff = 0
        niter = 0
        lconverged = .false.

        do i=1, n
            icell = tab_index_cell(i)
            xtot(i) = sum(x(:,icell))
        enddo
        imax = locate(xtot,maxval(xtot))
        !replace one equation with mass conservation of all cells ?? 

        do while (.not. lconverged)
            niter = niter + 1

            ! Calculate the new x values using the damped Jacobi method

            suma = 0
            non_empty_loop: do i=1,n !loop over non-empty cells representing the populations of that cell
                if (i==imax) cycle non_empty_loop
                icell = tab_index_cell(i)

                y(:) = 0.0_dp 
                neighb_loop : do j=1,m!loop over (potential) neighbours of the cell i
                    icell_n = tab_index_neighb(i,j) !neighbour cell icell_n (j) of cell icell (i)
                    if (icell_n==0) cycle neighb_loop
                    y(:) = y(:) - matmul(nnz(:,:,i,j+1),x(:,icell_n))
                enddo neighb_loop

                x_new(:,i) = omega * (b(:,i) + y(:))
                ap(:,:) = nnz(:,:,i,1)
                call gaussslv(ap, x_new(:,i), nl)
                x_new(:,i) = x_new(:,i) + (1.0_dp - omega) * x(:,icell)
                
                write(*,*) "cell", i, icell
                write(*,*) "x_old:", x(:,icell)
                write(*,*) "x_new:", x_new(:,i)
                write(*,*) "ntot=", sum(x_new(:,i))

                suma = suma + sum(x(:,icell))
            enddo non_empty_loop
            icell = tab_index_cell(imax)
            y(:) = 0.0_dp
            do j=1,m!loop over (potential) neighbours of the cell i
                icell_n = tab_index_neighb(imax,j) !neighbour cell icell_n (j) of cell icell (i)
                if (icell_n==0) cycle
                y(:) = y(:) - matmul(nnz(:,:,imax,j+1),x(:,icell_n))
            enddo
            !replace the equation of conservation of mass of levels at cell imax with a global
            !mass conservation for all cells
            imax_l = locate(b(:,imax),maxval(b(:,imax)))
            b(imax_l,imax) = sum(xtot)
            diag(imax_l) = 1.0_dp
            y(imax_l) = y(imax_l) - suma
            write(*,*) 'imax_l:', imax_l
            write(*,*) 'b', b(:,imax)
            write(*,*) 'y=',y

            x_new(:,imax) = omega * (b(:,imax) + y(:))
            ap(:,:) = nnz(:,:,imax,1)
            call gaussslv(ap, x_new(:,imax), nl)
            x_new(:,imax) = x_new(:,imax) + (1.0_dp - omega) * x(:,icell)

            write(*,*) "cell", imax, icell
            write(*,*) "x_old:", x(:,icell)
            write(*,*) "x_new:", x_new(:,imax)
            write(*,*) "ntot=", sum(x_new(:,imax)), " nglob=", suma + x_new(imax_l,imax), sum(xtot)

            ! Calculate the error and check for convergence
            ! diff = maxval(abs(x_new - x)/x,mask=x>0)
            diff = 0.0
            do i=1,n
                diff = max( diff, maxval( abs( 1.0_dp - x(:,tab_index_cell(i)) / x_new(:,i) ) ) )
            enddo
            lconverged = (diff < tol)

            ! Update x for the next iteration
            do i=1,n
                x(:,tab_index_cell(i)) = x_new(:,i)
            enddo
            ! x = x_new
            if (niter > nIterMax) exit

            write(*, *) "Iteration ", niter, ": Error = ", diff
read(*,*)
        end do

        write(*, *) "(JACOBI) Iteration ", niter, ": Error = ", diff
stop
        return
    end subroutine jacobi_sparse_nlocal_see


    subroutine jacobi_sparse_nlocal_see_back(nl,n,m,nnz,b,x,niter,diff)
    ! Note: In principle could go to utils.f90. However, this routine is so specific
    ! to how we solve the non-local multi-level populations, that it cannot be used
    ! for general Ax = b systems.
    !
    ! Solve a set of linear equations Ax = b with the Jacobi  method. Where A is
    ! is a sparse matrix subset of B (nl*n, nl*n), of size (nl*n,nl*(m+1)) or (nl,nl,n,m+1).
    ! Here, x is the population vectors (nlevel,n_cells) and A, is a large sparse, 
    ! block matrix (nlevel,nlevel,n_non_empty_cells,n_neighbours_max+1) that contain
    ! a set of SEE for each cell taking into account the coupling with the neighbours.
    ! matrix vector multiplication of the type Ax, are therefore fine tuned for that problem.
    ! The diagonal is the concatenated diagonal of each SEE for each cell. that is:
    ! do i=1,n
    !     diag(1+(i-1)*nl:i*nl) = matdiag(nnz(:,:,i,1),nl)
    ! enddo
    !for a given cell, the diagonal is matdiag(nnz(:,:,i,1),nl). The off-diagonal elements (nnz(:,:,i,j>1))
    !have matdiag(nnz(:,:,i,j>1),nl) = 0.0 by definition.

        real(kind=dp), parameter :: tol = 1d-9 ! Convergence tolerance
        !eventually, omega would be defined as a function of the eigen values of A
        real(kind=dp), parameter :: omega = 1.0 ! Damping factor (0 < omega < 1)
        integer, parameter :: nIterMax = 500000

        integer, intent(in) :: nl, n, m
        real(kind=dp), intent(in) :: nnz(nl,nl,n,m+1)
        real(kind=dp), intent(inout) :: x(nl,*), b(nl,n)
        !to monitor convergence or not
        integer, intent(out) :: niter
        real(kind=dp), intent(out) :: diff

        logical :: lconverged
        integer :: i, j, icell, icell_n, l, imax
        real(kind=dp) :: diag(nl), x_new(nl,n), y(nl)


        diff = 0
        niter = 0
        lconverged = .false.

        do while (.not. lconverged)
            niter = niter + 1

            ! Calculate the new x values using the damped Jacobi method

            non_empty_loop: do i=1,n !loop over non-empty cells representing the populations of that cell
                if (i==imax) cycle non_empty_loop
                icell = tab_index_cell(i)

                !can be computed before hand, it is constant here.
                diag(:) = matdiag(nnz(:,:,i,1),nl)

                y(:) = diag(:) * x(:,icell) - matmul(nnz(:,:,i,1),x(:,icell))

                neighb_loop : do j=1,m!loop over (potential) neighbours of the cell i
                    icell_n = tab_index_neighb(i,j) !neighbour cell icell_n (j) of cell icell (i)
                    if (icell_n==0) cycle neighb_loop
                    y(:) = y(:) - matmul(nnz(:,:,i,j+1),x(:,icell_n))
                enddo neighb_loop

                x_new(:,i) = omega * (b(:,i) + y(:)) / diag(:) + (1.0_dp - omega) * x(:,icell)

            enddo non_empty_loop


            ! Calculate the error and check for convergence
            ! diff = maxval(abs(x_new - x)/x,mask=x>0)
            diff = 0.0
            do i=1,n
                diff = max( diff, maxval( abs( 1.0_dp - x(:,tab_index_cell(i)) / x_new(:,i) ) ) )
            !     do l=1, nl
            !         if (x_new(l,i) < frac_limit_pops * maxval(b(:,i))) then
            !             x_new(l,i) = abs(x_new(l,i))
            !         else
            !             diff = max( diff, abs( 1.0_dp - x(l,tab_index_cell(i)) / x_new(l,i) ) )
            !         endif
            !     enddo
            enddo
            lconverged = (diff < tol)

            ! Update x for the next iteration
            do i=1,n
                x(:,tab_index_cell(i)) = x_new(:,i)
            enddo
            ! x = x_new
            if (niter > nIterMax) exit

            write(*, *) "Iteration ", niter, ": Error = ", diff
        end do

        write(*, *) "(JACOBI) Iteration ", niter, ": Error = ", diff
stop
        return
    end subroutine jacobi_sparse_nlocal_see_back

    subroutine update_populations_nonlocal(iterate_ne)
    ! --------------------------------------------------------------------!
    ! Performs a solution of SEE for each atom and simultneously solves
    ! the charge conservation equations for hydrogen and helium.          !
    ! To do: implements Ng's acceleration iteration here
    ! --------------------------------------------------------------------!
        logical, intent(in) :: iterate_ne
        type(AtomType), pointer :: at
        integer :: nact

        !-> not yet
        ! if (iterate_ne) then
        !     !all atoms + ne at the same time
        !     ! call nlocal_see_atoms_ne()
        ! else
            do nact=1, Nactiveatoms
                at => activeatoms(nact)%p
                !the rate matrix is already filled
                call nlocal_SEE_atom(at)
            enddo
            at => null()
        ! endif

        return
    end subroutine update_populations_nonlocal

    subroutine rate_matrix(id)
        integer, intent(in) :: id
        integer :: nact

        do nact=1, NactiveAtoms
            call rate_matrix_atom(id, ActiveAtoms(nact)%p)
        enddo

        return
    end subroutine rate_matrix

    subroutine rate_matrix_atom(id, atom)
    !see Hubeny & Mihalas 2014 eq. 14.8a to 14.8c for the rate matrix elements.
    !Rij and Rji are 0 if lte.
    !atom%cswitch is 1.0 by default
    ! If lforce_lte (-see_lte) radiative rates are initialized at 0 and never updated.
    !TO DO:
    ! occupation probability
        integer, intent(in) :: id
        type (AtomType), intent(inout) :: atom
        integer :: kr, l, i, j

        atom%Gamma(:,:,id) = 0.0_dp

        do kr=1, atom%Nline

            i = atom%lines(kr)%i; j = atom%lines(kr)%j
            ! write(*,*) i, j
            ! write(*,*) "Rij=",atom%lines(kr)%Rij(id),"Cij=",atom%lines(kr)%Cij(id)
            ! write(*,*)"Rji=", atom%lines(kr)%Rji(id),"Cji=",atom%lines(kr)%Cji(id)

            atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - (atom%lines(kr)%Rij(id)+atom%cswitch*atom%lines(kr)%Cij(id))
            atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - (atom%lines(kr)%Rji(id)+atom%cswitch*atom%lines(kr)%Cji(id))

        enddo

        do kr=1, atom%Ncont

            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            ! write(*,*) i, j
            ! write(*,*) "Rij=",atom%continua(kr)%Rij(id),"Cij=",atom%continua(kr)%Cij(id)
            ! write(*,*) "Rji=",atom%continua(kr)%Rji(id),"Cji=",atom%continua(kr)%Cji(id)

            atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - (atom%continua(kr)%Rij(id)+atom%cswitch*atom%continua(kr)%Cij(id))
            atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - (atom%continua(kr)%Rji(id)+atom%cswitch*atom%continua(kr)%Cji(id))

        enddo

        do l = 1, atom%Nlevel
            atom%Gamma(l,l,id) = 0.0_dp
            !Gamma(j,i) = -(Rij + Cij); Gamma(i,j) = -(Rji + Cji)
            !Gamma(i,i) = Rij + Cij = -sum(Gamma(j,i))
            atom%Gamma(l,l,id) = sum(-atom%Gamma(:,l,id)) !positive
        end do

        !fill the diagonal of the non-zero version of the rate matrix
        atom%g_nnz(:,:,tab_index_i(cells_id(id)),1) = atom%Gamma(:,:,id)
        !j = 1 for the local point right ?
        ! do j=1, atom%Nlevel
        !     do i=1, atom%Nlevel
        !         atom%g_nnz(innz(atom%Nlevel,tab_index_i(cells_id(id)),j,1,i))
        !     enddo
        ! enddo

! stop
        return
    end subroutine rate_matrix_atom

    subroutine init_rates(id,icell)
        integer, intent(in) :: id, icell
        integer :: n

        if (alo_order > 0) then
            psi_odiag(:,:,id) = 0.0_dp
            Ujdown_odiag(:,:,:,:,id) = 0.0_dp
        endif

        do n=1,NactiveAtoms
            call init_radrates_atom(id,icell,ActiveAtoms(n)%p)

            !x ne included. Derivatives to ne not included.
            if (activeatoms(n)%p%id=='H') then
            !-> error if missing transitions, the rate are not added
                ! call init_colrates_atom(id,ActiveAtoms(n)%p)
                call collision_rates_hydrogen_loc(id,icell)
            else
                call init_colrates_atom(id,ActiveAtoms(n)%p)
            !-> error if missing transitions, the rate are not added
                call collision_rates_atom_loc(id,icell,ActiveAtoms(n)%p)
            endif

            activeatoms(n)%p%g_nnz(:,:,tab_index_i(cells_id(id)),:) = 0.0_dp
        enddo

        return
    end subroutine init_rates

    subroutine init_colrates_atom(id,atom)
        integer, intent(in) :: id
        type (atomtype), intent(inout) :: atom
        integer :: kr

        do kr=1,atom%Ncont
            atom%continua(kr)%Cij(id) = 0.0_dp
            atom%continua(kr)%Cji(id) = 0.0_dp
        enddo

        do kr=1,atom%Nline
            atom%lines(kr)%Cij(id) = 0.0_dp
            atom%lines(kr)%Cji(id) = 0.0_dp
        enddo

        return
    end subroutine init_colrates_atom

    subroutine init_radrates_atom(id,icell,atom)
        integer, intent(in) :: id, icell
        type (atomtype), intent(inout) :: atom
        integer :: kr, i, j
        real(kind=dp) :: nlte_fact
        nlte_fact = 1.0_dp
        if (lforce_lte) nlte_fact = 0.0_dp

        do kr=1,atom%Ncont
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            atom%continua(kr)%Rij(id) = 0.0_dp
            !updated value of ni and nj!
            !-> 0 if cont does not contribute to opac.
            atom%continua(kr)%Rji(id) = nlte_fact * tab_Aji_cont(kr,atom%activeindex,icell) * &
                 atom%nstar(i,icell)/atom%nstar(j,icell)
        enddo

        do kr=1,atom%Nline
            atom%lines(kr)%Rij(id) = 0.0_dp
            if (.not.atom%lines(kr)%lcontrib) then
                atom%lines(kr)%Rji(id) = 0.0_dp
                cycle
            endif
            atom%lines(kr)%Rji(id) = nlte_fact * atom%lines(kr)%Aji
        enddo

        return

    end subroutine init_radrates_atom

    subroutine accumulate_radrates(id, icell, iray, domega)
    ! --------------------------------------------------------- !
    ! Accumulate radiative rates for the MALI method.
    !  Rates are integrated in frequency and accumulated
    !  ray-by-ray.
    !
    ! BEWARE: Psi independent of iray (always index 1 in dim=2).
    ! 
    ! TO DO: - higher order frequency integration ---> need change
    !        in the cross-coupling and the init_radrates too.
    !        IT IS IMPORTANT to conserve the accuracy of all 
    !        frequency integration for stability.
    !        - level's dissolution.
    ! --------------------------------------------------------- !
        integer, intent(in) :: id, icell, iray
        real(kind=dp), intent(in) :: dOmega
        real(kind=dp), dimension(Nlambda_max_trans) :: ehnukt!, Ieff
        real(kind=dp), dimension(Nlambda_max_line) :: phi0
        real(kind=dp), dimension(n_lambda) :: Ieff
        type(AtomType), pointer :: at
        integer :: kr, i, j, Nl, Nr, Nb, nact
        integer :: ip, jp, Nrp, Nbp, krr
        real(kind=dp) :: jbar_up, jbar_down, xcc_down, xcc_up
        real(kind=dp) :: wphi, ni_on_nj_star, gij
        integer :: m, k

! TO DO: check the sum with at%g_nnz(j,i,k,m+1) as the access of neighbours is not working properly if icell_n = 0

        !ehnukt(:) = exp(hnu_k/T(icell)) --> defined only on tab_lambda_cont /= tab_lambda_nm !

        k = tab_index_i(icell)
        if (minval(Itot(:,iray,id)) < 0.0_dp) call error("(accumulate_radrates) Itot < 0!")

        atom_loop : do nact = 1, Nactiveatoms
            at => ActiveAtoms(nact)%p

            !Ieff = (Psi - Psi^\ast) eta
            Ieff(:) = Itot(:,iray,id)! - Psi(:,1,id) * eta_atoms(:,nact,id)
            do j=1, at%Nlevel
                Ieff(:) = max(Ieff(:) - psi(:,1,id) * Uji_down(:,j,nact,id) * at%n(j,icell), 0.0_dp)
            enddo
            !adding the neighbours' contribution
            if (alo_order > 0) Then
                do m=1, n_neighbours_max
                    if (tab_index_neighb(k,m)==0) cycle
                    do j=1, at%Nlevel
                        Ieff(:) = max(Ieff(:) - psi_odiag(:,m,id) * Ujdown_odiag(:,j,nact,m,id)*at%n(j,tab_index_neighb(k,m)),0.0_dp)
                    enddo
                enddo
            endif

            line_loop : do kr=1, at%Nline

                if (.not.at%lines(kr)%lcontrib) cycle line_loop

                Nr = at%lines(kr)%Nr; Nb = at%lines(kr)%Nb
                i = at%lines(kr)%i
                j = at%lines(kr)%j

                Nl = Nr - Nb + 1

                !get local profile for that cell id in the direction iray
                phi0(1:Nl) = phi_loc(1:Nl,kr,nact,iray,id)
                !-> otherwise phi_loc is 0 and there are normalisation issues with wphi
                ! due to a condition in the opacit y routine.
                if ((at%n(i,icell) - at%n(j,icell)*at%lines(kr)%gij) <= 0.0_dp) cycle line_loop

                Jbar_up = 0.0
                xcc_down = 0.0
                wphi = 0.0

                ! !(Psi - Psi^\ast) eta
                ! Ieff(1:Nl) = Itot(Nb:Nr,iray,id) - Psi(Nb:Nr,1,id) * eta_atoms(Nb:Nr,nact,id)

                Jbar_up = sum(Ieff(nb:nr)*phi0(1:Nl)*at%lines(kr)%wei)
                wphi = sum(phi0(1:Nl)*at%lines(kr)%wei)
                xcc_down = sum(chi_up(Nb:Nr,i,nact,id)*psi(Nb:Nr,1,id)*Uji_down(Nb:Nr,j,nact,id))

                if (wphi <= 0.0) then
                    call error("(accumulate_radrates) critical! wphi <= 0!")
                endif
                Jbar_up = Jbar_up / wphi

                jbar_down = jbar_up

                !init at Aji
                at%lines(kr)%Rji(id) = at%lines(kr)%Rji(id) + dOmega * Jbar_down * at%lines(kr)%Bji
                at%lines(kr)%Rij(id) = at%lines(kr)%Rij(id) + dOmega * Jbar_up * at%lines(kr)%Bij

                !accumulate off-diagonal terms
                if (alo_order > 0) then
                    do m=1, n_neighbours_max
                        if (tab_index_neighb(k,m)==0) cycle
                        ! + because G(i,j) = -R_ji but Rji = Aji - xcc_down
                        at%g_nnz(i,j,k,m+1) = at%g_nnz(i,j,k,m+1) + dOmega * &
                            sum(chi_up(Nb:Nr,i,nact,id)*psi_odiag(Nb:Nr,m,id)*Ujdown_odiag(Nb:Nr,j,nact,m,id))
                    enddo
                endif

                !cross-coupling terms + off-diagonal terms
                xcc_up = 0.0_dp
                do krr=1, at%Nline
                    ip = at%lines(krr)%i
                    jp = at%lines(krr)%j
                    Nbp = at%lines(krr)%Nb
                    Nrp = at%lines(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                        if (alo_order > 0) then
                            !Gji = - Rij = -( Rij + xcc_up * dOmega )
                            do m=1, n_neighbours_max
                                if (tab_index_neighb(k,m)==0) cycle
                                at%g_nnz(j,i,k,m+1) = at%g_nnz(j,i,k,m+1) - dOmega * &
                                    sum(chi_down(Nbp:Nrp,j,nact,id)*psi_odiag(Nbp:Nrp,m,id)*Ujdown_odiag(Nbp:Nrp,i,nact,m,id))
                            enddo
                        endif
                    endif
                enddo
                do krr=1, at%Ncont
                    ip = at%continua(krr)%i
                    jp = at%continua(krr)%j
                    Nbp = at%continua(krr)%Nb
                    Nrp = at%continua(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                        if (alo_order > 0) then
                            do m=1, n_neighbours_max
                                if (tab_index_neighb(k,m)==0) cycle
                                at%g_nnz(j,i,k,m+1) = at%g_nnz(j,i,k,m+1) - dOmega * &
                                    sum(chi_down(Nbp:Nrp,j,nact,id)*psi_odiag(Nbp:Nrp,m,id)*Ujdown_odiag(Nbp:Nrp,i,nact,m,id))
                            enddo
                        endif
                    endif
                enddo

                at%lines(kr)%Rji(id) = at%lines(kr)%Rji(id) - xcc_down * dOmega
                at%lines(kr)%Rij(id) = at%lines(kr)%Rij(id) + xcc_up * dOmega

                !debug
                if ((allocated(lambda_ij)).and.(kr==3)) then
                    lambda_ij(icell,icell) = lambda_ij(icell,icell) + dOmega * sum(psi(Nb:Nr,1,id))
                    do m=1, n_neighbours_max
                        if (tab_index_neighb(k,m)==0) cycle
                    !    write(*,*) "icell", icell, " k=", k
                    !    write(*,*) "icell_n", tab_index_neighb(k,m)
                        Lambda_ij(icell,tab_index_neighb(k,m)) = Lambda_ij(icell,tab_index_neighb(k,m)) + dOmega * sum(psi_odiag(Nb:Nr,m,id))
                    enddo
                endif

            end do line_loop

            cont_loop : do kr = 1, at%Ncont
                if (.not.at%continua(kr)%lcontrib) cycle cont_loop

                i = at%continua(kr)%i; j = at%continua(kr)%j

                ! ni_on_nj_star = ne(icell) * phi_T(icell, at%g(i)/at%g(j), at%E(j)-at%E(i))
                ni_on_nj_star = at%nstar(i,icell)/at%nstar(j,icell)
                gij = ni_on_nj_star * exp(-hc_k/T(icell)/at%continua(kr)%lambda0)
                if ((at%n(i,icell) - at%n(j,icell) * gij) <= 0.0_dp) cycle cont_loop

                Nb = at%continua(kr)%Nb; Nr = at%continua(kr)%Nr
                Nl = Nr - Nb + 1

                Jbar_down = 0.0
                Jbar_up = 0.0
                xcc_down = 0.0

                ! !(Psi - Psi^\ast) eta
                ! Ieff(1:Nl) = Itot(Nb:Nr,iray,id) - Psi(Nb:Nr,1,id) * eta_atoms(Nb:Nr,nact,id)
                ehnukt(1:Nl) = exp(-hc_k/tab_lambda_nm(nb:nr)/T(icell))

                !Note tab_Vij_cont contains the photoionisation cross-section for each continuum
                !on the non-LTE grid (tab_lambda_nm and not tab_lambda_cont). But, is has a size
                !of Nlambda_max_cont > Nl < n_lambda
                Jbar_up = sum(tab_Vij_cont(1:Nl,kr,nact)*Ieff(nb:nr)*at%continua(kr)%wei(:))
                !the term proportional to 2hnu^3/c2 already included in %Rji = "Aji" for continua.
                !exp(-hc_k/T(icell)/tab_lambda_nm(Nb:Nr))
                Jbar_down = sum(ehnukt(1:nl)*tab_Vij_cont(1:Nl,kr,nact)*Ieff(nb:nr)*at%continua(kr)%wei(:))
                xcc_down = sum(chi_up(Nb:Nr,i,nact,id)*Uji_down(Nb:Nr,j,nact,id)*psi(Nb:Nr,1,id))


                at%continua(kr)%Rij(id) = at%continua(kr)%Rij(id) + dOmega*fourpi_h * Jbar_up
                !init at tab_Aji_cont(kr,nact,icell) <=> Aji
                at%continua(kr)%Rji(id) = at%continua(kr)%Rji(id) + dOmega*fourpi_h * Jbar_down * ni_on_nj_star

                !accumulate off-diagonal terms
                if (alo_order > 0) then
                    do m=1, n_neighbours_max
                        if (tab_index_neighb(k,m)==0) cycle
                        ! + because G(i,j) = -R_ji but Rji = Aji - xcc_down
                        at%g_nnz(i,j,k,m+1) = at%g_nnz(i,j,k,m+1) + &
                            dOmega * sum(chi_up(Nb:Nr,i,nact,id)*psi_odiag(Nb:Nr,m,id)*Ujdown_odiag(Nb:Nr,j,nact,m,id))
                    enddo
                endif

                !cross-coupling terms
                xcc_up = 0.0_dp
                do krr=1, at%Nline
                    ip = at%lines(krr)%i
                    jp = at%lines(krr)%j
                    Nbp = at%lines(krr)%Nb
                    Nrp = at%lines(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                        if (alo_order > 0) then
                            do m=1, n_neighbours_max
                                if (tab_index_neighb(k,m)==0) cycle
                                at%g_nnz(j,i,k,m+1) = at%g_nnz(j,i,k,m+1) - dOmega * &
                                    sum(chi_down(Nbp:Nrp,j,nact,id)*psi_odiag(Nbp:Nrp,m,id)*Ujdown_odiag(Nbp:Nrp,i,nact,m,id))
                            enddo
                        endif
                    endif
                enddo
                do krr=1, at%Ncont
                    ip = at%continua(krr)%i
                    jp = at%continua(krr)%j
                    Nbp = at%continua(krr)%Nb
                    Nrp = at%continua(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                        if (alo_order > 0) then
                            do m=1, n_neighbours_max
                                if (tab_index_neighb(k,m)==0) cycle
                                at%g_nnz(j,i,k,m+1) = at%g_nnz(j,i,k,m+1) - dOmega * &
                                    sum(chi_down(Nbp:Nrp,j,nact,id)*psi_odiag(Nbp:Nrp,m,id)*Ujdown_odiag(Nbp:Nrp,i,nact,m,id))
                            enddo
                        endif
                    endif
                enddo

                at%continua(kr)%Rji(id) = at%continua(kr)%Rji(id) - xcc_down * dOmega
                at%continua(kr)%Rij(id) = at%continua(kr)%Rij(id) + xcc_up * dOmega

            enddo cont_loop
            at => NULL()

        end do atom_loop

        return
    end subroutine accumulate_radrates

    ! -> old, local and working routine.
    subroutine accumulate_radrates_mali(id, icell, iray, domega)
    ! --------------------------------------------------------- !
    ! Accumulate radiative rates for the MALI method.
    !  Rates are integrated in frequency and accumulated
    !  ray-by-ray.
    !
    !
    ! BEWARE: Psi independent of iray (so iray can be different of 1 during sub-iter)
    !
    ! TO DO: higher order frequency integ, level dissolution etc..
    ! --------------------------------------------------------- !
        integer, intent(in) :: id, icell, iray
        real(kind=dp), intent(in) :: dOmega
        real(kind=dp), dimension(Nlambda_max_trans) :: Ieff
        real(kind=dp), dimension(Nlambda_max_line) :: phi0
        type(AtomType), pointer :: atom, atom2
        integer :: kr, i, j, Nl, Nr, Nb, ip, jp, Nrp, Nbp
        integer :: i0, l, nact, krr
        real(kind=dp) :: jbar_up, jbar_down, xcc_down, xcc_up
        real(kind=dp) :: wl, wphi, anu, anu1, ni_on_nj_star, gij
        real(kind=dp) :: ehnukt, ehnukt1
! write(*,*) icell, T(icell)
!TO DO: check Itot >= 0 and ensure Ieff >= 0
        atom_loop : do nact = 1, Nactiveatoms
            atom => ActiveAtoms(nact)%p

            line_loop : do kr=1, atom%Nline

                if (.not.atom%lines(kr)%lcontrib) cycle line_loop

                Nr = atom%lines(kr)%Nr;Nb = atom%lines(kr)%Nb
                i = atom%lines(kr)%i
                j = atom%lines(kr)%j

                Nl = Nr - Nb + 1
                i0 = Nb - 1

                phi0(1:Nl) = phi_loc(1:Nl,kr,nact,iray,id)
                !-> otherwise phi_loc is 0 and there are normalisation issues with wphi
                if ((atom%n(i,icell) - atom%n(j,icell)*atom%lines(kr)%gij) <= 0.0_dp) cycle line_loop

                Jbar_up = 0.0
                xcc_down = 0.0
                wphi = 0.0

                !(Psi - Psi^\ast) eta
                Ieff(1:Nl) = Itot(Nb:Nr,iray,id) - Psi(Nb:Nr,1,id) * eta_atoms(Nb:Nr,nact,id)

                ! do l=2, nl
                !     wl = c_light * (tab_lambda_nm(i0+l) - tab_lambda_nm(i0+l-1))/atom%lines(kr)%lambda0
                !     Jbar_up = Jbar_up + 0.5 * (Ieff(l)*phi0(l)+Ieff(l-1)*phi0(l-1)) * wl
                !     wphi = wphi + 0.5*(phi0(l)+phi0(l-1)) * wl
                !     !xcc_down = xcc_down + 0.5 * ...
                ! enddo
                ! xcc_down = sum(chi_up(Nb:Nr,i,nact,id)*psi(Nb:nR,1,id)*Uji_down(Nb:Nr,j,nact,id))
                do l=1, nl
                    if (l==1) then
                        wl = 0.5*(tab_lambda_nm(Nb+1)-tab_lambda_nm(Nb)) * c_light / atom%lines(kr)%lambda0
                    elseif (l==nl) then
                        wl = 0.5*(tab_lambda_nm(Nr)-tab_lambda_nm(Nr-1)) * c_light / atom%lines(kr)%lambda0
                    else
                        wl = 0.5*(tab_lambda_nm(i0+l+1)-tab_lambda_nm(i0+l-1)) * c_light / atom%lines(kr)%lambda0
                    endif
                    Jbar_up = Jbar_up + Ieff(l)*phi0(l) * wl
                    wphi = wphi + phi0(l) * wl
                    xcc_down = xcc_down + chi_up(i0+l,i,nact,id)*psi(i0+l,1,id)*Uji_down(i0+l,j,nact,id)
                enddo

                if (wphi <= 0.0) then
                    call error("critical! wphi in accumulate_rates_mali!")
                endif
                Jbar_up = Jbar_up / wphi

                jbar_down = jbar_up

                !init at Aji
                atom%lines(kr)%Rji(id) = atom%lines(kr)%Rji(id) + dOmega * Jbar_down * atom%lines(kr)%Bji
                atom%lines(kr)%Rij(id) = atom%lines(kr)%Rij(id) + dOmega * Jbar_up * atom%lines(kr)%Bij

                !cross-coupling terms
                xcc_up = 0.0_dp
                do krr=1, atom%Nline
                    ip = atom%lines(krr)%i
                    jp = atom%lines(krr)%j
                    Nbp = atom%lines(krr)%Nb
                    Nrp = atom%lines(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo
                do krr=1, atom%Ncont
                    ip = atom%continua(krr)%i
                    jp = atom%continua(krr)%j
                    Nbp = atom%continua(krr)%Nb
                    Nrp = atom%continua(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo

                atom%lines(kr)%Rji(id) = atom%lines(kr)%Rji(id) - xcc_down * dOmega
                atom%lines(kr)%Rij(id) = atom%lines(kr)%Rij(id) + xcc_up * dOmega
! write(*,*) "line", kr, "x->down",xcc_down, " x->up", xcc_up
! write(*,*) " Jup", Jbar_up*domega, " Jdown", Jbar_down*domega
! if (kr==1) write(*,*) atom%lines(kr)%Rji(id), atom%lines(kr)%Aji, xcc_down * dOmega, Jbar_down * dOmega * atom%lines(kr)%Bji
            end do line_loop

            cont_loop : do kr = 1, atom%Ncont
                if (.not.atom%continua(kr)%lcontrib) cycle cont_loop

                i = atom%continua(kr)%i; j = atom%continua(kr)%j

          !ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
                ni_on_nj_star = atom%nstar(i,icell)/atom%nstar(j,icell)
                gij = ni_on_nj_star * exp(-hc_k/T(icell)/atom%continua(kr)%lambda0)
                if ((atom%n(i,icell) - atom%n(j,icell) * gij) <= 0.0_dp) cycle cont_loop

                Nb = atom%continua(kr)%Nb; Nr = atom%continua(kr)%Nr
                Nl = Nr-Nb + 1
                i0 = Nb - 1

                Jbar_down = 0.0
                Jbar_up = 0.0
                xcc_down = 0.0

                Ieff(1:Nl) = Itot(Nb:Nr,iray,id) - Psi(Nb:Nr,1,id) * eta_atoms(Nb:Nr,nact,id)
                ! write(*,*) Itot(Nb:Nr,iray,id)
                ! write(*,*) Psi(Nb:Nr,1,id)*eta_atoms(Nb:Nr,nact,id)
                ! write(*,*) Psi(Nb:Nr,1,id)
                ! write(*,*) eta_atoms(Nb:Nr,nact,id)
                ! stop

                ! do l=2, Nl

                !     anu = tab_Vij_cont(l,kr,nact)
                !     anu1 = tab_Vij_cont(l-1,kr,nact)

                !     wl = (tab_lambda_nm(i0+l) - tab_lambda_nm(i0+l-1))

                !     Jbar_up = Jbar_up + 0.5 * (anu*Ieff(l)/tab_lambda_nm(i0+l)+anu1*Ieff(l-1)/tab_lambda_nm(i0+l-1)) * wl

                !     ehnukt = exp(-hc_k/T(icell)/tab_lambda_nm(i0+l))/tab_lambda_nm(i0+l)
                !     ehnukt1 = exp(-hc_k/T(icell)/tab_lambda_nm(i0+l-1))/tab_lambda_nm(i0+l-1)

                !     Jbar_down = jbar_down + 0.5 * ( anu*Ieff(l)*ehnukt + anu1*Ieff(l-1)*ehnukt1 ) * wl

                ! enddo
                ! xcc_down = sum(chi_up(Nb:Nr,i,nact,id)*Uji_down(Nb:Nr,j,nact,id)*psi(Nb:Nr,1,id))

                do l=1, Nl
                    if (l==1) then
                        wl = 0.5*(tab_lambda_nm(Nb+1)-tab_lambda_nm(Nb)) / tab_lambda_nm(Nb)
                    elseif (l==n_lambda) then
                        wl = 0.5*(tab_lambda_nm(Nr)-tab_lambda_nm(Nr-1)) / tab_lambda_nm(Nr)
                    else
                        wl = 0.5*(tab_lambda_nm(i0+l+1)-tab_lambda_nm(i0+l-1)) / tab_lambda_nm(i0+l)
                    endif
                    anu = tab_Vij_cont(l,kr,nact)

                    Jbar_up = Jbar_up + anu*Ieff(l)*wl

                    ehnukt = exp(-hc_k/T(icell)/tab_lambda_nm(i0+l))

                    Jbar_down = jbar_down + anu*Ieff(l)*ehnukt*wl

                    xcc_down = xcc_down + chi_up(i0+l,i,nact,id)*Uji_down(i0+l,j,nact,id)*psi(i0+l,1,id)
                enddo


                atom%continua(kr)%Rij(id) = atom%continua(kr)%Rij(id) + dOmega*fourpi_h * Jbar_up
                !init at tab_Aji_cont(kr,nact,icell) <=> Aji
                atom%continua(kr)%Rji(id) = atom%continua(kr)%Rji(id) + dOmega*fourpi_h * Jbar_down * ni_on_nj_star


                !cross-coupling terms
                xcc_up = 0.0_dp
                do krr=1, atom%Nline
                    ip = atom%lines(krr)%i
                    jp = atom%lines(krr)%j
                    Nbp = atom%lines(krr)%Nb
                    Nrp = atom%lines(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                    ! if (kr==3) then
                    !     write(*,*) "xcc with line", krr, ip, jp
                    !     write(*,*) "val=", sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    !     write(*,*) tab_lambda_nm(Nbp), tab_lambda_nm(Nrp), atom%lines(krr)%voigt
                    !     write(*,*) sum(chi_down(Nbp:Nrp,j,nact,id)), maxval(Uji_down(Nbp:Nrp,i,nact,id)), maxval(psi(nbp:nrp,1,id))
                    !     write(*,*) Nbp, Nrp
                    !     ! do l=Nbp,Nrp
                    !     !     write(*,*) tab_lambda_nm(l),psi(l,1,id)
                    !     ! enddo
                    !     ! stop
                    ! endif
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo
                do krr=1, atom%Ncont
                    ip = atom%continua(krr)%i
                    jp = atom%continua(krr)%j
                    Nbp = atom%continua(krr)%Nb
                    Nrp = atom%continua(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                    ! if (kr==3) then
                    !     write(*,*) "xcc with cont", krr, ip, jp
                    ! endif
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo

                atom%continua(kr)%Rji(id) = atom%continua(kr)%Rji(id) - xcc_down * dOmega
                atom%continua(kr)%Rij(id) = atom%continua(kr)%Rij(id) + xcc_up * dOmega
! write(*,*) "cont", kr, "x->down",xcc_down, " x->up", xcc_up
! write(*,*) " Jup", Jbar_up, " Jdown", Jbar_down, "uji = ", tab_Aji_cont(kr,nact,icell)

            enddo cont_loop
            atom => NULL()

        end do atom_loop
! stop
        return
    end subroutine accumulate_radrates_mali

    subroutine update_populations(id, icell, iterate_ne, delta)
    ! --------------------------------------------------------------------!
    ! Performs a solution of SEE for each atom and simultneously solves
    ! the charge conservation equations for hydrogen and helium.          !
    ! To do: implements Ng's acceleration iteration here
    ! --------------------------------------------------------------------!
        integer, intent(in) :: id, icell
        logical, intent(in) :: iterate_ne
        type(AtomType), pointer :: atom
        integer :: nact
        real(kind=dp) :: dM, dne
        real(kind=dp), intent(out) :: delta

        delta = 0.0
        if (iterate_ne) then
            !all atoms + ne at the same time
            call see_atoms_ne(id,icell,delta,dne)
        else
            do nact=1, Nactiveatoms
                atom => activeatoms(nact)%p
                call rate_matrix_atom(id, atom)
                call SEE_atom(id, icell, atom, dM)
                !max for all atoms
                delta = max(delta, dM)
                atom => null()
            enddo
        endif
        !do nact=1,NactiveAtoms
        ! call calc_Tex_atom(icell, ActiveAtoms(nact)%p, dT, Tex, Tion)
        !enddo

        return
    end subroutine update_populations

   subroutine check_radrates(id,icell,verbose,at)
        integer, intent(in) :: id,icell
        logical, intent(in) :: verbose
        integer :: l
        type(AtomType), pointer, intent(inout) :: at
        real(kind=dp), parameter :: rate_limit = 1d-15

        do l=1, at%Nline
            if (at%lines(l)%Rij(id)<rate_limit) then
                if (verbose) then
                    write(*,*) "WARNING at cell", icell, id, " for atom", at%ID
                    write(*,*) "Rij of line", at%lines(l)%i, at%lines(l)%j, " is negative!"
                endif
                at%lines(l)%Rij(id) = 0.0_dp
            endif
            if (at%lines(l)%Rji(id)<rate_limit) then
                if (verbose) then
                    write(*,*) "WARNING at cell", icell, id
                    write(*,*) "Rji of line", at%lines(l)%i, at%lines(l)%j, " is negative!"
                endif
                at%lines(l)%Rji(id) = 0.0_dp
            endif
        enddo
        do l=1, at%Ncont
            if (at%continua(l)%Rij(id)<rate_limit) then
                if (verbose) then
                    write(*,*) "WARNING at cell", icell, id, " for atom", at%ID
                    write(*,*) "Rij of cont", at%continua(l)%i, at%continua(l)%j, " is negative!"
                endif
                at%continua(l)%Rij(id) = 0.0_dp
            endif
            if (at%continua(l)%Rji(id)<rate_limit) then
                if (verbose) then
                    write(*,*) "WARNING at cell", icell, id
                    write(*,*) "Rji of cont", at%continua(l)%i, at%continua(l)%j, " is negative!"
                endif
                at%continua(l)%Rji(id) = 0.0_dp
            endif
        enddo

        return
    end subroutine check_radrates

    subroutine see_atoms_ne(id,icell,dM,dne)
    !
    ! Solve coupled sets of SEE for all atoms including
    !  charge conservation equation.
    !
        integer, intent(in) :: id, icell
        real(kind=dp), intent(out) :: dm, dne

        integer, parameter :: maxIter = 1000
        integer :: n_iter, max_iter
        integer, parameter :: damp_scaling = 10
        real(kind=dp), parameter :: damp_char = 5.0, precision = 1d-5
        logical ::verbose, lconverged, rest_damping, lsecond_try, neg_pops
        real(kind=dp) :: d_damp, tmp_fact

        type (AtomType), pointer :: at
        integer :: n, nlev, i, j, kr, ip, jp
        real(kind=dp) ::delta_f, dfpop,dfne

        verbose = .false.!debug mode
        tmp_fact = 1.0_dp
        if (lforce_lte) tmp_fact = 0.0_dp

        nlev = Neq_ne - 1

        dM = 0.0
        dne = 0.0
        lconverged = .false.
        rest_damping = .false.
        d_damp = damp_char
        lsecond_try = .false. !set to .true. to avoid starting again if the first iterative scheme failed

        dgrdne(:,:,id) = 0.0_dp
        gr_save(:,:,:,id) = 0.0_dp
        npop_dag(:,id) = 0.0_dp

        !gather informations for all active atoms in the same arrays
        !store populations and densities of the cell icell of this MALI iteration.
        npop_dag(Neq_ne,id) = ne(icell)
        i = 1
        do n=1,nactiveatoms
            at => ActiveAtoms(n)%p
            npop_dag(i:(i-1)+at%Nlevel,id) = at%n(:,icell)
            ! !-> should not be needed...
            ! call check_radrates(id,icell,verbose,at)
            call radrate_matrix_on_ne_atom(id,icell,at,at%gamma(:,:,id),at%dgdne(:,:,id))
            dgrdne(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel,id) = tmp_fact*at%dgdne(:,:,id)
            gr_save(1:at%Nlevel,1:at%Nlevel,n,id) = tmp_fact*at%gamma(:,:,id)
            i = i + at%Nlevel
        enddo
        at => null()

        !start Newton-Raphson iterations
        if (verbose) write(*,*) "T = ", T(icell), " nHtot = ", nHtot(icell), " ne=", ne(icell)

        n_iter = 0
        max_iter = maxIter

        iterative_loop : do while (.not.lconverged)

            n_iter = n_iter + 1
            xvar(1:nlev,id) = 0.0
            xvar(Neq_ne,id) = ne(icell) !current value
            fvar(:,id) = 0.0
            dfvar(:,:,id) = 0.0
            pops_ion(:,:,id) = 0.0
            gtot(:,:,id) = 0.0_dp
            dgcdne(:,:,id) = 0.0_dp

            !initiate with CURRENT values
            !get the total populations of each ion except the neutrals!
            !ionisation fraction is
            !do j=1,nstage-1
            ! j * pops_ion(j) (pops_ion(1)=singly ionised, pops_ion(2)=twice ionised)
            !enddo
            i = 1
            do n=1,nactiveatoms
                at => ActiveAtoms(n)%p
                xvar(i:(i-1)+at%Nlevel,id) = at%n(:,icell) !current values
                ! write(*,*) ""
                ! write(*,*) at%Id, at%Nstage, gs_ion(1,n)
                do j=2, at%Nstage-1
                    ! write(*,*) gs_ion(j,n),gs_ion(j+1,n)
                    pops_ion(j-1,n,id) = sum(at%n(gs_ion(j,n):gs_ion(j+1,n)-1,icell))
                    ! write(*,*) pops_ion(j-1,n,id)
                enddo
                pops_ion(at%Nstage-1,n,id) = at%n(at%Nlevel,icell)
                ! write(*,*) pops_ion(at%Nstage-1,n,id),  at%n(at%Nlevel,icell)
                if (at%ID=='H') then
                    ! call init_colrates_atom(id,at)
                    call collision_rates_hydrogen_loc(id,icell)
                else
                    call init_colrates_atom(id,at)
                    call collision_rates_atom_loc(id,icell,at)
                endif
!*****!
                !check that or need a temporary array here ??
!******!
                !add current ne estimate to continua radiative rates
                gtot(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel,id) = gr_save(1:at%Nlevel,1:at%Nlevel,n,id)
                do kr=at%Ntr_line+1,at%Ntr
                    ip = at%i_trans(kr); jp = at%j_trans(kr)
                    gtot(i+ip-1,i+jp-1,id) = gtot(i+ip-1,i+jp-1,id) * xvar(neq_ne,id)
                enddo
                !recompute diagonal for that atom now that ne is taken intout account in gam_rad
                !can be done better
                do j=1,at%Nlevel
                    gtot(i-1+j,i-1+j,id) = 0.0_dp !set to 0 because the diag id already filled
                    gtot(i-1+j,i-1+j,id) = -sum(gtot(i:(i-1)+at%Nlevel,i-1+j,id))
                enddo

                call colrate_matrix_atom(id, icell, at, at%gamma(:,:,id),at%dgdne(:,:,id))
                dgcdne(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel,id) = at%dgdne(:,:,id)
                gtot(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel,id) = gtot(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel,id) + at%Gamma(:,:,id)
                i = i + at%nlevel
            enddo
            at=>null()
            !rate equation for this atom stored in f and df !
            !an equation per level dn_i/dt = 0 = sum_lp n_l * Gamma_lp_l
            call rate_equations(id, icell, Neq_ne, gtot(:,:,id), dgrdne(:,:,id), dgcdne(:,:,id), xvar(:,id), &
                 fvar(:,id), dfvar(:,:,id))

            !charge conservation!
            call non_lte_charge_conservation (id,icell, neq_ne, xvar(:,id), fvar(:,id), dfvar(:,:,id))

            !replace one equation of SEE by particle number conservation
            !particule conservation!
            call particle_conservation (icell, Neq_ne, xvar(:,id), fvar(:,id), dfvar(:,:,id))

            !newton raphson!
            call multivariate_newton_raphson (neq_ne, dfvar(:,:,id), fvar(:,id), xvar(:,id))

            !update atomic populations and ne
            neg_pops = .false.
            i = 1
            do n=1,NactiveAtoms
                at => ActiveAtoms(n)%p
                do j=1,at%Nlevel
                    at%n(j,icell) = at%n(j,icell) * ( 1.0_dp + fvar((i-1)+j,id)/(1.0_dp + d_damp * abs(fvar((i-1)+j,id))) )
                    if (at%n(j,icell) < 0.0) neg_pops = .true.
                    ! if (at%n(j,icell) < frac_limit_pops * at%Abund*nHtot(icell) )then
                    !     write(*,*) "small pops for level ", j
                    ! endif
                enddo
                i = i + at%Nlevel
            enddo
            at => null()
            if (verbose) write(*,*) ne(icell)
            ne(icell) = ne(icell) * ( 1.0 + fvar(neq_ne,id)/(1.0_dp + d_damp * abs(fvar(neq_ne,id))) )
            if (verbose)write(*,*) ( 1.0 + fvar(neq_ne,id)/(1.0_dp + d_damp * abs(fvar(neq_ne,id))) )
            if (verbose)write(*,*) fvar(neq_ne,id), d_damp
            !                       1d-16
            ! if (ne(icell) < frac_limit_pops * nhtot(icell)) write(*,*) "** small ne at cell ", icell
            ! if ( (ne(icell) < frac_limit_pops * nHtot(icell)).or.(neg_pops) ) rest_damping = .true.
            if ( (ne(icell) < frac_limit_pops * sum(pops_ion(:,:,id))).or.(neg_pops) ) rest_damping = .true.
            !-> here pops_ion has the old values !
            !restart with more iterations and larger damping (more stable, slower convergence)
            if (rest_damping .and. d_damp < (damp_char + 1.0)) then
                xvar(:,id) = npop_dag(:,id)
                d_damp = damp_char * damp_scaling
                rest_damping = .false.
                n_iter = 0
                max_iter = damp_scaling * maxIter !can be large !

            elseif (rest_damping) then
                neg_pops = .false.
                lconverged = .false.
            endif

            !should be fractional
            delta_f = maxval(abs(fvar(:,id)))
            dfpop = maxval(abs(fvar(1:Neq_ne-1,id)))
            dfne = abs(fvar(neq_ne,id))

            if (verbose) then
                write(*,'("niter #"(1I5))') n_iter
                write(*,'("non-LTE ionisation delta="(1ES17.8E3)" dfpop="(1ES17.8E3)" dfne="(1ES17.8E3))') delta_f, dfpop, dfne
                write(*,'("non-LTE ionisation ne="(1ES17.8E3)" m^-3; nedag="(1ES17.8E3)" m^-3")') ne(icell), npop_dag(Neq_ne,id)
                ! stop
            endif

            if (n_iter > max_iter) then
                if (lsecond_try) then
                    lconverged = .true.
                    !! TO DO: add a warning here when verbose is .false.
                    !! otherwise we don't know if there is a problem...
                    !! and it is dangerous/unethical to hide that(?).
                    if (verbose) then
                        write(*,*) ""
                        write(*,*) "-> (non-LTE ionisation): Not enough iterations to converge after second try"
                        write(*,'("cell #"(1I7)"; proc #"(1I3)", niter #"(1I5))') icell, id, n_iter
                        write(*,'("maxIter: "(1I5)"; damp="(1I4))')  max_iter, nint(d_damp)
                        write(*,'("non-LTE ionisation delta="(1ES17.8E3)" dfpop="(1ES17.8E3)" dfne="(1ES17.8E3))') &
                             delta_f, dfpop, dfne
                        write(*,*) ""
                    endif
                    if (verbose) stop
                else
                    lsecond_try = .true.
                    lconverged = .false.
                    n_iter = 0
                    max_iter = damp_scaling * maxIter
                    d_damp = real(damp_scaling) * damp_char
                    xvar(:,id) = npop_dag(:,id)
                    if (verbose) write(*,*) "second iter", max_iter, d_damp
                endif
            endif

            if ((delta_f < precision).and..not.rest_damping) then
                lconverged = .true.
                exit iterative_loop
            endif

        enddo iterative_loop !on convergence

        !compute global convergence rate and reset values (ne, n) to
        ! the value of the MALI iteration.
        ! the new ne and n are updated only after the end of a MALI iteration (loop over all cells) to avoid
        ! introducing to much non-linearity.

        dne = abs(1.0_dp - npop_dag(Neq_ne,id)/ne(icell))
        dM = 0
        i = 1
        do n=1,NactiveAtoms
            at => Activeatoms(n)%p
            dM = max(dM,maxval(abs(1.0 - npop_dag(i:(i-1)+at%Nlevel,id)/at%n(:,icell))))
            ngpop(1:at%Nlevel,at%activeindex,icell,1) = at%n(:,icell)
            !reset
            at%n(:,icell) = npop_dag(i:(i-1)+at%Nlevel,id) !the first value before iterations
            i = i + at%Nlevel
        enddo
        at => null()

        if (verbose) write(*,'("(DELTA) non-LTE ionisation dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne

        ! ne_new(icell) = ne(icell)
        ngpop(1,NactiveAtoms+1,icell,1) = ne(icell)
        ne(icell) = npop_dag(Neq_ne,id)


        return
    end subroutine see_atoms_ne

  subroutine ionisation_frac_lte(elem, k, ne, fjk, dfjk)
  !special version that forces ionisation fraction at LTE.
    real(kind=dp), intent(in) :: ne
    integer, intent(in) :: k
    type (Element), intent(in) :: Elem
    real(kind=dp), dimension(:), intent(inout) :: fjk, dfjk
    real(kind=dp) :: Uk, Ukp1, sum1, sum2
    integer :: j,  Nstage
    !return Nj / Ntot
    !for neutral j * Nj/Ntot = 0

    !fjk(1) is N(1) / ntot !N(1) = sum_j=1 sum_i=1^N(j) n(i)
    fjk(1)=1.
    dfjk(1)=0.
    sum1 = 1.
    sum2 = 0.
    Uk = get_pf(elem,1,T(k))
    do j=2,Elem%Nstage
       Ukp1 = get_pf(elem,j,T(k))
       fjk(j) = Sahaeq(T(k),fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)
       !fjk(j) = ne * cste, der = cste  fjk(j) / ne

       dfjk(j) = -(j-1)*fjk(j)/ne

       sum1 = sum1 + fjk(j)
       sum2 = sum2 + dfjk(j)

       Uk = Ukp1
    end do

    fjk(:)=fjk(:)/sum1
    dfjk(:)=(dfjk(:)-fjk(:)*sum2)/sum1


    return
  end subroutine ionisation_frac_lte

  function jions(id)
  !sum of ions times charge for nlte ions
    real(kind=dp) :: jions
    integer, intent(in) :: id
    real(kind=dp) :: sum_ions
    integer :: n, j

    jions = 0
    do n=1,Nactiveatoms
        do j=1, activeatoms(n)%p%Nstage-1
            jions = jions + j * pops_ion(j,n,id)
        enddo
    enddo

    return
  endfunction jions

    subroutine non_lte_charge_conservation (id,icell, neq, x, f, df)
    !F_cc = 1.0 - 1/ne * (np + nHeII + 2*nHeIII) = 0
    !F_cc = 1.0 - 1/ne * (sum_atom pops_ion(atom)*stage(ion)) = 0
    !also computes the derivative of the charge conservation equation with
    !ne and total population of stages with stage > 0 (like nH(Nlevel)=nHII)
        integer, intent(in) :: id, icell, neq
        real(kind=dp), intent(in) :: x(neq)
        real(kind=dp), intent(inout) :: df(neq,neq), f(neq)
        integer :: n, j, i, k
        type (element) :: elem
        real(kind=dp) :: akj, sum_ions, st
        real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk

        ! sum_ions = 0
        ! do n=1,Nactiveatoms
        !     do j=1, activeatoms(n)%p%Nstage-1
        !         sum_ions = sum_ions + j * pops_ion(j,n,id)
        !     enddo
        ! enddo
        sum_ions = jions(id)

        !derivative of CC to ne
        f(neq) = 1.0_dp - sum_ions / x(neq) !factor 2 for nHeIII included (1 for np and nheII)
        df(neq,neq) = sum_ions / x(neq)**2

        !derivative of CC to ions and individual levels making up to the ions pops.
        !in otherword, derivative to x(1:neq-1) if x belongs to an ion (stage > 0, not neutral)
        !taking into account the ionisation stage.
        !e.g. dne/dnp = -1/ne
        !      dne/dnHeIII = -2/ne
        !     dne/nHeII_i = -1/ne
        ! which means that ions with the same ionisation stage gives the same derivative.
        !I can do better here
        i = 1
        do n=1,NactiveAtoms
            do k=1, activeatoms(n)%p%Nlevel
                st = activeatoms(n)%p%stage(k)
                !for neutrals it adds 0 !
                df(neq,(i-1)+k) = -st / x(neq)
                !for debug only
                if ((i-1)+k == neq) call error("non_lte_charge_conservation!")
            enddo
            i = i + activeatoms(n)%p%Nlevel
        enddo


        ! f(neq) = 1.0
        ! df(neq,neq) = 0
        ! df(neq,:) = 0.0
        !now contribution from bacgrkound atoms.
        lte_elem_loop : do n=1, 26 !for 26 elements, can go up to Nelem(size(Elems)) in principle!
            elem = Elems(n)
            if (elem%nm>0) then
                if (atoms(elem%nm)%p%active) cycle lte_elem_loop
            endif

            call ionisation_frac_lte(elem, icell, x(neq), fjk, dfjk)

            do j=2, elem%Nstage
                !pure LTE term j = 1 corresponds to the first ionisation stage always
                !unlike in solve_ne where this loop is also for non-LTE models (with non-LTE ionisation fraction, still)
                akj = (j-1) * elem%abund * nhtot(icell) !(j-1) = 0 for neutrals and 1 for singly ionised
                f(neq) = f(neq) - akj * fjk(j) / x(neq)
                df(neq,neq) = df(neq,neq) + akj * fjk(j) / x(neq)**2 - akj / x(neq) * dfjk(j)
            end do

        enddo lte_elem_loop

        return
    end subroutine non_lte_charge_conservation

    subroutine radrate_matrix_on_ne_atom(id, icell, atom, gr, dgrdne)

    !Rji continua are divided by ne !

        integer, intent(in) :: icell, id
        type (AtomType), intent(in) :: atom
        real(kind=dp), intent(out), dimension(atom%Nlevel, atom%Nlevel) :: gr, dgrdne
        integer :: kr, i, j

        gr(:,:) = 0.0_dp
        dgrdne(:,:) = 0.0_dp !derivative of recombinaison rate matrix to dne

        do kr=1, atom%Nline
            i = atom%lines(kr)%i; j = atom%lines(kr)%j
            gr(j,i) = gr(j,i) - atom%lines(kr)%Rij(id)
            gr(i,j) = gr(i,j) - atom%lines(kr)%Rji(id)
            ! write(*,*) "line", i, j, "Rij=",atom%lines(kr)%Rij(id), "Rji=",atom%lines(kr)%Rji(id)
        enddo

        do kr=1, atom%Ncont
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            gr(j,i) = gr(j,i) - atom%continua(kr)%Rij(id)
            gr(i,j) = gr(i,j) - atom%continua(kr)%Rji(id) / ne(icell)
            dgrdne(i,j) = dgrdne(i,j) - atom%continua(kr)%Rji(id) / ne(icell)
            ! write(*,*) "cont", i, j,  "Rij=",atom%continua(kr)%Rij(id),  "Rji=",atom%continua(kr)%Rji(id)
        enddo

        do i=1, atom%Nlevel
            gr(i,i) = -sum(gr(:,i))
            dgrdne(i,i) = -sum(dgrdne(:,i))
        enddo

        return
    end subroutine radrate_matrix_on_ne_atom

    subroutine colrate_matrix_atom(id, icell, atom, gc, dgcdne)

    !see Hubeny & Mihalas 2014 eq. 14.8b

        integer, intent(in) :: id, icell
        type (AtomType), intent(in) :: atom!, pointer
        real(kind=dp), intent(out), dimension(atom%Nlevel,atom%Nlevel) :: gc, dgcdne
        integer :: kr, i, j

        gc(:,:) = 0.0_dp
        dgcdne(:,:) = 0.0_dp

        do kr=1, atom%Nline
            i = atom%lines(kr)%i; j = atom%lines(kr)%j
            gc(j,i) = gc(j,i) - atom%lines(kr)%Cij(id)
            gc(i,j) = gc(i,j) - atom%lines(kr)%Cji(id)
            dgcdne(i,j) = dgcdne(i,j) - atom%lines(kr)%Cji(id) / ne(icell) !no factor 2 here, bound-bound
            dgcdne(j,i) = dgcdne(j,i) - atom%lines(kr)%Cij(id) / ne(icell)
        enddo

        do kr=1, atom%Ncont
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            gc(j,i) = gc(j,i) - atom%continua(kr)%Cij(id)
            gc(i,j) = gc(i,j) - atom%continua(kr)%Cji(id)
            dgcdne(j,i) = dgcdne(j,i) - atom%continua(kr)%Cij(id) / ne(icell)
            dgcdne(i,j) = dgcdne(i,j) - 2.0 * atom%continua(kr)%Cji(id) / ne(icell)
        enddo

        do i=1, atom%Nlevel
            gc(i,i) = -sum(gc(:,i))
            dgcdne(i,i) = -sum(dgcdne(:,i))
        enddo

        if (atom%cswitch>1.0) then
            gc(:,:) = gc(:,:) * atom%cswitch
            dgcdne(:,:) = dgcdne(:,:) * atom%cswitch
        endif

        return
    end subroutine colrate_matrix_atom

    subroutine multivariate_newton_raphson (neq, df, f, x)

    !matmul(Jacobian, delta) = - f_equation_vector

        integer, intent(in) :: neq
        real(kind=dp), intent(in) :: x(neq)
        real(kind=dp), intent(inout) :: df(neq,neq), f(neq)
        integer :: ieq, jvar
        ! real(kind=dp) :: Adag(neq,neq), bdag(neq), res(neq)

        do ieq=1, neq
            f(ieq) = f(ieq) * -1.0_dp
            do jvar=1, neq
                df(ieq,jvar) = df(ieq,jvar) * x(jvar)
            enddo
        enddo

        ! *********************** !
        ! Adag(:,:) = df(:,:)
        ! bdag(:) = f(:)
        ! *********************** !

        ! call GaussSlv(df,f,neq)
        call solve_lin(df,f,neq)
        if (is_nan_infinity_vector(f)>0) then
        !  write(*,*) "fdag=", bdag
        !  write(*,*) "f=", f
         call error("(Newton-raphson) nan in (f,df) after GSLV")
        endif

        ! *********************** !
        !compute difference between new solution and old one
        ! res(:) = bdag(:)
        ! do ieq=1,Neq
        !     do jvar=1,Neq
        !         res(ieq) = res(ieq) - Adag(ieq,jvar)*f(jvar)
        !     enddo
        ! enddo
        ! res(:) = bdag(:) - matmul(Adag,f)

        ! !solve for the residual
        ! call Gaussslv(Adag, res, Neq)

        ! f(:) = f(:) + res(:)
        ! *********************** !

        return
    end subroutine multivariate_newton_raphson

    subroutine rate_equations(id,icell,neq,gam,dgrdne,dgcdne,x,f,df)
        integer, intent(in) :: id, icell, neq
        real(kind=dp), intent(in) :: x(neq)
        real(kind=dp), intent(in), dimension(neq-1,neq-1) :: gam,dgrdne, dgcdne
        real(kind=dp), intent(out) ::  f(neq), df(neq,neq)
        integer :: n, i
        type (AtomType), pointer :: at

        !derivative of SEE to ne
        !dCij/dne = Cij/ne
        !dCji/dne = 2*Cji/ne
        !dRji_cont/dne = Rji_cont/ne

        Nlevels = neq - 1
        !derivative of the equation of all levels to that level
        df(1:Nlevels,1:Nlevels) = gam(:,:)
        i = 1
        do n=1, Nactiveatoms
            at => ActiveAtoms(n)%p
            f(i:(i-1)+at%Nlevel) = matmul( gam(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel),x(i:(i-1)+at%Nlevel) )
        !derivatives of all levels to ne.
            !-> compiler errors here write(*,*) matmul( dgcdne(:,:)+dgrdne(:,:), x(1:Nlevels) )
            df(i:(i-1)+at%Nlevel,neq) = matmul( dgcdne(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel), x(i:(i-1)+at%Nlevel) ) + &
                            matmul( dgrdne(i:(i-1)+at%Nlevel,i:(i-1)+at%Nlevel), x(i:(i-1)+at%Nlevel) )
            i = i + at%Nlevel
        enddo
        at => null()

        return
    end subroutine rate_equations

    subroutine particle_conservation (icell, Neq, x, f, df)
    !replace each SEE of each atom with the mass conservation for that atom
        integer, intent(in) :: icell, Neq
        real(kind=dp), intent(in) :: x(neq)
        real(kind=dp), intent(inout) :: f(neq), df(neq, neq)
        real(kind=dp) :: ntotal
        integer :: imaxpop, nlev, i, n

        i = 1
        do n=1, NactiveAtoms
            nlev = activeatoms(n)%p%Nlevel
            ntotal = nHtot(icell) * activeatoms(n)%p%abund
            !relative to sub array of x
            imaxpop = (i-1) + locate(x(i:(i-1)+nlev), maxval(x(i:(i-1)+nlev)))
            f(imaxpop) = 1.0_dp -sum(x(i:(i-1)+nlev)) / ntotal !~0, better to set to 0?
            df(imaxpop,neq) = 0.0_dp
            df(imaxpop,i:(i-1)+nlev) = -1.0_dp / ntotal
            i = i + nlev
        enddo

        return
    end subroutine particle_conservation


    subroutine write_rates()
    !From the rates, the collisional, radiative and total rates matrices can be easily
    ! built.
        integer :: n

        do n=1, NactiveAtoms
            call write_collisional_rates_atom(ActiveAtoms(n)%p)
            !call write_radiative_rates_atom(ActiveAtoms(n)%p)
        enddo

        return
    end subroutine write_rates

    !TO DO
    subroutine write_radiative_rates_atom(atom)
        type(AtomType), intent(inout) :: atom
        integer :: icell, id, kr, i, j, kr_start
        real(kind=dp), dimension(:,:,:), allocatable :: rates
        integer :: naxis, naxes(3), unit, status, nelements
        real(kind=dp) :: nu0

        return
    end subroutine write_radiative_rates_atom

    subroutine write_collisional_rates_atom(atom)
    !change axes order ?
        type(AtomType), intent(inout) :: atom
        integer :: icell, id, kr, i, j, kr_start
        real(kind=dp), dimension(:,:,:), allocatable :: rates
        integer :: naxis, naxes(4), unit, status, nelements
        real(kind=dp) :: nu0

        allocate(rates(atom%Ntr,n_cells,2),stat=status)
        rates(:,:,:) = 0.0_dp
        if (status>0) call error("(write_collisional_rates_atom) Cannot allocate rates!")

        id = 1
        !could be para
        do icell=1, n_cells
            if (icompute_atomRT(icell) > 0) then
                if (atom%id=='H') then
                    call collision_rates_hydrogen_loc(id,icell)
                else
                    call init_colrates_atom(id,atom)
                    call collision_rates_atom_loc(id,icell,atom)
                endif
                do kr=1, atom%Nline
                    rates(kr,icell,1) = atom%lines(kr)%Cij(id)
                    rates(kr,icell,2) = atom%lines(kr)%Cji(id)
                enddo
                do kr=atom%Ntr_line+1, atom%Ntr
                    rates(kr,icell,1) = atom%continua(kr-atom%Ntr_line)%Cij(id)
                    rates(kr,icell,2) = atom%continua(kr-atom%Ntr_line)%Cji(id)
                enddo
            endif
        enddo

        !get unique unit number
        call ftgiou(unit,status)
        if (status > 0) call print_error(status)

        if (lVoronoi) then
            naxis = 2
            naxes(1) = n_cells
            naxes(2) = 2
            nelements = naxes(1) * naxes(2)
        else
            if (l3D) then
                naxis = 4
                naxes(1) = n_rad
                naxes(2) = 2*nz
                naxes(3) = n_az
                naxes(4) = 2
                nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
            else
                naxis = 3
                naxes(1) = n_rad
                naxes(2) = nz
                naxes(3) = 2
                nelements = naxes(1) * naxes(2) * naxes(3)
            end if
        end if


        !separate the rates in individual HDU for more clarity

        call ftinit(unit,trim(atom%ID)//"_colrate.fits.gz",1,status)
        if (status > 0) call print_error(status)
        call ftphpr(unit,.true.,-64,naxis,naxes,0,1,.true.,status)
        if (status > 0) call print_error(status)
        kr = 1
        i = atom%i_trans(kr)
        j = atom%j_trans(kr)
        call ftpkyj(unit, "j", j,'', status)
        call ftpkyj(unit, "i", i, ' ', status)
        if (atom%Nline > 1) then
            nu0 = c_light / nm_to_m / atom%lines(kr)%lambda0
            kr_start = atom%Ntr_line + 1
        else
            nu0 = c_light / nm_to_m / atom%continua(kr)%lambda0
            kr_start = 2
        endif
        call ftpkyd(unit, "nu0", nu0, -14, "Hz", status)
        call ftpprd(unit,1,1,nelements,rates(kr,:,:),status)

        do kr=2, atom%Nline
            call ftcrhd(unit, status)
            if (status > 0) call print_error(status)
            call ftphpr(unit,.true.,-64,naxis,naxes,0,1,.true.,status)
            if (status > 0) call print_error(status)
            i = atom%lines(kr)%i
            j = atom%lines(kr)%j
            call ftpkyj(unit, "j", j,'', status)
            call ftpkyj(unit, "i", i, ' ', status)
            nu0 = c_light / nm_to_m / atom%lines(kr)%lambda0
            call ftpkyd(unit, "nu0", nu0, -14, "Hz", status)
            call ftpprd(unit,1,1,nelements,rates(kr,:,:),status)
        enddo

        do kr=kr_start, atom%Ntr
            call ftcrhd(unit, status)
            if (status > 0) call print_error(status)
            call ftphpr(unit,.true.,-64,naxis,naxes,0,1,.true.,status)
            if (status > 0) call print_error(status)
            i = atom%continua(kr - atom%Ntr_line)%i
            j = atom%continua(kr - atom%Ntr_line)%j
            call ftpkyj(unit, "j", j,'', status)
            call ftpkyj(unit, "i", i, ' ', status)
            nu0 = c_light / nm_to_m / atom%continua(kr-atom%Ntr_line)%lambda0
            call ftpkyd(unit, "nu0", nu0, -14, "Hz", status)
            call ftpprd(unit,1,1,nelements,rates(kr,:,:),status)
        enddo

        call ftclos(unit, status)
        call ftfiou(unit, status)

        if (status > 0) call print_error(status)

        deallocate(rates)
        return
    end subroutine write_collisional_rates_atom



end module see
