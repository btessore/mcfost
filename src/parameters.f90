module parametres
!****************************
! Parametres generaux
!****************************

  use mcfost_env, only : dp

  implicit none
  save

  real :: para_version

  logical :: lpara, lstop_after_init
  integer :: indice_etape, etape_i, etape_f

  ! Nombre de photons lances
  logical :: ldust_transfer
  integer :: nbre_photons_loop, nbre_photons_eq_th, nbre_photons_lambda, nbre_photons_image, nbre_photons_spectre
  real :: nbre_photons_lim = 1.e4 ! combien de fois plus on aurait recu sans disque
  integer :: nnfot1
  real :: nbre_photons_tot
  real(kind=dp) :: E_paquet
  integer :: n_dif_max_eq_th = 100000 ! Nbre max de dif autorises dans calcul eq. th OUTDATED
  real :: tau_dark_zone_eq_th = 1500 !1500.   15000 pour benchmark a tau=1e6
  real :: tau_dark_zone_obs = 100 ! idem que 1000 si 1500. ci-dessus
  integer :: n_Stokes

  ! Nbre d'angles pour echantillonner la fct de phase
  integer, parameter :: nang_scatt = 180  ! TODO : ca bug si c'est pas 180

  ! lvariable_dust = true si les proprites de la poussiere sont variables dans chaque cellule
  logical :: lvariable_dust, lmigration, lhydrostatic, ldust_sublimation
  integer :: settling_type ! 1 = Parametric, 2 = Dubrulle or 3 = Fromang

  logical :: lRE_LTE, lRE_nLTE, lnRE, lonly_LTE, lonly_nLTE
  logical :: loutput_J, loutput_J_step1, loutput_UV_field, lxJ_abs, lxJ_abs_step1

  ! Methode de calcul de la diffusion : a choisir pour optimiser taille memoire et temps cpu
  ! 0 -> automatique
  ! 1 -> choix taille du grain diffuseur + matrice Mueller par grain
  ! 2 -> matrice de Mueller moyenne par cellule (benchmark)
  integer :: scattering_method, scattering_method0
  logical :: lscattering_method1

  ! Theorie de Mie ou HG
  integer :: aniso_method ! 1 = full phase function, 2 = HG
  logical :: lmethod_aniso1

  integer :: RT_sed_method ! cf routine dust_map pour def

  ! Etapes de l'�mission thermique
  logical :: ltemp, lsed, lsed_complete, l_em_disk_image, lchauff_int, lextra_heating, lno_internal_energy
  character(len=512), dimension(:), allocatable :: indices
  character(len=512) :: tab_wavelength

  ! Emission moleculaire
  logical :: lemission_mol,  lpop, lprecise_pop, lmol_LTE, ldust_mol, lonly_top, lonly_bottom

  ! Atomic line radiative transfer
  !!lstore_opac futur deprecation
  logical :: lemission_atom, lelectron_scattering, lvacuum_to_air, lcontrib_function_ray, &
  				lmagnetoaccr, lforce_lte, lspherical_velocity, lstop_after_jnu, &
       			ldissolve, laccurate_integ, loutput_rates, lorigin_atom, lzeeman_polarisation!lcontrib_function
  integer :: Nrays_atom_transfer, istep_start
  
  !HEALpix
  integer :: healpix_lorder, healpix_lmin, healpix_lmax !lmin and lmax not yet (for local evaluation)
  
  logical :: llimit_mem, lfix_backgrnd_opac
  logical :: lcheckpoint, lsafe_stop
  !Convergence relative errors
  real :: dpops_max_error, dpops_sub_max_error, art_hv, safe_stop_time
  integer :: checkpoint_period
  
  !Ng's acceleration
  logical :: lng_acceleration
  integer :: iNg_Norder, iNg_ndelay, iNg_Nperiod
  
  !electron density
  logical :: lsolve_for_ne, lno_iterate_ne_mc = .true. !.true.==no electron iteration during step 2
  integer :: ndelay_iterate_ne, n_iterate_ne !0 means once SEE is solved. Otherwise, > 1, iterated every n_iterate_ne during the nlte_loop
  
  !Wavelength table for spectrally resolved images and spectra
  character(len=50) :: tab_wavelength_image, jnu_atom_file
  logical :: ltab_wavelength_image, lread_jnu_atom
  !lmhd_voronoi enables the voronoi tesselation for either a pluto-formatted model (lpluto_file)
  !or an ascii file ready for tesselation (lmodel_ascii)
  integer :: n_pluto_files
  logical :: lpluto_file, lmodel_ascii, lmhd_voronoi

  ! Decomposition image
  logical :: lsepar_contrib, lsepar_pola, lonly_capt_interet
  integer :: N_type_flux
  ! les flux sont I, (Q,U,V), (star, scatt, disk th, disk th scatt.)

  ! rotation du plan du disque en deg., sens trigo.
  real(kind=dp) :: ang_disque

  ! Production d'images symetriques
  ! La symetrie est effectuee avant de choisir les pixels
  ! le syst�me est-il centrosymetrique
  ! le systeme a-t-il une symetrie axiale (ne compte que si N_phi > 1)
  logical :: l_sym_ima, l_sym_centrale, l_sym_axiale

  ! Parametres des cartes
  integer :: N_thet, N_incl, N_phi, capt_interet, delta_capt, capt_inf, capt_sup, capt_debut, capt_fin
  integer ::  npix_x, npix_y, npix_x_save, npix_y_save
  real :: angle_interet, zoom, tau_seuil, wl_seuil

  real  :: cutoff = 7.0
  
  !must be initialized to 0
  integer(kind=8) :: mem_alloc_tot !total memory allocated dynamically or not in bytes

  ! R�solution de la grille de densit�
  ! Nombre de cellules dans la direction r (echantillonage log)
  integer :: grid_type ! 1 = cylindrical, 2 = spherical
  integer :: n_rad, n_rad_in  ! subdivision de la premiere cellule
  ! Nombre de couches verticales ( = de stratifications)
  integer :: nz, p_n_rad, p_nz, p_n_az, p_n_lambda_pos, p_n_lambda_grain
  ! Nombre de cellules azimuthales
  integer :: n_az, j_start, pj_start
  ! Nombre de cellules totale
  integer :: n_cells, nrz, p_n_cells, icell_ref

  logical :: letape_th, limg, lorigine, laggregate, l3D, lremove, lwarp, lcavity, ltilt, lwall
  logical :: lopacite_only, lseed, ldust_prop, ldisk_struct, loptical_depth_map, lreemission_stats
  logical :: lapprox_diffusion, lcylindrical, lspherical, lVoronoi, is_there_disk, lno_backup, lonly_diff_approx, lforce_diff_approx
  logical :: laverage_grain_size, lisotropic, lno_scattering, lqsca_equal_qabs
  logical :: ldensity_file, lsigma_file, lvelocity_file, lvfield_cyl_coord, lphantom_file, lphantom_multi, lphantom_avg
  logical :: lgadget2_file, lascii_SPH_file, llimits_file, lforce_SPH_amin, lforce_SPH_amax, lmcfost_lib
  logical :: lweight_emission, lcorrect_density, lProDiMo2mcfost, lProDiMo2mcfost_test, lastrochem, lML
  logical :: lspot, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium, lchange_Tmax_PAH, lISM_heating, lcasa
  integer :: ISR_model ! 0 : no ISM radiation field, 1 : ProDiMo, 2 : Bate & Keto

  ! benchmarks
  logical :: lbenchmark_Pascucci, lbenchmark_vanZadelhoff1, lbenchmark_vanZadelhoff2, lDutrey94, lHH30mol
  logical :: lbenchmark_water1, lbenchmark_water2, lbenchmark_water3, lbenchmark_SHG, lMathis_field
  real :: Mathis_field

  ! Prodimo
  logical :: lprodimo, lprodimo_input_dir, lforce_ProDiMo_PAH

  logical, parameter :: ltest_rt3 = .false. ! marche pas
  logical, parameter :: ltest_rt4 = .false.  ! marche pas non plus

  logical :: lSeb_Charnoz, lread_Seb_Charnoz, lread_Seb_Charnoz2, lread_Misselt, lread_DustEM
  logical :: lread_grain_size_distrib, lphase_function_file,ltau1_surface, lwrite_column_density, lwrite_mol_column_density

  ! Phantom
  logical :: ldudt_implicit, lscale_length_units, lscale_mass_units, lignore_dust
  logical :: ldelete_Hill_sphere, lrandomize_Voronoi, lrandomize_azimuth, lrandomize_gap, lrandomize_outside_gap, lcentre_on_sink
  real(kind=dp) :: ufac_implicit,scale_length_units_factor,scale_mass_units_factor,correct_density_factor_elongated_cells
  real(kind=dp) :: SPH_amin, SPH_amax, fluffyness, gap_factor
  logical :: lupdate_velocities, lno_vr, lno_vz, lvphi_Kep, lfluffy
  integer :: isink_centre

  ! Disk parameters
  real :: distance ! Distance du disque en pc
  real(kind=dp) :: map_size

  integer :: n_zones, n_regions

  type disk_zone_type
     real(kind=dp) :: Rin, Rmin, Rc, Rout, Rmax, Rref, edge, exp_beta, surf
     real(kind=dp) :: moins_gamma_exp, sclht, diskmass, gas_to_dust, vert_exponent
     integer :: geometry ! 1=disk, 2=tappered-disk, 3=envelope
     integer :: region
  end type disk_zone_type

  type disk_region_type
     integer :: n_zones
     real(kind=dp) :: Rmin, Rmax
     integer :: iRmin, iRmax
     integer, dimension(:), allocatable :: zones
  end type disk_region_type

  type cavity_type
     real(kind=dp) ::  exp_beta, sclht, Rref
  end type cavity_type

  type(disk_zone_type), dimension(:), allocatable, target :: disk_zone
  type(disk_region_type), dimension(:), allocatable, target :: regions
  type(cavity_type) :: cavity

  real(kind=dp) :: Rmin, Rmax, Rmax2, diskmass, correct_Rsub

  real :: exp_strat, a_strat
  real :: alpha

  ! Description analytique du puffed-up inner rim
  logical :: lpuffed_rim
  real :: puffed_rim_h, puffed_rim_r, puffed_rim_delta_r

  real :: z_warp, tilt_angle

  ! SPH
  real :: SPH_keep_particles, planet_az
  logical :: lplanet_az, lfix_star, lcorrect_density_elongated_cells, lturn_off_planets, lturn_off_Lacc, lforce_Mdot
  integer :: which_planet

  logical :: lgap_Gaussian
  real :: f_gap_Gaussian, r_gap_Gaussian, sigma_gap_Gaussian

  ! Correction locale de la densite (dans un anneau)
  real :: correct_density_factor, correct_density_Rin, correct_density_Rout

  ! Vertical scaling of the envelope
  real :: z_scaling_env

  character(len=512) :: density_file, sigma_file, grain_size_file, limits_file
  character(len=512), dimension(:), allocatable :: density_files
  integer :: n_phantom_files

  !Need to know the position on the star of the spot and its extent
  type surface_brightness_type
     real(kind=dp) :: r(3), muo, mui, phio, phii
     real(kind=dp) :: T ! Temperature of the region
  end type surface_brightness_type

  ! Stars
  ! 27/02/2019, adding spots coordinates
  type star_type
     ! todo : indicate all units
     real :: T, M, fUV, slope_UV, othin_sublimation_radius
     real :: Mdot ! Msun
     real(kind=dp) :: r, x,y,z, vx,vy,vz
     logical :: lb_body, out_model, find_spectrum, force_Mdot
     character(len=512) :: spectre
     integer :: icell
     integer :: Nr = 0
     type(surface_brightness_type), dimension(:), allocatable :: SurfB
  end type star_type

  integer :: n_etoiles
  type(star_type), dimension(:), allocatable :: etoile

  ! Spot
  real :: T_spot, surf_fraction_spot, theta_spot, phi_spot

  ! Limb darkening
  logical :: llimb_darkening
  character(len=512) :: limb_darkening_file
  real, dimension(:), allocatable :: mu_limb_darkening, limb_darkening, pola_limb_darkening

  ! ISM radiation following ProDiMo's definitions
  real, parameter :: Wdil =  9.85357e-17
  real, parameter :: TCmb = 2.73
  real, parameter :: T_ISM_stars = 20000.
  real :: chi_ISM = 1.0

  character(len=8) :: system_age

  logical :: lscatt_ray_tracing, lscatt_ray_tracing1, lscatt_ray_tracing2, loutput_mc
  ! ray-tracing 1 : sauve le champ de radiation diffuse
  ! ray-tracing 2 : sauve l'intensite specifique

  integer :: n_pop


  real :: T_max, T_min, Tmax_PAH ! Temp_sublimation et Temp nuage
  integer  :: n_T

  real, dimension(:), allocatable :: s11_file

  ! inclinaisons
  real :: RT_imin, RT_imax, RT_az_min, RT_az_max
  integer ::  RT_n_incl, RT_n_az
  logical :: lRT_i_centered


  integer :: nLevels
  real(kind=dp) :: largeur_profile


end module parametres
