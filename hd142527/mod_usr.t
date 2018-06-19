!> Customization module
! content :
!   * initial conditions
!   * gravity routine

module mod_usr

  use mod_hd
  use mod_constants

  implicit none
  double precision :: au2cm  = 1.49597870691d13 ! a.u. to cm         conversion factor
  double precision :: msun2g = 1.988d33         ! solar mass to gram conversion factor

  ! &usr_list
  double precision :: density_slope, cavity_radius, cavity_width
  double precision :: rhozero, rhomin
  double precision :: aspect_ratio

  ! &perturbation_list
  logical :: pert_noise = .false.
  integer :: pert_moment = 1
  double precision :: pert_amp = one

  ! &usr_dust_list
  ! custom dust parameters
  double precision :: grain_density_gcm3 = one  ! (g/cm^3)
  double precision :: gas2dust_ratio = 1d2
  double precision, allocatable :: grain_size_cm(:)

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_disk_parameters, only: read_disk_parameters
    use mod_disk_phys, only: central_gravity
    use mod_disk_boundaries, only: wave_killing_parabolic

    ! Choose coordinate system according to user input at setup
    {^IFONED call set_coordinate_system("polar_1.5D")}
    {^IFTWOD call set_coordinate_system("polar_2D")}
    {^IFTHREED call mpistop("3D case not implemented")}

    ! A routine for initial conditions is always required
    usr_init_one_grid    => initial_conditions
    usr_set_parameters   => parameters
    usr_gravity          => central_gravity
    usr_special_bc       => cst_bound
    usr_process_adv_grid => wave_killing_parabolic

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = au2cm                     ! 1au            (cm)
    unit_density = msun2g / au2cm**2         ! 1M_sun / au^2  (g/cm^2)
    unit_time    = 3.1536d7 / (2.0d0*dpi)    ! (2pi)^-1 yr    (s)

    ! Activate the physics module
    call hd_activate()
    call read_disk_parameters(par_files)
    call read_usr_parameters(par_files)
  end subroutine usr_init


  subroutine read_usr_parameters(files)
    !> Read parameters from .par files
    use mod_dust, only: dust_n_species
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer n

    namelist /usr_list/ aspect_ratio,&
         cavity_radius, cavity_width,&
         rhozero, rhomin, density_slope

    namelist /perturbation_list/ pert_noise,&
         pert_moment, pert_amp

    namelist /usr_dust_list/ gas2dust_ratio,&
         grain_density_gcm3, grain_size_cm

    allocate(grain_size_cm(dust_n_species))
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    rewind(unitpar)
       read(unitpar, perturbation_list, end=112)
112    rewind(unitpar)
       read(unitpar, usr_dust_list, end=113)
113    close(unitpar)
    end do
  end subroutine read_usr_parameters


  subroutine parameters()
    ! Overwrite some default parameters.
    use mod_dust, only: dust_n_species, dust_density, dust_size
    use mod_global_parameters
    ! .. local ..
    double precision :: norm_density
    integer i

    hd_adiab = aspect_ratio**2/hd_gamma * rhozero**(one-hd_gamma)
    ! dust ----------------------------------
    norm_density = msun2g/au2cm**3
    if (hd_dust) then
       do i = 1, dust_n_species
          !(au2cm)**-1 is 1cm in code units
          dust_size(i) = grain_size_cm(i) / au2cm

          !1g/cm^3 in code unit
          dust_density(i) = grain_density_gcm3 / norm_density
       end do
    end if

    if (mype==0) then
       print*,'====================================================='
       print*, 'User specific parameters'
       write(*,'(a,f17.3,a)') ' HD_ADIAB = ', hd_adiab
       if (hd_dust) then
          write(*,*), 'using ', dust_n_species, 'dust bins'
          write(*,'(a,f17.3,a)') ' gas to dust ratio = ', gas2dust_ratio
       end if
       print*,'====================================================='
    end if
  end subroutine parameters


  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    ! Set up initial conditions
    use mod_global_parameters
    use mod_disk_phys, only: disk_density_cavity,&
         set_centrifugal_eq_angular_vel,&
         set_keplerian_angular_motion
    use mod_disk_check, only: check_radial_range
    use mod_dust, only: dust_n_species, dust_rho, dust_mom, dust_size
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pressure_term(ixI^S), gradp(ixI^S)
    double precision :: dust2gas_frac0, partial_dust2gas_fracs(dust_n_species)
    integer n

    call check_radial_range(ixI^L, x)
    ! proper init ---------------------------
    w(ixO^S, mom(1)) = zero
    w(ixO^S, mom(2)) = zero
    w(ixO^S, rho_)   = rhomin

    ! density : make a cavity ---------------
    call disk_density_cavity(ixI^L, ixO^L, rhozero, density_slope, cavity_radius, cavity_width, x, w, gradp)
    w(ixO^S, rho_) = max(w(ixO^S, rho_), rhomin) ! clip to floor value
    ! azimuthal velocity --------------------
    pressure_term(ixO^S) = x(ixO^S, r_) * gradp(ixO^S) / w(ixO^S, rho_)
    call set_centrifugal_eq_angular_vel(ixI^L, ixO^L, w, x, pressure_term)

    ! dust ----------------------------------
    ! we compute partial dust to gas fractions assuming the total
    ! fraction = 1/gas2dust_ratio and a size distribution of
    ! n(s) \proto s**-3.5 (standard collisional equilibriuum
    ! assumption from debris disks)

    if (hd_dust) then
       dust2gas_frac0 = one / gas2dust_ratio / sum((dust_size(:)/dust_size(1))**(-0.5))
       partial_dust2gas_fracs = dust2gas_frac0 * (dust_size(:)/dust_size(1))**(-0.5)
       if (sum(partial_dust2gas_fracs(:)) - one / gas2dust_ratio > 1e-14) then
          call mpistop("error in dust init: total dust2gas fraction does not match user parameter")
       end if

       do n=1, dust_n_species
          w(ixO^S, dust_rho(n)) = w(ixO^S, rho_) * partial_dust2gas_fracs(n)
          w(ixO^S, dust_mom(1,n)) = zero
          call set_keplerian_angular_motion(ixI^L, ixO^L, w, x, dust_rho(n), dust_mom(2,n))
       end do
    end if

    ! add perturbations ---------------------
    if (it == 0 .and. pert_noise) then
       call pert_random_noise(ixI^L, ixO^L, w, x, pert_moment, pert_amp)
    end if
  end subroutine initial_conditions


  ! Boundaries
  ! ----------
  subroutine cst_bound(qt,ixG^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_disk_boundaries, only: constant_boundaries
    !embed a function defined in mod_disk_boundaries.t
    integer, intent(in)             :: ixG^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    call constant_boundaries(qt,ixG^L,ixB^L,iB,w,x,rho_)
    call constant_boundaries(qt,ixG^L,ixB^L,iB,w,x,mom(2))
  end subroutine cst_bound

  ! Optional perturbations to the initial state
  ! ------------------------------------------

  subroutine pert_random_noise(ixI^L, ixO^L, w, x, mflag, scale_amp)
    ! random perturbations meant to trigger the RWI.
    use mod_global_parameters
    integer, intent(in)             :: mflag
    double precision, intent(in)    :: scale_amp
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    ! .. local ..
    double precision, dimension(ixI^S) :: amps, norms
    integer :: seed_size, block_index, i
    integer, allocatable :: seed(:)

    {^IFONED call mpistop("Perturbations are meant for ndim >= 2")}

    ! seeding. We use 'mype' (the proc index) to avoid repetitions accross processes,
    ! and a calculated (azimuthal) 'block_index' of sorts to further avoid repetitions accross blocks.
    {^IFTWOD
    block_index = int(x(ixomin1, ixomin2, phi_)/(xprobmax2-xprobmin2) * domain_nx2/block_nx2)
    }!this allow us to use variables specific to 2D without breaking compilation in 1D
    call random_seed(size=seed_size)
    allocate(seed(seed_size), source=(mype*37+(401*block_index))*[(i, i=0, seed_size-1)])
    call random_seed(put=seed)
    deallocate(seed)

    call random_number(amps(ixO^S))

    norms(ixO^S) = sqrt(w(ixO^S, mom(1))**2 + w(ixO^S, mom(2))**2)
    amps(ixO^S) = norms(ixO^S)*(two*amps(ixO^S)-one)*scale_amp

    w(ixO^S, mom(mflag)) = w(ixO^S, mom(mflag)) &
         + amps(ixO^S)*exp(-(x(ixO^S, r_)-cavity_radius)**2/(10*cavity_width**2))
  end subroutine pert_random_noise

end module mod_usr
