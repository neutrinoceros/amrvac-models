!> Customization module
! content :
!   * initial conditions
!   * gravity routine

module mod_usr

    use mod_hd
    use mod_constants

    implicit none
    ! conversion factors
    double precision :: base_length_au = 1d2
    double precision :: au2cm  = 1.49597870691d13 ! a.u. to cm         conversion factor
    double precision :: msun2g = 1.988d33         ! solar mass to gram conversion factor
    double precision :: yr2s   = 3.1536d7         ! year to seconds

    ! &usr_list
    ! "pert" prefix stands for "perturbation"
    double precision :: rhomin, cavity_radius, cavity_width
    integer :: pert_moment = 1  ! number of the moment w(:, pert_moment) to perturb
    double precision :: pert_amp = 0d0  ! relative amplitude of perturbations to initial conditions

    ! &usr_dust_list
    ! custom dust parameters
    double precision :: grain_density_gcm3 = 1d0  ! (g/cm^3)
    double precision :: gas2dust_ratio = 1d2
    double precision, allocatable :: grain_size_cm(:)

contains

   !> This routine should set user methods, and activate the physics module
   subroutine usr_init()
      use mod_disk, only: disk_activate
      use mod_disk_parameters, only: central_mass, ref_radius, rho0, aspect_ratio, temperature_exponent
      use mod_disk_boundaries, only: wave_killing_parabolic

      call read_usr_parameters(par_files)
      ! Choose coordinate system according to user input at setup
      {^IFONED call set_coordinate_system("polar_1.5D")
      if (abs(pert_amp) > 0d0) &
            call mpistop("Error: perturbation is meant for ndim >= 2")
      }
      {^IFTWOD
         call set_coordinate_system("polar_2D")
      }
      {^IFTHREED call mpistop("3D case not implemented")}

      ! A routine for initial conditions is always required
      usr_init_one_grid    => initial_conditions
      usr_set_parameters   => parameters
      usr_process_adv_grid => wave_killing_parabolic

      ! those normalization factors are only used to set up physical grain sizes
      unit_length  = base_length_au * au2cm ! 100au (cm)
      unit_density = msun2g / unit_length**3

      ! devnote: I think this one is completely useless
      ! orbital period at ref_radius (s)
      ! unit_time = yr2s * central_mass**(-0.5) * (base_length_au*ref_radius)**3/2

      call hd_activate()
      call read_usr_dust_parameters(par_files)
      call disk_activate(hd_gravity)
   end subroutine usr_init


   subroutine read_usr_parameters(files)
      !> Read parameters from .par files
      use mod_dust, only: dust_n_species
      use mod_global_parameters, only: unitpar
      character(len=*), intent(in) :: files(:)
      integer n

      namelist /usr_list/ &
           rhomin, cavity_radius, cavity_width, &
           pert_moment, pert_amp


      do n = 1, size(files)
         open(unitpar, file=trim(files(n)), status="old")
         read(unitpar, usr_list, end=111)
111      rewind(unitpar)
      end do
   end subroutine read_usr_parameters

   subroutine read_usr_dust_parameters(files)
      !> Read parameters from .par files
      use mod_dust, only: dust_n_species
      use mod_global_parameters, only: unitpar
      character(len=*), intent(in) :: files(:)
      integer n

      namelist /usr_dust_list/ gas2dust_ratio, grain_density_gcm3, grain_size_cm

      allocate(grain_size_cm(dust_n_species))
      do n = 1, size(files)
         open(unitpar, file=trim(files(n)), status="old")
         read(unitpar, usr_dust_list, end=111)
111      close(unitpar)
      end do
   end subroutine read_usr_dust_parameters


   subroutine parameters()
      ! Overwrite some default parameters.
      use mod_dust, only: dust_n_species, dust_density, dust_size
      use mod_disk_parameters, only: G
      ! .. local ..
      integer :: idust = -1

      ! dust ----------------------------------
      if (hd_dust) then
         do idust = 1, dust_n_species
            !(au2cm)**-1 is 1cm in code units
            dust_size(idust) = grain_size_cm(idust) / unit_length

            !1g/cm^3 in code unit
            dust_density(idust) = grain_density_gcm3 / unit_density
         end do
      end if

      ! note: I'm forcing this value computationnaly to ensure it is always set with maximal precision
      hd_gamma = 5d0/3d0

      if (mype==0) then
         print*,'mod_usr messages ===================================='
         print*, 'Warning : forcing hd_gamma = 5/3'
         print*, 'G/4pi^2 = ', G/(4*dpi**2)
         print*, 'hd_adiab = ', hd_adiab

         if (hd_dust) then
            print*, 'using ', dust_n_species, 'dust bins'
            write(*,'(a,f17.3,a)') ' gas to dust ratio = ', gas2dust_ratio
         end if
         print*,'====================================================='
      end if
   end subroutine parameters


   subroutine initial_conditions(ixI^L, ixO^L, w, x)
      use mod_disk_phys, only: set_keplerian_angular_motion
      use mod_disk_parameters, only: rho_slope, rho0, aspect_ratio, temperature_exponent, central_mass, G, lisoth_eos
      use mod_hd, only: hd_get_pthermal

      use mod_dust, only: dust_n_species, dust_rho, dust_mom, dust_size
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w(ixI^S,1:nw)

      double precision, dimension(ixI^S) :: pth, pressure_term
      double precision :: dust2gas_frac0, sumfrac
      double precision :: partial_dust2gas_fracs(dust_n_species)
      integer :: idust = -1

      ! init all fields to avoid uninitialized values
      w(ixI^S, 1:nw) = 0d0

      ! init density
      w(ixI^S, rho_) = rho0 * x(ixI^S, r_)**rho_slope * 0.5d0 * (1d0 + tanh((x(ixI^S, r_) - cavity_radius) / cavity_width))
      w(ixO^S, rho_) = max(w(ixO^S, rho_), rhomin) ! clip to floor value

      ! set azimuthal moment (rotational equilibrium)
      pressure_term(ixI^S) = 0d0
      pth(ixI^S) = 0d0
      if (hd_energy) call mpistop("Can not use energy equation (not implemented).")

      call hd_get_pthermal(w, x, ixI^L, ixI^L, pth)
      call gradient(pth, ixI^L, ixO^L, r_, pressure_term)
      pressure_term(ixO^S) = pressure_term(ixO^S) * x(ixO^S, r_) / w(ixO^S, rho_)
      w(ixO^S, mom(phi_)) = w(ixO^S, rho_) * dsqrt(G*central_mass / x(ixO^S, r_) + pressure_term(ixO^S))

      ! add perturbations to initial conditions
      if (it == 0 .and. pert_amp > 0.0) then
         call pert_random_noise(ixI^L, ixO^L, w, x, pert_moment, pert_amp)
      end if

      ! dust ----------------------------------
      !     we compute partial dust to gas fractions assuming the total
      !     fraction = 1/gas2dust_ratio and a size distribution of
      !     n(s) \proto s**-3.5 (standard collisional equilibriuum
      !     assumption from debris disks)

      if (hd_dust) then
         dust2gas_frac0 = 1.0d0 / gas2dust_ratio / sum((dust_size(:) / &
         dust_size(1))**(-0.5d0))

         partial_dust2gas_fracs = dust2gas_frac0 * &
         (dust_size(:)/dust_size(1))**(-0.5d0)

         sumfrac = sum(partial_dust2gas_fracs(:))
         if (sumfrac - 1.0d0 / gas2dust_ratio > 1e-14) then
            call mpistop("error in dust init: total dust2gas fraction does not match user parameter")
         end if

         do idust = 1, dust_n_species
            w(ixO^S, dust_rho(idust)) = w(ixO^S, rho_) * partial_dust2gas_fracs(idust)
            w(ixO^S, dust_mom(1, idust)) = 0d0
            call set_keplerian_angular_motion(ixI^L, ixO^L, w, x, dust_rho(idust), dust_mom(phi_, idust))
         end do
      end if
   end subroutine initial_conditions


   subroutine pert_random_noise(ixI^L, ixO^L, w, x, mflag, noise_scale)
      ! random perturbations meant to trigger the RWI.
      integer, intent(in)             :: mflag
      double precision, intent(in)    :: noise_scale
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w(ixI^S,1:nw)

      double precision, dimension(ixI^S) :: random_factor, noise
      integer :: seed_size, block_index, i
      integer, allocatable :: seed(:)

      ! rng seeding. We use 'mype' (the proc index) to avoid repetitions accross processes,
      ! and a calculated (azimuthal) 'block_index' of sorts to further avoid repetitions accross blocks.
      {^IFTWOD
         ! this condition block allows to use variables specific to 2D without breaking compilation in 1D
         block_index = int(x(ixomin1, ixomin2, phi_)/(xprobmax2-xprobmin2) * domain_nx2/block_nx2)
      }

      call random_seed(size=seed_size)
      allocate(seed(seed_size), source=(mype*37+(401*block_index))*[(i, i=0, seed_size-1)])
      call random_seed(put=seed)
      deallocate(seed)

      call hd_get_csound2(w, x, ixI^L, ixO^L, noise)
      noise(ixO^S) = noise_scale * dsqrt(noise(ixO^S)) * w(ixO^S, rho_)

      ! randomize
      call random_number(random_factor(ixO^S))
      noise(ixO^S) = noise(ixO^S) * random_factor(ixO^S)

      ! spatially localize the effect
      noise(ixO^S) = noise(ixO^S) * exp(-(x(ixO^S, r_) - cavity_radius)**2 / (2*cavity_width**2))

      ! application
      w(ixO^S, mom(mflag)) = w(ixO^S, mom(mflag)) + noise(ixO^S) 
   end subroutine pert_random_noise


end module mod_usr
