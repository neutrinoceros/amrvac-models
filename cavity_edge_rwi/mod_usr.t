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
  double precision :: unit_mass

  ! &usr_list
  double precision :: rhomin, cavity_radius, cavity_width
  logical :: constant_pressure = .false. !useful for debuging

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
      use mod_disk, only: disk_activate
      use mod_disk_parameters, only: central_mass, ref_radius, rho0, aspect_ratio, temperature_exponent
      use mod_disk_boundaries, only: wave_killing_parabolic

      call read_usr_parameters(par_files)
      ! Choose coordinate system according to user input at setup
      {^IFONED call set_coordinate_system("polar_1.5D")
      if (pert_noise) &
            call mpistop("Error: pert_noise=.true. is meant for ndim >= 2")
      }
      {^IFTWOD
         call set_coordinate_system("polar_2D")
      }
      {^IFTHREED call mpistop("3D case not implemented")}

      ! A routine for initial conditions is always required
      usr_init_one_grid    => initial_conditions
      usr_set_parameters   => parameters
      usr_process_adv_grid => wave_killing_parabolic
      !usr_special_bc => cst_bound (unused, deprecated)

      ! Choose independent normalization units if using dimensionless variables.
      unit_mass    = msun2g ! NOT A STANDARD AMRVAC VARIABLE
      unit_length  = base_length_au * au2cm ! 100au (cm)
      unit_density = unit_mass / unit_length**2

      ! orbital period at ref_radius (s)
      unit_time = yr2s * central_mass**(-0.5) * (base_length_au*ref_radius)**3/2

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

      namelist /usr_list/ rhomin, cavity_radius, cavity_width, constant_pressure

      namelist /perturbation_list/ pert_noise, pert_moment, pert_amp

      do n = 1, size(files)
         open(unitpar, file=trim(files(n)), status="old")
         read(unitpar, usr_list, end=111)
111      rewind(unitpar)
         read(unitpar, perturbation_list, end=112)
112      close(unitpar)
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
      double precision :: norm_density
      integer i

      ! dust ----------------------------------
      norm_density = unit_mass / unit_length**3
      if (hd_dust) then
         do i = 1, dust_n_species
            !(au2cm)**-1 is 1cm in code units
            dust_size(i) = grain_size_cm(i) / unit_length

            !1g/cm^3 in code unit
            dust_density(i) = grain_density_gcm3 / norm_density
         end do
      end if

      ! note: I'm forcing this value computationnaly to ensure it is always set with maximal precision
      hd_gamma = 5d0/3d0

      if (mype==0) then
         print*,'mod_user messages ======================================='
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

      double precision, dimension(ixI^S) :: gradp_r, pth, pressure_term, tanh_term
      double precision :: dust2gas_frac0, sumfrac
      double precision :: partial_dust2gas_fracs(dust_n_species)
      integer :: idust  = -1

      if (hd_energy) call mpistop("Can not use energy equation (not implemented).")

      ! proper init ---------------------------
      pressure_term(ixI^S) = 0d0
      pth(ixI^S) = 0d0
      gradp_r(ixI^S) = 0d0

      w(ixI^S, 1:nw) = 0.0d0
      w(ixI^S, rho_) = rho0 * x(ixI^S, r_)**rho_slope * 0.5d0 * (1d0 + tanh((x(ixI^S, r_) - cavity_radius) / cavity_width))
      w(ixO^S, rho_) = max(w(ixO^S, rho_), rhomin) ! clip to floor value

      ! Set rotational equilibrium

      ! analytic pressure gradient (for barotropic case only, deprecated)
      if (lisoth_eos) then
         call hd_get_pthermal(w, x, ixI^L, ixI^L, pth)
         call gradient(pth, ixI^L, ixO^L, r_, gradp_r)
      else
         tanh_term(ixO^S) = tanh(((x(ixO^S, r_) - cavity_radius) / cavity_width))
         gradp_r(ixO^S) = 0.5d0*rho0*(rho_slope* x(ixO^S, r_)**(rho_slope-1.0d0) * (1.0d0 + tanh_term(ixO^S)) + x(ixO^S, r_)**(rho_slope) * (1.0d0 - tanh_term(ixO^S)**2) / cavity_width)
         gradp_r(ixO^S) = hd_adiab * hd_gamma * w(ixO^S, rho_)**(hd_gamma-1.0d0) * gradp_r(ixO^S)
      end if
      pressure_term(ixO^S) = x(ixO^S, r_) * gradp_r(ixO^S) / w(ixO^S, rho_)

      ! azimuthal velocity --------------------
      w(ixO^S, mom(phi_)) = w(ixO^S, rho_) * dsqrt(G*central_mass / x(ixO^S,r_) + pressure_term(ixO^S))

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

      ! add perturbations ---------------------
      if (it == 0 .and. pert_noise) then
         call pert_random_noise(ixI^L, ixO^L, w, x, pert_moment, pert_amp)
      end if
   end subroutine initial_conditions


!     Optional perturbations to the initial state
!     ------------------------------------------

   subroutine pert_random_noise(ixI^L, ixO^L, w, x, mflag, scale_amp)
!     random perturbations meant to trigger the RWI.
      integer, intent(in)             :: mflag
      double precision, intent(in)    :: scale_amp
      integer, intent(in)             :: ixI^L, ixO^L
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w(ixI^S,1:nw)
!     .. local ..
      double precision, dimension(ixI^S) :: amps, norms
      integer :: seed_size, block_index, i
      integer, allocatable :: seed(:)

!     seeding. We use 'mype' (the proc index) to avoid repetitions accross processes,
!     and a calculated (azimuthal) 'block_index' of sorts to further avoid repetitions accross blocks.
      {^IFTWOD
      block_index = int(x(ixomin1, ixomin2, phi_)/(xprobmax2-xprobmin2)&
      * domain_nx2/block_nx2)
      }
!     this condition block allows to use variables specific to 2D without breaking compilation in 1D
      call random_seed(size=seed_size)
      allocate(seed(seed_size), source=(mype*37+(401*block_index))*[(i, i=0, seed_size-1)])
      call random_seed(put=seed)
      deallocate(seed)

      call random_number(amps(ixO^S))

      norms(ixO^S) = sqrt(w(ixO^S, mom(r_))**2 + w(ixO^S, mom(phi_))**2)
      amps(ixO^S) = norms(ixO^S) * (2d0*amps(ixO^S)-1d0) * scale_amp

      w(ixO^S, mom(mflag)) = w(ixO^S, mom(mflag)) + amps(ixO^S) * exp(-(x(ixO^S, r_) - cavity_radius)**2 / (10*cavity_width**2))
   end subroutine pert_random_noise

   !  unused (to be removed)
   subroutine cst_bound(qt,ixG^L,ixB^L,iB,w,x)
      use mod_disk_boundaries, only: constant_boundaries
      integer, intent(in)             :: ixG^L, ixB^L, iB
      double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
      double precision, intent(inout) :: w(ixG^S,1:nw)
      call constant_boundaries(qt,ixG^L,ixB^L,iB,w,x,rho_)
      call constant_boundaries(qt,ixG^L,ixB^L,iB,w,x,mom(phi_))
   end subroutine cst_bound

end module mod_usr
