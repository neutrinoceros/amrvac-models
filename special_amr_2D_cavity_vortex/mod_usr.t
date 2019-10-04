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
  character(len=std_len) :: usr_geometry
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
    use mod_global_parameters
    use mod_usr_methods
    use mod_disk, only: disk_activate
    use mod_disk_parameters, only: central_mass, ref_radius
    use mod_disk_boundaries, only: wave_killing_parabolic

    call read_usr_parameters(par_files)
    ! Choose coordinate system according to user input at setup
    {^IFONED call set_coordinate_system("polar_1.5D")
    if (pert_noise) &
         call mpistop("Error: pert_noise=.true. is meant for ndim >= 2")
    }
    {^IFTWOD
    select case(trim(usr_geometry))
    case('rphi')
       call set_coordinate_system("polar_2D")
    case('rz')
       if (pert_noise) &
            call mpistop("Error: pert_noise=.true. is not compatible with usr_geometry='rz'")
       call set_coordinate_system("cylindrical_2.5D")
    case default
       call mpistop("Error: usr_geometry is not set. Choose 'rz' or 'rphi'.")
    end select
    }
    {^IFTHREED call mpistop("3D case not implemented")}! set "cylindrical_3D" here

    ! A routine for initial conditions is always required
    usr_init_one_grid    => initial_conditions
    usr_set_parameters   => parameters
    usr_special_bc       => cst_bound
    usr_process_adv_grid => wave_killing_parabolic

    usr_refine_grid     => specialrefine_grid
    usr_var_for_errest  => curl_for_errest

    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output


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

    namelist /usr_list/ usr_geometry, rhomin,&
         cavity_radius, cavity_width,&
         constant_pressure

    namelist /perturbation_list/ pert_noise,&
         pert_moment, pert_amp

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_list, end=111)
111    rewind(unitpar)
       read(unitpar, perturbation_list, end=112)
112    close(unitpar)
    end do
  end subroutine read_usr_parameters

  subroutine read_usr_dust_parameters(files)
    !> Read parameters from .par files
    use mod_dust, only: dust_n_species
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer n

    namelist /usr_dust_list/ gas2dust_ratio,&
         grain_density_gcm3, grain_size_cm

    allocate(grain_size_cm(dust_n_species))
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, usr_dust_list, end=111)
111    close(unitpar)
    end do
  end subroutine read_usr_dust_parameters


  subroutine parameters()
    ! Overwrite some default parameters.
    use mod_dust, only: dust_n_species, dust_density, dust_size
    use mod_disk_parameters, only: G
    use mod_global_parameters
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

    if (mype==0) then
       print*,'User messages ======================================='
       write(*,*), 'G/4pi^2 = ', G/(4*dpi**2)
       write(*,*), 'hd_adiab = ', hd_adiab
       write(*,*), 'usr_geometry = ', usr_geometry

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
    use mod_disk_phys, only: set_keplerian_angular_motion
    use mod_disk_parameters, only: rho_slope, rho0, aspect_ratio, central_mass, G
    use mod_dust, only: dust_n_species, dust_rho, dust_mom, dust_size
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: tanh_term(ixI^S), gradp_r(ixI^S)
    double precision :: dust2gas_frac0, partial_dust2gas_fracs(dust_n_species)
    integer n

    if (hd_energy) &
       call mpistop("Can not use energy equation (case not implemented).")

    ! proper init ---------------------------
    w(ixO^S, 1:nw) = 0.0d0

    w(ixO^S, rho_) = rho0 * x(ixO^S, r_)**rho_slope * 0.5d0 * (1d0 + tanh((x(ixO^S, r_) - cavity_radius) / cavity_width))
    if (z_ > 0) & !vertical hydrostatic equilibrium
         w(ixO^S, rho_) =  w(ixO^S, rho_) * exp(-x(ixO^S, z_)**2 / (2d0*aspect_ratio * x(ixO^S, r_))**2)

    w(ixO^S, rho_) = max(w(ixO^S, rho_), rhomin) ! clip to floor value

    ! Set rotational equilibrium
    tanh_term(ixO^S) = tanh(((x(ixO^S, r_) - cavity_radius) / cavity_width))
    ! compute the radial component of grad(P) (from grad(rho))
    gradp_r(ixO^S) = rho0 * 0.5d0 *(rho_slope * x(ixO^S, r_)**(rho_slope-1.0d0) * (1.0d0 + tanh_term(ixO^S)) &
                                 + x(ixO^S, r_)**(rho_slope) * (1.0d0 - tanh_term(ixO^S)**2) / cavity_width)

    if (z_ > 0) &
         gradp_r(ixO^S) = gradp_r(ixO^S) * exp(-x(ixO^S, z_)**2 / (2d0*aspect_ratio * x(ixO^S, r_))**2)

    gradp_r(ixO^S) = hd_adiab * hd_gamma * w(ixO^S, rho_)**(hd_gamma-1.0d0) * gradp_r(ixO^S)
    if (constant_pressure) gradp_r(ixO^S) = 0.0d0 !dbg

    if (z_ > 0) then
       w(ixO^S, mom(phi_)) = w(ixO^S, rho_) * dsqrt( &
            + G*central_mass * x(ixO^S, r_)**2 / (x(ixO^S, r_)**2 + x(ixO^S, z_)**2)**1.5d0 &
            + x(ixO^S, r_) / w(ixO^S, rho_) * gradp_r(ixO^S))
    else
       w(ixO^S, mom(phi_)) = w(ixO^S, rho_) * dsqrt( &
            + G*central_mass / x(ixO^S, r_) &
            + x(ixO^S, r_) / w(ixO^S, rho_) * gradp_r(ixO^S))
    end if

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
          w(ixO^S, dust_mom(1,n)) = 0d0
          call set_keplerian_angular_motion(ixI^L, ixO^L, w, x, dust_rho(n), dust_mom(phi_, n))
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
    call constant_boundaries(qt,ixG^L,ixB^L,iB,w,x,mom(phi_))
    !if (z_ > 0) &
    !     call constant_boundaries(qt,ixG^L,ixB^L,iB,w,x,mom(z_))
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

    ! seeding. We use 'mype' (the proc index) to avoid repetitions accross processes,
    ! and a calculated (azimuthal) 'block_index' of sorts to further avoid repetitions accross blocks.
    {^IFTWOD
    block_index = int(x(ixomin1, ixomin2, phi_)/(xprobmax2-xprobmin2) * domain_nx2/block_nx2)
    }!this condition block allows to use variables specific to 2D without breaking compilation in 1D
    call random_seed(size=seed_size)
    allocate(seed(seed_size), source=(mype*37+(401*block_index))*[(i, i=0, seed_size-1)])
    call random_seed(put=seed)
    deallocate(seed)

    call random_number(amps(ixO^S))

    norms(ixO^S) = sqrt(w(ixO^S, mom(r_))**2 + w(ixO^S, mom(phi_))**2)
    amps(ixO^S) = norms(ixO^S) * (2d0*amps(ixO^S)-1d0) * scale_amp

    w(ixO^S, mom(mflag)) = w(ixO^S, mom(mflag)) &
         + amps(ixO^S) * exp(-(x(ixO^S, r_) - cavity_radius)**2 / (10*cavity_width**2))
  end subroutine pert_random_noise

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: wlocal(ixI^S,1:nw)
    double precision :: vvec(ixI^S,1:ndir)
    ! For ndir=2 only 3rd component of curl can exist
    double precision :: vorticity(ixI^S,7-2*ndir:3)
    integer          :: idirmin,idirmin0,idir

    idirmin0 = 7-2*ndir

    ! wlocal only to get to vector of length 1:nw, which can use precoded subroutines
    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! compute velocity
    do idir=1,ndir
       vvec(ixI^S,idir)=wlocal(ixI^S,mom(idir))/wlocal(ixI^S,rho_)
    enddo

    call curlvector(vvec,ixI^L,ixO^L,vorticity,idirmin,idirmin0,ndir)
    w(ixO^S,nw+1)=vorticity(ixO^S,3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the varnames/primnames string
  character(len=*) :: varnames
  varnames='vortz'

  end subroutine specialvarnames_output

  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    if (any((x(ix^S,1) <= 2.0d0) .and. (x(ix^S,1)>= 1.5d0))) then
        refine=1
        coarsen=-1
    endif

  end subroutine specialrefine_grid

  subroutine curl_for_errest(ixI^L,ixO^L,iflag,w,x,var)
    integer, intent(in)           :: ixI^L,ixO^L,iflag
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: var(ixI^S)

    double precision :: vvec(ixI^S,1:ndir)
    ! For ndir=2 only 3rd component of curl can exist
    double precision :: vorticity(ixI^S,7-2*ndir:3)
    integer          :: idirmin,idirmin0,idir

    ! compute velocity
    do idir=1,ndir
       vvec(ixI^S,idir)=w(ixI^S,mom(idir))/w(ixI^S,rho_)
    enddo

    !call curlvector(vvec,ixI^L,ixO^L,vorticity,idirmin,idirmin0,ndir)
    !var(ixO^S)=vorticity(ixO^S,3)
    var(ixO^S)=vvec(ixO^S,1)
   
  end subroutine curl_for_errest

end module mod_usr
