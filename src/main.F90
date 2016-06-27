module m_config
  use m_eos
  use m_eos_romenski
  use m_plastic
  use m_plastic_null
  use m_plastic_mises_huber
  use m_ic
  use m_ic_RP
  use m_bc
  use m_bc_offset_periodic
  integer, parameter :: MAX_NAME_LEN=100
  integer, parameter :: MAX_FILENAME_LEN=1000

! maximum characteristic speed in the domain (over the whole
! simulation) This will be highly problem/material dependent: since
! not using an RP based method, should be seen as a configuration
! option, and just an alternative way of setting the (fixed) timestep.
  real max_wave_speed
! CFL coefficient
  real cfl
! timestep length
  real dt
! output every this many steps
  integer outstep
! final time
  real tmax
! whether and when to apply the divergence constraint
  logical div_constraint
  integer div_constraint_step
! to specify the equation of state to use
  character(len=MAX_NAME_LEN) eos_name, plastic_model_name, ic_type, bc_type
  character(len=MAX_FILENAME_LEN) output_filename_stem

  namelist/CONFIG/max_wave_speed,outstep,cfl,tmax,dt, &
       & div_constraint_step,eos_name,plastic_model_name,ic_type,&
       & bc_type,output_filename_stem

! equation of state (determined by eos_name, and set by init_config)
  class(eos), allocatable :: eq
  class(plastic_model), allocatable :: plmodel
  class(ic), allocatable :: initial_conditions
  class(bc), allocatable :: boundary_conditions

contains
  subroutine init_config(u)
    use m_error
    use m_string, only: to_lower
    integer, intent(in) :: u
!   set some default values (or deliberately nonsense values to
!   trigger a sanity check below)

!   only need to specify two of dt, max_wave_speed and cfl
    max_wave_speed = 0
    cfl = 0.9
!   output every timestep by default
    outstep = 1
    tmax = 0
    div_constraint = .true.
!   0 ==> don't apply the divergence constraint 
!   (div_constraint set to .false. further down in this case)
    div_constraint_step = 0
    eos_name = 'none specified'
    plastic_model_name = 'none specified'
    ic_type = 'none specified'
    bc_type = 'none specified'
    output_filename_stem = 'output_'

    read (u, nml=CONFIG)
    rewind(u)

!   sanity checks on configuration
    call assert(max_wave_speed.gt.0, "max_wave_speed must be set >0")
    if (tmax.le.0) call warn("warning: tmax set to <=0")
    call assert(outstep.gt.0, "outstep must be positive")
    if (div_constraint_step.le.0) div_constraint=.false.

!   handle equation of state, plasticity model, initial conditions, boundary conditions
    select case (eos_name)
    case ('Romenski')
       allocate(eos_romenski :: eq)
    case default
       call panic('unknown equation of state requested: ' // eos_name)
    end select
    call eq%init_from_config(u)
    rewind(u)

    select case (to_lower(plastic_model_name))
    case ('mises-huber','mises huber')
       allocate(plastic_model_mises_huber :: plmodel)
    case default
       allocate(plastic_model_null :: plmodel)
    end select
    call plmodel%init_from_config(u)
    rewind(u)

    select case (to_lower(ic_type))
    case ('rp','riemann','riemann problem')
       allocate(ic_RP :: initial_conditions)
    case default
       call panic('unknown initial conditions requested: ' // ic_type)
    end select
    call initial_conditions%init_from_config(u)
    rewind(u)

    select case (to_lower(bc_type))
    case ('offset periodic', 'offset_periodic')
       allocate(bc_offset_periodic :: boundary_conditions)
    case default
       call panic('unknown boundary conditions requested: ' // bc_type)
    end select
    call boundary_conditions%init_from_config(u)
    rewind(u)

  end subroutine init_config
end module m_config

module m_simulation_data
!  use m_config
!  use m_state
!  use m_domain
!  use m_eos
! simulation time
  real t
! timestep counter
  integer it
! solution (conserved and primitive variables), and its transpose (for
! the y sweep).  A separate 'sol_next' array and copy step is needed
! because the x and y sweeps are not independent of the solution in
! the corresponding direction (via L).
  real, allocatable :: sol(:,:,:), solt(:,:,:), sol_next(:,:,:)
  real, allocatable :: solp(:,:,:), solpt(:,:,:)
! second differences of the solution
  real, allocatable :: d(:,:)
! left and right reconstructed solution (precomputed)
  real, allocatable :: uLp(:,:,:), uRp(:,:,:)
  real, allocatable :: uLpt(:,:,:), uRpt(:,:,:)

contains
  ! init_domain MUST be called before init_simulation_data
  subroutine init_simulation_data
    use m_domain, only: nx,ny
    use m_state, only: nq
    allocate(sol(nq,nx,ny), solt(nq,ny,nx), sol_next(nq,nx,ny))
    allocate(solp(nq,nx,ny), solpt(nq,ny,nx))
    allocate(d(nq,max(nx,ny)))
    allocate(uLp(nq,nx,ny), uRp(nq,nx,ny))
    allocate(uLpt(nq,ny,nx), uRpt(nq,ny,nx))
    it = 0; t = 0
  end subroutine init_simulation_data

  subroutine store_prim_solution
    use m_domain, only: nx,ny
    use m_config, only: eq
    use m_eos, only: cons_to_prim
    integer ix, iy
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
    do iy=1,ny; do ix=1,nx
       solp(:,ix,iy) = cons_to_prim(eq, sol(:,ix,iy))
    enddo; enddo
!$OMP END PARALLEL DO
  end subroutine store_prim_solution
end module m_simulation_data

program main
  use m_time
  use omp_lib, only: omp_get_wtime
  use m_state
  use m_ssprk3
  use m_flux
  use m_weno
  use m_materials
  use m_domain
  use m_config
  use m_simulation_data
  use m_L
  use m_poisson

! used to exit the main timestepping loop
  logical stopflag
! wall clock time
  real time

  call init_domain(5)
  call init_config(5)
  call init_simulation_data

  call apply_IC
  call apply_BC
  call rescale_rhoF
  if (div_constraint) call divergence_constraint
  call output

  stopflag = .false.
  call clock(clock_start)
  do 
     dt = (cfl/max_wave_speed) * dx;
     if (t + dt > tmax) then
        dt = tmax - t
        if (dt < 1e-5) exit
        stopflag = .true.
     end if
     select case (mod(it,2))
     case (0)
        call x_sweep
        call apply_BC
        call y_sweep
        call plastic_src
     case (1)
        call y_sweep
        call apply_BC
        call x_sweep
        call plastic_src
     end select
     call apply_BC
     call rescale_rhoF
     if (div_constraint) then
        if (mod(it,div_constraint_step).eq.0) call divergence_constraint
     end if
     t = t + dt
     it = it + 1
     if (stopflag) exit
     if (mod(it,outstep).eq.0) call output
     call clock(time)
     write (0,*) it, t, time-clock_start
  end do

  call rescale_rhoF
  call apply_BC
  if (div_constraint) call divergence_constraint

  call output

contains
  subroutine clock(t)
    real, intent(out) :: t
#ifdef _OPENMP
    t = omp_get_wtime()
#else
    call cpu_time(t)
#endif
  end subroutine clock

  subroutine x_sweep
    integer iy, ix, iq
    call store_prim_solution
    do iy=1, ny
       do ix = 2, nx-1
          d(:,ix) = solp(:,ix+1,iy) - 2*solp(:,ix,iy) + solp(:,ix-1,iy)
       end do
       do ix = 1+2, nx-2
          do iq = 1,nq
             call wenom_reconstruct(solp(iq,ix-2:ix+2,iy),uLp(iq,ix,iy),uRp(iq,ix,iy))
          end do
          uLp(:,ix,iy) = median_state(d(:,ix+1:ix-1:-1), solp(:,ix+1:ix-1:-1,iy), uLp(:,ix,iy))
          uRp(:,ix,iy) = median_state(d(:,ix-1:ix+1),    solp(:,ix-1:ix+1,   iy), uRp(:,ix,iy))
       end do
    end do
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
    do iy=1+2,ny-2
       sol_next(:,:,iy) = advance_solution(eq, sol(:,:,iy-2:iy+2), &
            & solp(:,:,iy-2:iy+2), uLp(:,:,iy-2:iy+2), uRp(:,:,iy-2:iy+2), &
            & dx/dt, dirn=1)
    end do
!$OMP END PARALLEL DO
    sol(:,:,:) = sol_next
  end subroutine x_sweep

  subroutine y_sweep
    integer ix, iy, iq
    call store_prim_solution
    ! transpose the spatial dimensions: the first dimension should
    ! go parallel to the fluxes (in this case y)
    solt(:,:,:)  = reshape(sol,  [nq,ny,nx], order=[1,3,2])
    solpt(:,:,:) = reshape(solp, [nq,ny,nx], order=[1,3,2])
    do ix=1, nx
       do iy = 2, ny-1
          d(:,iy) = solpt(:,iy+1,ix) - 2*solpt(:,iy,ix) + solpt(:,iy-1,ix)
       end do
       do iy = 1+2, ny-2
          do iq = 1,nq
             call wenom_reconstruct(solpt(iq,iy-2:iy+2,ix),uLpt(iq,iy,ix),uRpt(iq,iy,ix))
          end do
          uLpt(:,iy,ix) = median_state(d(:,iy+1:iy-1:-1), solpt(:,iy+1:iy-1:-1,ix), uLpt(:,iy,ix))
          uRpt(:,iy,ix) = median_state(d(:,iy-1:iy+1),    solpt(:,iy-1:iy+1,   ix), uRpt(:,iy,ix))
       end do
    end do
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
    do ix=1+2,nx-2
        sol_next(:,ix,:) = advance_solution(eq, solt(:,:,ix-2:ix+2), &
             & solpt(:,:,ix-2:ix+2), uLpt(:,:,ix-2:ix+2), uRpt(:,:,ix-2:ix+2), &
             & dx/dt, dirn=2)
    end do
!$OMP END PARALLEL DO
    sol(:,:,:) = sol_next(:,:,:)
  end subroutine y_sweep

  ! loop over each cell, relax each one to the yield surface
  subroutine plastic_src
    real F(3,3), rhoF(3,3)
    integer ix,iy
    call store_prim_solution
    do ix=1,nx
       do iy=1,ny
          F = reshape(solp(prim_F,ix,iy),[3,3])
          call plastic_relax(eq, plmodel, solp(prim_S,ix,iy), F)
          ! prim_F of solp now stores the relaxed deformation gradient
          ! (the prim_S component is unchanged)
          solp(prim_F,ix,iy) = reshape(F,[9])
          sol(:,ix,iy) = prim_to_cons(eq,solp(:,ix,iy))
       end do
    end do
  end subroutine plastic_src

  subroutine apply_IC
    use m_ic
    call initial_conditions%apply(eq,sol)
    sol_next(:,:,:) = sol(:,:,:)
  end subroutine apply_IC

  ! boundary conditions
  subroutine apply_BC
    use m_bc
    call boundary_conditions%apply(eq,sol,nb)
    !call offset_periodic(sol,nb)
  end subroutine apply_BC

  ! rescale rhoF so that rho0*det(F) is relaxed to the conserved density
  subroutine rescale_rhoF
    use m_matutil, only: det3
    real scale_factor, target_det, actual_det
    integer ix,iy
    do iy=1,ny; do ix=1,nx
       target_det = eq%rho0 * sol(cons_rho,ix,iy)**2
       actual_det = det3(cons_get_rhoF(sol(:,ix,iy)))
       ! Miller-Collela (2001): stability of the rescaling depends on
       ! a relaxation time of six timesteps
       scale_factor = 1 + ((target_det/actual_det)**(1.0/3.0) - 1)/6.0
       ! remember: reducing detF increases the density, so if the
       ! density computed from rhoF is too small, want to scale rhoF
       ! down.
       sol(cons_rhoF,ix,iy) = scale_factor * sol(cons_rhoF,ix,iy)
    end do; end do
  end subroutine rescale_rhoF

  ! applies the constraint \partial_k (\rho F_kj) = 0
  subroutine divergence_constraint
    use m_matutil, only: div_2d, grad_2d
    ! first component is 'j' in the above expression (i.e., remaining
    ! index on F)
    integer i,j,ix,iy
    real div_rhoF(3,2:nx-1,2:ny-1), rhoFv(3,3,-1:1,-1:1), rhoF(3,3)
    real phi(3,2:nx-1,2:ny-1), rhoF_corrected(3,3)

    div_rhoF = 0

    ! compute div (rho F), interpreted as in the comment above
    do iy = 3,ny-2
       do ix = 3,nx-2
          do j=-1,1
             do i=-1,1
                rhoFv(:,:,i,j) = cons_get_rhoF(sol(:,ix+i,iy+j))
             end do
          end do

          do j=1,3
             div_rhoF(j,ix,iy) = div_2d(rhoFv(:,j,:,:))/dx
          end do
       end do
    end do

    do j=1,3
       call poisson_2d(div_rhoF(j,:,:), phi(j,:,:), dx)
    end do

    ! size of boundary here?
    do iy = 5,ny-4
       do ix = 5,nx-4
          rhoF = cons_get_rhoF(sol(:,ix,iy))
          do j=1,3
             rhoF_corrected(:,j) = rhoF(:,j) + [grad_2d(phi(j,ix-1:ix+1,iy-1:iy+1)), 0.0]/dx
          end do
          sol(cons_rhoF,ix,iy) = reshape(rhoF_corrected,[9])
       end do
    end do
  end subroutine divergence_constraint

  subroutine output
    use m_output
    character(len=MAX_FILENAME_LEN) fname
    call store_prim_solution
    write (fname,'(A,I0,A)') trim(output_filename_stem), it, '.vtk'
    open(unit=16,file=fname)
    call visit_output(16,it,t,eq,sol,solp)
    close(16)
  end subroutine output
end program main
