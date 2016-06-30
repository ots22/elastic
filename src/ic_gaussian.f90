! Gaussian pulse initial conditions
module m_ic_gaussian
  use m_ic

  private
  
  type, public, extends(ic) :: ic_gaussian
     real stddev
   contains
     procedure apply
     procedure init_from_config
  end type ic_gaussian

contains
  subroutine apply(this,eq,u)
    use m_state
    use m_eos, only: eos
    use m_domain
    use m_matutil, only: pi, identity

    class(ic_gaussian) :: this
    class(eos), intent(in) :: eq
    real, intent(out) :: u(:,:,:)
    real omega, rho, rhoF(3,3), F(3,3), mom(3), rhoE, E0, E1, r
    integer ix,iy

    ! just 1-d at the moment (all y slices set equally)
    
    do ix=nb+1,nx-nb
       r = dx*(ix-nb-1) + 0.5*dx
       omega = 1.0/(this%stddev*sqrt(2*pi)) * exp(-r**2/(2*this%stddev**2))
       ! Miller-Collela elastic test:
       ! do iy=1,ny
       !    F = identity(3) / 1.1
       !    rho = F_density(eq%rho0,F)
       !    rhoF = rho*F

       !    mom = [0.0, 1.0, 0.0]

       !    E0 = eq%E(0.0,rhoF/rho)
       !    E1 = 10.0*E0

       !    rhoE = rho*((1.0-omega)*E0 + omega*E1)

       !    u(:,ix,iy) = [rho, mom, rhoF, rhoE]
       ! end do
       do iy=1,ny
          F = 0.0
          F(1,1) = 1.0 / 1.1
          F(2,2) = 1.0/((1.0 + 9.0*omega)*1.1)
          F(3,3) = (1.0 + 9.0*omega)/1.1

          rho = F_density(eq%rho0,F)
          rhoF = rho*F

          mom = [0.0, 1.0, 0.0]

          E0 = eq%E(0.0,rhoF/rho)
!         E1 = 10.0*E0

!         rhoE = rho*((1.0-omega)*E0 + omega*E1)
          rhoE = rho*E0

          u(:,ix,iy) = [rho, mom, rhoF, rhoE]
       end do
    end do
  end subroutine apply

  subroutine init_from_config(this,u)
    class(ic_gaussian) this
    integer, intent(in) :: u ! unit number
    real stddev
    namelist/ic_gaussian/stddev

    stddev = 0

    read (u,nml=ic_gaussian)
    this%stddev = stddev
  end subroutine init_from_config

end module m_ic_gaussian
