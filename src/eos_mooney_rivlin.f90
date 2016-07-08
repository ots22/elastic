module m_eos_mooney_rivlin
  use m_eos

  private

  type, public, extends(eos) :: eos_mooney_rivlin
     real lambda_0, lambda_s, mu_0, mu_s, theta_0, theta_1
   contains
     procedure E
     procedure S
     procedure stress
     procedure hardening
     procedure dhardening_dkappa
!!!  commented out to use the finite-difference formula in the base class
     procedure dstress_dkappa_E
     procedure init_from_config
  end type eos_mooney_rivlin

contains
    function E(this, S, F, kappa) result(e_internal)
    use m_state, only: Cauchy_Green_right
    class(eos_mooney_rivlin) :: this
    real, intent(in) :: S, F(3,3), kappa
    real e_internal
    real C(3,3)
    !! change this to work-hardening parameter when eos supports it
    C = Cauchy_Green_right(F)
    associate (&
         C11 => C(1,1),  &
         C12 => C(1,2),  &
         C13 => C(1,3),  &
         C22 => C(2,2),  &
         C23 => C(2,3),  &
         C33 => C(3,3),  &
         rho0   => this%rho0,     &
         lambda_0 => this%lambda_0, &
         lambda_s => this%lambda_s, &
         mu_0 => this%mu_0, &
         mu_s => this%mu_s, &
         theta_0 => this%theta_0, &
         theta_1 => this%theta_1)
      include 'eos/mooney_rivlin_energy.inc'
    end associate
  end function E

  function S(this, E, F, kappa) result(entropy)
    use m_state, only: Cauchy_Green_right
    class(eos_mooney_rivlin) :: this
    real, intent(in) :: E, F(3,3), kappa
    real entropy
    real C(3,3)
    C = Cauchy_Green_right(F)
    associate (&
         C11 => C(1,1),  &
         C12 => C(1,2),  &
         C13 => C(1,3),  &
         C22 => C(2,2),  &
         C23 => C(2,3),  &
         C33 => C(3,3),  &
         rho0   => this%rho0,     &
         lambda_0 => this%lambda_0, &
         lambda_s => this%lambda_s, &
         mu_0 => this%mu_0, &
         mu_s => this%mu_s, &
         theta_0 => this%theta_0, &
         theta_1 => this%theta_1)
      include 'eos/mooney_rivlin_entropy.inc'
    end associate
  end function S

  function stress(this, S, F, kappa)
    use m_state, only: Cauchy_Green_right
    class(eos_mooney_rivlin) :: this
    real, intent(in) :: S, F(3,3), kappa
    real stress(3,3)
    real C(3,3)
    C = Cauchy_Green_right(F)
    associate (&
         C11 => C(1,1),  &
         C12 => C(1,2),  &
         C13 => C(1,3),  &
         C22 => C(2,2),  &
         C23 => C(2,3),  &
         C33 => C(3,3),  &
         F11 => F(1,1),  &
         F12 => F(1,2),  &
         F13 => F(1,3),  &
         F21 => F(2,1),  &
         F22 => F(2,2),  &
         F23 => F(2,3),  &
         F31 => F(3,1),  &
         F32 => F(3,2),  &
         F33 => F(3,3),  &
         rho0   => this%rho0,     &
         lambda_0 => this%lambda_0, &
         lambda_s => this%lambda_s, &
         mu_0 => this%mu_0, &
         mu_s => this%mu_s, &
         theta_0 => this%theta_0, &
         theta_1 => this%theta_1)
      include 'eos/mooney_rivlin_stress.inc'
    end associate
  end function stress

! \rho * \partial E/\partial kappa
  function hardening(this, S, F, kappa)
    class(eos_mooney_rivlin) this
    real, intent(in) :: S, F(3,3), kappa
    real hardening
    hardening = this%theta_0 * (1 - exp(-this%theta_1*kappa))
  end function hardening

  function dhardening_dkappa(this, S, F, kappa)
    class(eos_mooney_rivlin) this
    real, intent(in) :: S, F(3,3), kappa
    real dhardening_dkappa
    dhardening_dkappa = this%theta_0 * this%theta_1 * exp(-this%theta_1*kappa)
  end function dhardening_dkappa

  function dstress_dkappa_E(this, S, F, kappa)
    use m_state, only: Cauchy_Green_right
    class(eos_mooney_rivlin) :: this
    real, intent(in) :: S, F(3,3), kappa
    real dstress_dkappa_E(3,3), E
    real C(3,3)
    C = Cauchy_Green_right(F)
    E = this%E(S,F,kappa)
    associate (&
         C11 => C(1,1),  &
         C12 => C(1,2),  &
         C13 => C(1,3),  &
         C22 => C(2,2),  &
         C23 => C(2,3),  &
         C33 => C(3,3),  &
         F11 => F(1,1),  &
         F12 => F(1,2),  &
         F13 => F(1,3),  &
         F21 => F(2,1),  &
         F22 => F(2,2),  &
         F23 => F(2,3),  &
         F31 => F(3,1),  &
         F32 => F(3,2),  &
         F33 => F(3,3),  &
         rho0   => this%rho0,     &
         lambda_0 => this%lambda_0, &
         lambda_s => this%lambda_s, &
         mu_0 => this%mu_0, &
         mu_s => this%mu_s, &
         theta_0 => this%theta_0, &
         theta_1 => this%theta_1)
      include 'eos/mooney_rivlin_dstress_dkappa_e.inc'
    end associate
  end function dstress_dkappa_E

  subroutine init_from_config(this,u)
    use m_error
    class(eos_mooney_rivlin) this
    integer, intent(in) :: u
    real rho0,lambda_0,lambda_s,mu_0,mu_s,theta_0,theta_1
    real, parameter :: sentinal=huge(1.0)
    namelist/EOS_MOONEY_RIVLIN/rho0,lambda_0,lambda_s,mu_0,mu_s,theta_0,theta_1
    rho0 = 0
    lambda_0 = sentinal
    lambda_s = sentinal
    mu_0 = sentinal
    mu_s = sentinal 
    theta_0 = sentinal
    theta_1 = sentinal
    
    read (u,nml=EOS_MOONEY_RIVLIN)

    call assert(rho0.gt.0, "eos_mooney_rivlin: density must be positive")

    this%rho0 = rho0
    this%lambda_0 =lambda_0
    this%lambda_s = lambda_s 
    this%mu_0 = mu_0 
    this%mu_s = mu_s 
    this%theta_0 = theta_0 
    this%theta_1 = theta_1
  end subroutine init_from_config
end module m_eos_mooney_rivlin
