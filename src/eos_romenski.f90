module m_eos_romenski
  use m_eos

  public :: eos_romenski
  
  private

  type, extends(eos) :: eos_romenski
     real K0, B0, alpha, beta, gamma, cv, T0
   contains
     procedure E
     procedure S
     procedure stress
     procedure init_from_config
  end type eos_romenski

contains
  function E(this, S, F, kappa) result(e_internal)
    use m_state, only: Cauchy_Green_left
    use m_matutil, only: inv3
    class(eos_romenski) :: this
    real, intent(in) :: S, F(3,3), kappa
    real e_internal
    real Binv(3,3)
    Binv = inv3(Cauchy_Green_left(F))
    associate (&
         Binv11 => Binv(1,1),     &
         Binv12 => Binv(1,2),     &
         Binv13 => Binv(1,3),     &
         Binv22 => Binv(2,2),     &
         Binv23 => Binv(2,3),     &
         Binv33 => Binv(3,3),     &
         rho0   => this%rho0,     &
         K0     => this%K0,       &
         B0     => this%B0,       &
         alpha  => this%alpha,    &
         beta   => this%beta,     &
         gamma  => this%gamma,    &
         cv     => this%cv,       &
         T0     => this%T0)
      include 'eos/romenski_energy.inc'
    end associate
  end function E

  function S(this, E, F, kappa) result(entropy)
    use m_state, only: Cauchy_Green_left
    use m_matutil, only: inv3
    class(eos_romenski) :: this
    real, intent(in) :: E, F(3,3), kappa
    real entropy
    real Binv(3,3)
    Binv = inv3(Cauchy_Green_left(F))
    associate (&
         Binv11 => Binv(1,1),     &
         Binv12 => Binv(1,2),     &
         Binv13 => Binv(1,3),     &
         Binv22 => Binv(2,2),     &
         Binv23 => Binv(2,3),     &
         Binv33 => Binv(3,3),     &
         rho0   => this%rho0,     &
         K0     => this%K0,       &
         B0     => this%B0,       &
         alpha  => this%alpha,    &
         beta   => this%beta,     &
         gamma  => this%gamma,    &
         cv     => this%cv,       &
         T0     => this%T0)
      include 'eos/romenski_entropy.inc'
    end associate
  end function S

  function stress(this, S, F, kappa)
    use m_state, only: Cauchy_Green_left
    use m_matutil, only: inv3
    class(eos_romenski) :: this
    real, intent(in) :: S, F(3,3), kappa
    real stress(3,3)
    real Binv(3,3)
    Binv = inv3(Cauchy_Green_left(F))
    associate (&
         Binv11 => Binv(1,1),     &
         Binv12 => Binv(1,2),     &
         Binv13 => Binv(1,3),     &
         Binv22 => Binv(2,2),     &
         Binv23 => Binv(2,3),     &
         Binv33 => Binv(3,3),     &
         rho0   => this%rho0,     &
         K0     => this%K0,       &
         B0     => this%B0,       &
         alpha  => this%alpha,    &
         beta   => this%beta,     &
         gamma  => this%gamma,    &
         cv     => this%cv,       &
         T0     => this%T0)
      include 'eos/romenski_stress.inc'
    end associate
  end function stress

  subroutine init_from_config(this,u)
    use m_error
    class(eos_romenski) this
    integer, intent(in) :: u
    real rho0,c0,l0,K0,B0,alpha,beta,gamma,cv,T0
    real, parameter :: sentinal=huge(1.0)
    namelist/EOS_ROMENSKI/rho0,c0,l0,K0,B0,alpha,beta,gamma,cv,T0
    rho0 = 0
    c0 = sentinal
    l0 = sentinal
    K0 = sentinal
    B0 = sentinal
    alpha = 0
    beta = 0
    gamma = 0
    cv = 0
    T0 = sentinal
    
    read (u,nml=EOS_ROMENSKI)

    call assert(rho0.gt.0, "eos_romenski: density must be positive")
    if (l0.ne.sentinal) B0 = l0**2
    if (l0.ne.sentinal.and.c0.ne.sentinal) K0 = c0**2 - (4.0/3.0)*l0**2
    call assert(B0.ge.0.and.B0.ne.sentinal, "eos_romenski%init_from_config: B0 must be set positive (can be set via l0)")
    call assert(K0.ge.0.and.K0.ne.sentinal, "eos_romenski%init_from_config: K0 must be set positive (can be set via l0 and c0)")
    call assert(cv.ge.0, "eos_romenski%init_from_config: heat capacity cv must be set positive")
    call assert(T0.ne.sentinal, "eos_romenski%init_from_config: T0 must be set")
    if (any([mod(alpha,1.0),mod(beta,1.0),mod(gamma,1.0)].eq.0)) then
       call warn("eos_romenski%init_from_config: exact integers for alpha, beta or gamma will& 
            & trigger a performance bug when using GLIBC")
    endif

    this%rho0 = rho0
    this%K0 = K0
    this%B0 = B0
    this%alpha = alpha
    this%beta = beta
    this%gamma = gamma
    this%cv = cv
    this%T0 = T0

  end subroutine init_from_config
end module m_eos_romenski
