module m_eos
  use m_configurable

  private
  public :: eos, cons_to_prim, prim_to_cons

  type, abstract, extends(configurable) :: eos
     real :: rho0
   contains
     procedure(E), deferred :: E
     procedure(S), deferred :: S
     procedure(stress), deferred :: stress
     procedure stress_E
     procedure hardening
     procedure dhardening_dkappa
     procedure dstress_dg_E
     procedure dstress_dkappa_E
  end type eos
  
  abstract interface
     function E(this, S, F, kappa)
       import eos
       intent(in) S, F, kappa
       class(eos) this
       real E, S, F(3,3), kappa
     end function E

     function S(this, E, F, kappa)
       import eos
       intent(in) E, F, kappa
       class(eos) this
       real S, E, F(3,3), kappa
     end function S

     function stress(this, S, F, kappa)
       import eos
       intent(in) S, F, kappa
       class(eos) this
       real stress(3,3), S, F(3,3), kappa
     end function stress
  end interface

contains
  function stress_E(this, E, F, kappa)
    intent(in) E, F, kappa
    class(eos) this
    real stress_E(3,3), S, F(3,3), kappa, E
    S = this%S(E,F,kappa)
    stress_E = this%stress(S,F,kappa)
  end function stress_E
  
  function hardening(this, S, F, kappa)
    use m_error
    class(eos) this
    real, intent(in) :: S, F(3,3), kappa
    real hardening
    hardening = 0.0
    call warn("hardening unimplemented, assuming 0.0")
  end function hardening

  function dhardening_dkappa(this, S, F, kappa)
    use m_error
    class(eos) this
    real, intent(in) :: S, F(3,3), kappa
    real dhardening_dkappa
    dhardening_dkappa = 0.0
    call warn("dhardening_dkappa unimplemented, assuming 0.0")
  end function dhardening_dkappa

  ! finite differences for now: implement something in the python script later
  function dstress_dg_E(this, S, F, kappa) result(result)
    use m_matutil, only: inv3
    class(eos) :: this
    real, intent(in) :: S, F(3,3), kappa
    real result(3,3,3,3), g(3,3), gplus(3,3), gminus(3,3), Fplus(3,3), Fminus(3,3), E
    real, parameter :: h=1.0E-7
    integer k,l ! loop counters
    g = inv3(F)
    E = this%E(S,F,kappa)
    do k=1,3; do l=1,3
       gplus  = g
       gminus = g
       gplus(k,l)  = g(k,l) + h
       gminus(k,l) = g(k,l) - h
       Fplus  = inv3(gplus)
       Fminus = inv3(gminus)
       result(:,:,k,l) = (this%stress_E(E,Fplus,kappa) - this%stress_E(E,Fminus,kappa))/(2*h)
    end do; end do
  end function dstress_dg_E

  function dstress_dkappa_E(this, S, F, kappa)
    use m_error
    class(eos) :: this
    real, intent(in) :: S, F(3,3), kappa
    real dstress_dkappa_E(3,3)
    real E
    real, parameter :: h=1.0E-7
    E = this%E(S,F,kappa)
    dstress_dkappa_E = (this%stress_E(E,F,kappa+h) - this%stress_E(E,F,kappa-h))/(2*h)
  end function dstress_dkappa_E

! convert a state vector c of conserved variables to a vector p of
! primitive variables
  function cons_to_prim(eq, c) result(p)
    use m_state
    intent(in) eq, c
    class(eos) eq
    real, dimension(nq) :: c, p
    real spfc_vol
    real rho, v(3), F(3,3), E, S, kappa

    spfc_vol = 1/rhoF_density(eq%rho0, cons_get_rhoF(c))

    rho = cons_get_rho(c)
    v = spfc_vol * cons_get_mom(c)
    F = spfc_vol * cons_get_rhoF(c)
    E = spfc_vol * cons_get_rhoE(c) - 0.5*dot_product(v,v)
    kappa = spfc_vol * cons_get_rhokappa(c)
    S = eq%S(E, F, kappa)

    p = [rho, v, F, S, kappa]
  end function cons_to_prim

  function prim_to_cons(eq, p) result(c)
    use m_state
    intent(in) eq, p
    class(eos) eq
    real, dimension(nq) :: p, c
    real density 
    real v(3), F(3,3), S, kappa
    real rho, mom(3), rhoF(3,3), rhoE, rhokappa
    
    rho = prim_get_rho(p)
    v = prim_get_v(p)
    F = prim_get_F(p)
    S = prim_get_S(p)
    kappa = prim_get_kappa(p)

    density = F_density(eq%rho0, F)

    mom  = density * v
    rhoF = density * F
    rhoE = density * (eq%E(S, F, kappa) + 0.5*dot_product(v,v))
    rhokappa = density * kappa

    c = [rho, mom, rhoF, rhoE, rhokappa]
  end function prim_to_cons
end module m_eos
