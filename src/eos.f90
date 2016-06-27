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
  end type eos
  
  abstract interface
     function E(this, S, F)
       import eos
       intent(in) S, F
       class(eos) this
       real E, S, F(3,3)
     end function E

     function S(this, E, F)
       import eos
       intent(in) E, F
       class(eos) this
       real S, E, F(3,3)
     end function S

     function stress(this, S, F)
       import eos
       intent(in) S, F
       class(eos) this
       real stress(3,3), S, F(3,3)
     end function stress      
  end interface

contains

! convert a state vector c of conserved variables to a vector p of
! primitive variables
  function cons_to_prim(eq, c) result(p)
    use m_state
    intent(in) eq, c
    class(eos) eq
    real, dimension(nq) :: c, p
    real spfc_vol
    real rho, v(3), F(3,3), E, S

    spfc_vol = 1/rhoF_density(eq%rho0, cons_get_rhoF(c))

    rho = cons_get_rho(c)
    v = spfc_vol * cons_get_mom(c)
    F = spfc_vol * cons_get_rhoF(c)
    E = spfc_vol * cons_get_rhoE(c) - 0.5*dot_product(v,v)
    S = eq%S(E, F)

    p = [rho, v, F, S]
  end function cons_to_prim

  function prim_to_cons(eq, p) result(c)
    use m_state
    intent(in) eq, p
    class(eos) eq
    real, dimension(nq) :: p, c
    real density 
    real v(3), F(3,3), S 
    real rho, mom(3), rhoF(3,3), rhoE
    
    rho = prim_get_rho(p)
    v = prim_get_v(p)
    F = prim_get_F(p)
    S = prim_get_S(p)

    density = F_density(eq%rho0, F)

    mom  = density * v
    rhoF = density * F
    rhoE = density * (eq%E(S, F) + 0.5*dot_product(v,v))

    c = [rho, mom, rhoF, rhoE]
  end function prim_to_cons
end module m_eos
