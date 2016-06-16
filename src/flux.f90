module m_flux
  use m_eos
  use m_state
contains
  ! cs is the cons state, ps is the (equivalent) prim state both are
  ! passed because normally they are both known anyway before the flux
  ! is computed.
  function flux(eq, cs, ps, dirn)
    use m_matutil, only: operator(.outer.)
    class(eos) eq
    real cs(nq), ps(nq), flux(nq)
    integer dirn
    intent(in) eq, cs, ps, dirn
    real v(3), F(3,3), S
    real mom(3), rhoF(3,3), rhoE
    real sigma(3,3)
    real mom_fl(3), rhoF_fl(3,3), rhoE_fl

    v = prim_get_v(ps)
    F = prim_get_F(ps)
    S = prim_get_S(ps)

    mom = cons_get_mom(cs)
    rhoF = cons_get_rhoF(cs)
    rhoE = cons_get_rhoE(cs)

    sigma = eq%stress(S, F)

    mom_fl = v(dirn) * mom - sigma(:,dirn)
    rhoF_fl = v(dirn) * rhoF - (v.outer.rhoF(dirn,:))
    rhoE_fl = v(dirn) * rhoE - dot_product(v, sigma(:,dirn))

    flux = [mom_fl, rhoF_fl, rhoE_fl]
  end function flux

  function lax_friedrichs_flux(L, R, fl_L, fl_R, dx_dt) result(fl)
    real, dimension(nq) :: L, R, fl_L, fl_R, fl
    real dx_dt
    fl = 0.5 * (fl_L + fl_R + dx_dt * (L - R))
  end function lax_friedrichs_flux

  function richtmeyer_flux(eq, L, R, fl_L, fl_R, dx_dt, dirn) result(fl)
    class(eos) :: eq
    real, dimension(nq) :: L, R, fl_L, fl_R, fl, ri_state, ri_pstate
    real dx_dt
    integer dirn
    ri_state = 0.5 * (L + R + (1.0/dx_dt) * (fl_L - fl_R))
    ri_pstate = cons_to_prim(eq, ri_state)
    fl = flux(eq, ri_state, ri_pstate, dirn)
  end function richtmeyer_flux

  function force_flux(eq, L, R, fl_L, fl_R, dx_dt, dirn) result(fl)
    class(eos) :: eq
    real, dimension(nq) :: L, R, fl_L, fl_R, fl
    real dx_dt
    integer dirn
    fl = 0.5 * (lax_friedrichs_flux(L,R,fl_L,fl_R,dx_dt) &
         &    + richtmeyer_flux (eq,L,R,fl_L,fl_R,dx_dt,dirn))
  end function force_flux
end module m_flux
