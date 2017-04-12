module m_state
  integer, parameter :: nq=15

  integer, parameter :: prim_rho=1, prim_v(3)=[2,3,4], prim_F(9)=[5,6,7,8,9,10,11,12,13], prim_S=14, prim_kappa=15
  integer, parameter :: cons_rho=1, cons_mom(3)=[2,3,4], cons_rhoF(9)=[5,6,7,8,9,10,11,12,13], cons_rhoE=14, cons_rhokappa=15
contains
  pure function Cauchy_Green_left(F) result (B)
    use m_matutil, only: inv3
    real, intent(in) :: F(3,3)
    real B(3,3)
    B = matmul(F, transpose(F))
  end function Cauchy_Green_left

  pure function Cauchy_Green_right(F) result (C)
    use m_matutil, only: inv3
    real, intent(in) :: F(3,3)
    real C(3,3)
    C = matmul(transpose(F), F)
  end function Cauchy_Green_right

  pure function dCinv_dF(p,q,Finv)
    integer, intent(in) :: q, p
    real, intent(in) :: Finv(3,3)
    real dCinv_dF(3,3), Cinv(3,3)
    integer i,j
    Cinv = matmul(Finv, transpose(Finv))
    forall(i=1:3,j=1:3) dCinv_dF(i,j) = -Cinv(i,q)*Finv(j,p) - Cinv(j,q)*Finv(i,p)
  end function dCinv_dF
  
  pure function prim_get_rho(u) result(rho)
    real, intent(in) :: u(nq)
    real rho
    rho = u(prim_rho)
  end function prim_get_rho

  pure function prim_get_v(u) result (v)
    real, intent(in) :: u(nq)
    real v(3)
    v = u(prim_v)
  end function prim_get_v

  pure function prim_get_F(u) result(F)
    real, intent(in) :: u(nq)
    real F(3,3)
    F = reshape(u(prim_F),[3,3])
  end function prim_get_F

  pure function prim_get_S(u) result(S)
    real, intent(in) :: u(nq)
    real S
    S = u(prim_S)
  end function prim_get_S
  
  pure function prim_get_kappa(u) result(kappa)
    real, intent(in) :: u(nq)
    real kappa
    kappa = u(prim_kappa)
  end function prim_get_kappa

  pure function cons_get_rho(u) result(rho)
    real, intent(in) :: u(nq)
    real rho
    rho = u(cons_rho)
  end function cons_get_rho

  pure function cons_get_mom(u) result (mom)
    real, intent(in) :: u(nq)
    real mom(3)
    mom = u(cons_mom)    
  end function cons_get_mom

  pure function cons_get_rhoF(u) result(rhoF)
    real, intent(in) :: u(nq)
    real rhoF(3,3)
    rhoF = reshape(u(cons_rhoF),[3,3])
  end function cons_get_rhoF

  pure function cons_get_rhoE(u) result(rhoE)
    real, intent(in) :: u(nq)
    real rhoE
    rhoE = u(cons_rhoE)
  end function cons_get_rhoE

  pure function cons_get_rhokappa(u) result(rhokappa)
    real, intent(in) :: u(nq)
    real rhokappa
    rhokappa = u(cons_rhokappa)
  end function cons_get_rhokappa

  pure function rhoF_density(rho0, rhoF) result(density)
    use m_matutil, only: det3
    real, intent(in) :: rho0, rhoF(3,3)
    real density
    density = sqrt(det3(rhoF) / rho0)
  end function rhoF_density
  
  pure function F_density(rho0, F) result(density)
    use m_matutil
    real, intent(in) :: rho0, F(3,3)
    real density
    density = rho0 / det3(F)
  end function F_density
end module m_state
