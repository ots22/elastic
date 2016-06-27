module m_state
  integer, parameter :: nq=14

  integer, parameter :: prim_rho=1, prim_v(3)=[2,3,4], prim_F(9)=[5,6,7,8,9,10,11,12,13], prim_S=14
  integer, parameter :: cons_rho=1, cons_mom(3)=[2,3,4], cons_rhoF(9)=[5,6,7,8,9,10,11,12,13], cons_rhoE=14
contains
  pure function Finger_G(F) result (G)
    use m_matutil, only: inv3
    real, intent(in) :: F(3,3)
    real G(3,3), F_inv(3,3)
    F_inv = inv3(F)
    G = matmul(F_inv, transpose(F_inv))
  end function Finger_G

  pure function Cauchy_c(F) result (c)
    use m_matutil, only: inv3
    real, intent(in) :: F(3,3)
    real c(3,3), F_inv(3,3)
    F_inv = inv3(F)
    c = matmul(transpose(F_inv), F_inv)
  end function Cauchy_c

  pure function dC_dF(q,p,F_inv)
    integer, intent(in) :: q, p
    real, intent(in) :: F_inv(3,3)
    real dC_dF(3,3), C(3,3)
    integer i,j
    C = matmul(transpose(F_inv), F_inv)
    forall(i=1:3,j=1:3) dC_dF(i,j) = -C(i,q)*F_inv(j,p) - C(j,q)*F_inv(i,p)
  end function dC_dF
  
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
