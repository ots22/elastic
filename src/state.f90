module m_state
  integer, parameter :: nq=13

  integer, parameter :: prim_v(3)=[1,2,3], prim_F(9)=[4,5,6,7,8,9,10,11,12], prim_S=13
  integer, parameter :: cons_mom(3)=[1,2,3], cons_rhoF(9)=[4,5,6,7,8,9,10,11,12], cons_rhoE=13
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
  
  pure function prim_get_v(u) result (v)
    real, intent(in) :: u(nq)
    real v(3)
    v = u(1:3)
  end function prim_get_v

  pure function prim_get_F(u) result(F)
    real, intent(in) :: u(nq)
    real F(3,3)
    F = reshape(u(4:12),[3,3])
  end function prim_get_F

  pure function prim_get_S(u) result(S)
    real, intent(in) :: u(nq)
    real S
    S = u(13)
  end function prim_get_S
  
  pure function cons_get_mom(u) result (mom)
    real, intent(in) :: u(nq)
    real mom(3)
    mom = u(1:3)    
  end function cons_get_mom

  pure function cons_get_rhoF(u) result(rhoF)
    real, intent(in) :: u(nq)
    real rhoF(3,3)
    rhoF = reshape(u(4:12),[3,3])
  end function cons_get_rhoF

  pure function cons_get_rhoE(u) result(rhoE)
    real, intent(in) :: u(nq)
    real rhoE
    rhoE = u(13)
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
