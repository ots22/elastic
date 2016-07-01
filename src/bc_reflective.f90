module m_bc_reflective
  use m_bc

  private

  type, public, extends(bc) :: bc_reflective
   contains
     procedure apply
  end type bc_reflective

contains
  subroutine apply(this,eq,sol,nb)
    use m_state
    use m_error
    use m_eos
    use m_matutil, only: identity

    class(bc_reflective) this
    class(eos), intent(in) :: eq
    real, intent(inout) :: sol(:,:,:)
    integer, intent(in) :: nb

    real Rx(3,3), Ry(3,3)
    integer ib, ib_refl, ib_right, ib_right_refl, ix, iy, nx, ny

    call assert(size(sol,1).eq.nq, "mismatch in size of solution array &
         & first dimension and number of state components")
    nx = size(sol,2)
    ny = size(sol,3)

    ! Reflection matrices about x and y
    Rx = identity(3); Rx(1,1) = -1
    Ry = identity(3); Ry(2,2) = -1

!                                |
!                                |
!  +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----
!  |     |     |     |     |     |     |     |2*nb |     |2*nb |
!  |   1 | ... | ib  | ... | nb  |nb+1 | ... |-ib+1| ... | +1  | ...
!  +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----
!                                |\_____________ _____________/
!                                |              v
!                                |              nb
    do ib=1,nb
       do iy=1,ny
!         ib and ib_right indexes a boundary cell (left or right), with
!         ib_refl and ib_right_refl the cells these copy data from
          ib_refl = 2*nb-ib+1
          ib_right = nx+1-ib
          ib_right_refl = nx+1-ib_refl
          
!         scalar components of the state are easy:
          sol(cons_rho,ib,iy)        = sol(cons_rho,ib_refl,iy)
          sol(cons_rho,ib_right,iy)  = sol(cons_rho,ib_right_refl,iy)
          
          sol(cons_rhoE,ib,iy)       = sol(cons_rhoE,ib_refl,iy)
          sol(cons_rhoE,ib_right,iy) = sol(cons_rhoE,ib_right_refl,iy)

!         the only vector component transforms as v'_{i} = Rx_{ij}*v_{j}
          sol(cons_mom,ib,iy)        = matmul(Rx,sol(cons_mom,ib_refl,iy))
          sol(cons_mom,ib_right,iy)  = matmul(Rx,sol(cons_mom,ib_right_refl,iy))
          
!         the tensor component transforms as F'_{ij} = Rx_{ik}*F_{kl}*Rx_{lj} (since Rx is self-inverse)
          sol(cons_rhoF,ib,iy)       = [matmul(Rx,matmul(cons_get_rhoF(sol(:,ib_refl,iy)),Rx))]
          sol(cons_rhoF,ib_right,iy) = [matmul(Rx,matmul(cons_get_rhoF(sol(:,ib_right_refl,iy)),Rx))]
       end do
    end do

! !   same again for the y direction
!     do ib=1,nb
!        do ix=1,nx
! !         ib and ib_right indexes a boundary cell (left or right), with
! !         ib_refl and ib_right_refl the cells these copy data from
!           ib_refl = 2*nb-ib+1
!           ib_right = ny+1-ib
!           ib_right_refl = ny+1-ib_refl
          
! !         scalar components of the state are easy:
!           sol(cons_rho,ix,ib)        = sol(cons_rho,ix,ib_refl)
!           sol(cons_rho,ix,ib_right)  = sol(cons_rho,ix,ib_right_refl)
          
!           sol(cons_rhoE,ix,ib)       = sol(cons_rhoE,ix,ib_refl)
!           sol(cons_rhoE,ix,ib_right) = sol(cons_rhoE,ix,ib_right_refl)

! !         the only vector component transforms as v'_{i} = Ry_{ij}*v_{j}
!           sol(cons_mom,ix,ib)        = matmul(Ry,sol(cons_mom,ix,ib_refl))
!           sol(cons_mom,ix,ib_right)  = matmul(Ry,sol(cons_mom,ix,ib_right_refl))
          
! !         the tensor component transforms as F'_{ij} = Ry_{ik}*F_{kl}*Ry_{lj} (since Ry is self-inverse)
!           sol(cons_rhoF,ix,ib)       = [matmul(Ry,matmul(cons_get_rhoF(sol(:,ix,ib_refl)),Ry))]
!           sol(cons_rhoF,ix,ib_right) = [matmul(Ry,matmul(cons_get_rhoF(sol(:,ix,ib_right_refl)),Ry))]
!       end do
!    end do
 
    forall (ib=1:nb)
       sol(:,:,ib) = sol(:,:,nb+1)
       sol(:,:,ny-ib+1) = sol(:,:,ny-nb)
    end forall

  end subroutine apply
end module m_bc_reflective
