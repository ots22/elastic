module m_ic
  use m_configurable
  
  private
  public ic

  type, abstract, extends(configurable) :: ic
   contains
     procedure(apply), deferred :: apply
  end type ic

  abstract interface
     subroutine apply(this,eq,u)
       use m_eos
       import ic
       class(ic) this
       class(eos), intent(in) :: eq
       real, intent(out) :: u(:,:,:)
     end subroutine apply
  end interface
end module m_ic

  ! subroutine miller_axisymmetric(eq,u,a,dx)
  !   use m_state, only: nq
  !   use m_matutil, only: identity, pi
  !   use m_eos, only: eos
  !   use m_error, only: assert
  !   class(eos), intent(in) :: eq
  !   real, intent(out) :: u(:,:,:)
  !   real, intent(in) :: a, dx
  !   real mom(3), rhoF(3,3), rhoE, r2
  !   integer nx, ny, ix, iy

  !   call assert(size(u,1).eq.nq, 'size of first dimension of u in &
  !        &miller_axisymmetric must be nq')
  !   nx = size(u,2)
  !   ny = size(u,3)    

  !   mom = 0
  !   rhoF = eq%rho0 * identity(3)

  !   do iy=1,ny
  !      do ix=1,nx
  !         r2 = ((dx*(ix-nx/2))**2 + (dx*(iy-ny/2))**2)
  !         rhoE = (0.05/(a*sqrt(pi))) * exp(-r2/(2*a**2))
  !         u(:,ix,iy) = [mom, rhoF, rhoE]
  !      end do
  !   end do
  ! end subroutine miller_axisymmetric
