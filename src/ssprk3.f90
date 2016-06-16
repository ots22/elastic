module m_ssprk3
  use m_error, only: assert
contains
! 'ssprk3' takes one integration step, on an nq component solution u
! of size nx, using SSPRK(3,3).  'L' is a function that takes the
! solution 'u', and computes the update for a forward Euler step into
! du (for a flux based method, this would typically be 
! dt/dx * (flux_{i} - flux_{i-1})).  The routine L should handle the
! boundary conditions.

  function ssprk3(L,u,work) result(u_next)
    interface
       subroutine L(u,du)
         real, intent(in) :: u(:,:)
         real, intent(out) :: du(:,:)
       end subroutine L
    end interface
    real, intent(in) :: u(:,:)
    real, intent(inout) :: work(:,:)
    real u_next(size(u,1),size(u,2))

!    call assert(all(shape(work).eq.shape(u)), &
!         & 'shape of work array to ssprk3 must match u')
    
    associate(du => work)
      call L(u,du)
      u_next = u + du
      call L(u_next,du)
      u_next = (3*u + u_next + du)/4
      call L(u_next,du)
      u_next = (u + 2*u_next + 2*du)/3
    end associate
  end function ssprk3
end module m_ssprk3
