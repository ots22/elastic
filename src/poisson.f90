module m_poisson
contains
  subroutine poisson_2d(uk0,uk,h)
    use m_error, only: assert
    use m_matutil, only: grad_2d, div_2d
    integer NX, NY, i, j, k
    integer, parameter :: NCHECK=10
    real, parameter :: TOL=1.0E-8
    integer, parameter :: MAXSTEPS=5000

    real, intent(in) :: h ! the grid resolution
    real, intent(in) :: uk0(:,:) ! source terms and boundary
    real, intent(out) :: uk(:,:) ! the solution (same shape as uk0)
    real ukp1(size(uk0,1),size(uk0,2))
    real diff, diffmax
    
    NX = size(uk0,1)
    NY = size(uk0,2)

    call assert(all(shape(uk0).eq.shape(uk)), &
         & "initial conditions and solution arrays must have the same shape")

    uk = 0
    uk(2:NX-1,2:NY-1) = uk0(2:NX-1,2:NY-1)
    diffmax = TOL + 1.0
    do k=1,MAXSTEPS
       do i = 2, NX-1
          do j = 2, NY-1
             ukp1(i,j) = 0.25*(H*H*uk0(i,j)  &
                    + uk(i-1,j) + uk(i,j-1)  &
                    + uk(i+1,j) + uk(i,j+1))
          enddo
       enddo
       if (mod(k,NCHECK) .eq. 0) then
          diffmax = 0.0
          do i = 2, NX-1
             do j = 2, NY-1
                diff = abs(ukp1(i,j) - uk(i,j))
                if (diff .gt. diffmax) diffmax = diff
             enddo
          enddo
       endif
       do i = 2, NX-1
          do j = 2, NY-1
             uk(i,j) = ukp1(i,j)
          enddo
       enddo
       if (diffmax .le. TOL) exit
    enddo
  end subroutine poisson_2d
end module m_poisson
