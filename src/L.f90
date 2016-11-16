module m_L
  use m_time
contains
  function advance_solution(eq, sol, solp, uLp, uRp, &
       &                    dx_dt, dirn) result(sol_next)
    use m_eos, only: eos
    use m_state, only: nq
    use m_domain, only: nquad
    use m_ssprk3, only: ssprk3
    use m_error, only: assert
    class(eos), intent(in)  :: eq
    real, intent(in)  :: sol(:,:,-2:)
    real sol_next(size(sol,1),size(sol,2))
    real, intent(in)  :: solp(:,:,-2:)
    real, intent(in)  :: uRp(:,:,-2:), uLp(:,:,-2:)
    real, intent(in) :: dx_dt
    integer, intent(in) :: dirn

    real :: work(size(sol,1),size(sol,2))
    real, dimension(size(solp,1),size(solp,2),-2:2) :: solp_work
    real, dimension(size(sol,1),size(sol,2),-2:2)   :: uLp_work, uRp_work
    real, dimension(size(sol,1),size(sol,2),nquad)  :: uLpq_work, uRpq_work
    real :: fl(size(sol,1),size(sol,2))
    real :: d(size(sol,1),size(sol,2))

    logical first_call_to_L

    integer nx
    nx = size(sol,2)

    ! call assert(size(d,1).eq.nq.and.size(d,2).ge.nx, '')
    ! call assert(all(lbound(sol).le.[1,1,-2].and.ubound(sol).ge.[nq,nx,2]),  'sol shape incorrect')
    ! call assert(all(lbound(solp).le.[1,1,-2].and.ubound(solp).ge.[nq,nx,2]),'solp shape incorrect')
    ! call assert(all(lbound(solp_work).le.[1,1,-2].and.ubound(solp_work).ge.[nq,nx,2]),'solp_work shape incorrect')
    ! call assert(all(lbound(work).le.[1,1].and.ubound(work).ge.[nq,nx]),'work shape incorrect')
    ! call assert(all(lbound(uLp).le.[1,1,-2,0].and.ubound(uLp).ge.[nq,nx,2,nquad]),'uLp shape incorrect')
    ! call assert(all(lbound(uRp).le.[1,1,-2,0].and.ubound(uRp).ge.[nq,nx,2,nquad]),'uRp shape incorrect')

    solp_work = solp
    uLp_work  = uLp
    uRp_work  = uRp

    first_call_to_L = .true.
    sol_next(:,:) = ssprk3(L, sol(:,:,0), work)

  contains
    subroutine L(u,q)
      use m_error, only: assert
!      use m_state, only: nq
      use m_weno, only: wenom_reconstruct, wenom_reconstruct_gqp, median_state
      use m_domain, only: quad_weight
      use m_eos, only: cons_to_prim, prim_to_cons
      use m_flux, only: flux, force_flux

      ! assumed shape for the ssprk3 interface
      real, intent(in)  :: u(:,:)
      real, intent(out) :: q(:,:)
      ! temporaries
      real, dimension(nq) :: Li, Ri, Lpi, Rpi, fl_Li, fl_Ri
      integer ix, iy, nx, iq, iquad

      real curvature(nq,-1:1)

      nx = size(u,2)

!      call assert(all(shape(u).ge.[nq,nx].and.shape(q).ge.[nq,nx]),   &
!           &     "Sizes of solution arrays passed to subroutine 'L' must have extents at least [nq,nx]")

!     some of the parameters don't change on repeated calls to L
!     within `advance' - don't recompute them
      if (.not.first_call_to_L) then
         do ix = 1, nx
            solp_work(:,ix,0) = cons_to_prim(eq, u(:,ix))
         end do
         do ix = 2, nx-1
            d(:,ix) = solp_work(:,ix+1,0) - 2*solp_work(:,ix,0) + solp_work(:,ix-1,0)
         end do
         do ix = 1+2, nx-2
            do iq = 1,nq
               call wenom_reconstruct(solp_work(iq,ix-2:ix+2,0),uLp_work(iq,ix,0),uRp_work(iq,ix,0))
            end do
            uLp_work(:,ix,0) = median_state(d(:,ix+1:ix-1:-1), &
                 solp_work(:,ix+1:ix-1:-1,0), uLp_work(:,ix,0))
            uRp_work(:,ix,0) = median_state(d(:,ix-1:ix+1), &
                 solp_work(:,ix-1:ix+1,0), uRp_work(:,ix,0))
         end do
      else
         first_call_to_L = .false.
      end if

!     reconstruct at the quadrature points
      do ix = 1+2, nx-2
         do iq = 1,nq
            call wenom_reconstruct_gqp(uLp_work(iq,ix,:), &
                 uLpq_work(iq,ix,1),uLpq_work(iq,ix,2))
            call wenom_reconstruct_gqp(uRp_work(iq,ix,:), &
                 uRpq_work(iq,ix,1),uRpq_work(iq,ix,2))
         end do
         ! do iy = -1,1
         !    curvature(:,iy) = uLp_work(:,ix,iy+1) - 2*uLp_work(:,ix,iy) + uLp_work(:,ix,iy-1)
         ! end do
         ! do iquad = 1,nquad
         !    uLpq_work(:,ix,iquad) = median_state(curvature,uLp_work(:,ix,-1:1),uLpq_work(:,ix,iquad))
         ! end do
         ! do iy = -1,1
         !    curvature(:,iy) = uRp_work(:,ix,iy+1) - 2*uRp_work(:,ix,iy) + uRp_work(:,ix,iy-1)
         ! end do
         ! do iquad = 1,nquad
         !    uRpq_work(:,ix,iquad) = median_state(curvature,uRp_work(:,ix,-1:1),uRpq_work(:,ix,iquad))
         ! end do
      end do

!     flux through the face is an integral average = a sum over the
!     quadrature points
      fl = 0
      do iquad = 1,nquad
         do ix = 1+2, nx-3
            ! at ix, the flux is computed for the interface between ix and
            ! ix+1 so the state to the left of the *interface* is the state
            ! at the right of the left *cell* and vice versa:

            Lpi = uRpq_work(:,ix,  iquad);     Li = prim_to_cons(eq,Lpi)
            Rpi = uLpq_work(:,ix+1,iquad);     Ri = prim_to_cons(eq,Rpi)

            fl_Li = flux(eq,Li,Lpi,dirn)

            fl_Ri = flux(eq,Ri,Rpi,dirn)
            fl(:,ix) = fl(:,ix) + quad_weight(iquad) * force_flux(eq, Li, Ri, fl_Li, fl_Ri, dx_dt, dirn)
         end do
      end do

!     the update to the solution to use
      q = 0
      do ix = 1+3, nx-3
         q(:,ix) = (fl(:,ix-1) - fl(:,ix)) / dx_dt
      end do
    end subroutine L
  end function advance_solution
end module m_L
