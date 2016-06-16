module m_bc_offset_periodic
  use m_bc

  private
  public bc_offset_periodic 

  type, extends(bc) :: bc_offset_periodic
   contains
     procedure apply
  end type bc_offset_periodic

contains
  subroutine apply(this,eq,sol,nb)
    use m_state, only: nq
    use m_error
    use m_eos

    class(bc_offset_periodic) this
    class(eos), intent(in) :: eq
    real, intent(inout) :: sol(:,:,:)
    integer, intent(in) :: nb ! number of boundary cells

    integer ib, ix
    integer nx, ny

    if (.false.) call no_op(this,eq)

    call assert(size(sol,1).eq.nq, "mismatch in size of solution array &
         & first dimension and number of state components")
    nx = size(sol,2)
    ny = size(sol,3)

    forall (ib=1:nb)
       sol(:,:,ib) = sol(:,:,nb+1)
       sol(:,:,ny-ib+1) = sol(:,:,ny-nb)
    end forall
    forall (ib=1:nb)
       sol(:,ib,:) = sol(:,nb+1,:)
       sol(:,nx-ib+1,:) = sol(:,nx-nb,:)
    end forall
    do ix=1,nx
       do ib=1,nb
          if (ix+ny-2*nb.le.nx)  sol(:, ix, ib) = sol(:, ix+ny-2*nb, ny+ib-2*nb)
          if (ix-ny+2*nb.ge.1)   sol(:, ix, ny-nb+ib) = sol(:, ix-ny+2*nb, nb+ib)
       end do
    end do
  end subroutine apply
end module m_bc_offset_periodic
