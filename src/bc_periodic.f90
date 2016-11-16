module m_bc_periodic
  use m_bc

  private
  public bc_periodic 

  type, extends(bc) :: bc_periodic
   contains
     procedure apply
  end type bc_periodic

contains
  subroutine apply(this,eq,sol,nb)
    use m_state, only: nq
    use m_error
    use m_eos

    class(bc_periodic) this
    class(eos), intent(in) :: eq
    real, intent(inout) :: sol(:,:,:)
    integer, intent(in) :: nb ! number of boundary cells

    integer ib
    integer nx, ny

    if (.false.) call no_op(this,eq)

    call assert(size(sol,1).eq.nq, "mismatch in size of solution array &
         & first dimension and number of state components")
    nx = size(sol,2)
    ny = size(sol,3)

    forall (ib=1:nb)
       sol(:,:,ib) = sol(:,:,ny-2*nb)
       sol(:,:,ny-ib+1) = sol(:,:,nb+ib)
    end forall
    forall (ib=1:nb)
       sol(:,ib,:) = sol(:,nx-2*nb,:)
       sol(:,nx-ib+1,:) = sol(:,nb+ib,:)
    end forall
  end subroutine apply
end module m_bc_periodic
