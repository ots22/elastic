module m_bc
  use m_configurable

  private
  public bc

  type, abstract, extends(configurable) :: bc
     contains
       procedure(apply), deferred :: apply
  end type bc

  abstract interface
     subroutine apply(this,eq,sol,nb)
       use m_eos
       import bc
       class(bc) this
       class(eos), intent(in) :: eq
       real, intent(inout) :: sol(:,:,:)
       integer, intent(in) :: nb
     end subroutine apply
  end interface
end module m_bc
