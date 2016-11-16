! 2-d top-hat initial conditions (i.e. square pulse)
module m_ic_tophat
  use m_ic
  use m_state, only: nq

  private
  
  type, public, extends(ic) :: ic_tophat
     ! coordinates of the box to set inside (to state1), ambient outside (state0)
     real lower_left(2), upper_right(2)
     real, dimension(nq) :: pstate0, pstate1
   contains
     procedure apply
     procedure init_from_config
  end type ic_tophat

contains
  subroutine apply(this,eq,u)
    use m_state
    use m_eos, only: eos, prim_to_cons
    use m_domain
    use m_matutil, only: pi, identity

    class(ic_tophat) :: this
    class(eos), intent(in) :: eq
    real, intent(out) :: u(:,:,:)
    real r(2) ! cell coordinate
    integer iq,ix,iy

    real, dimension(nq) :: state0, state1

    state0 = prim_to_cons(eq, this%pstate0)
    state1 = prim_to_cons(eq, this%pstate1)

    ! whole domain first set to ambient state (state0)
    forall (iq=1:nq) u(iq,:,:) = state0(iq)

    ! if inside, overwrite with state1
    do ix=1, nx
       do iy=1, ny
          r = [(ix-nb-1)*dx, (iy-nb-1)*dx] * domain_length
          if (all(r.gt.this%lower_left.and.r.le.this%upper_right)) u(:,ix,iy) = state1
       end do
    end do
  end subroutine apply

  subroutine init_from_config(this, u)
    use m_matutil, only: identity
    class(ic_tophat) this
    integer, intent(in) :: u ! unit number
    real lower_left(2), upper_right(2)
    real rho0,v0(3),F0(3,3),S0
    real rho1,v1(3),F1(3,3),S1

    namelist/ic_tophat/lower_left,upper_right,rho0,v0,F0,S0,rho1,v1,F1,S1

    ! lower_left > upper_right => no states set
    this%lower_left = 1
    this%upper_right = 0

    ! set them to something obviously invalid by default
    this%pstate0 = -huge(0.0)
    this%pstate1 = -huge(0.0)

    read (u,nml=ic_tophat)
    
    this%lower_left = lower_left
    this%upper_right = upper_right

    this%pstate0 = [rho0,v0,F0,S0]
    this%pstate1 = [rho1,v1,F1,S1]

  end subroutine init_from_config
end module m_ic_tophat
