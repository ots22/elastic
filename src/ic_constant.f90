! Constant initial conditions throughout the domain
! Useful for debugging plasticity
module m_ic_constant
  use m_ic
  
  private
  type, public, extends(ic) :: ic_constant
     real v(3), F(3,3), S, kappa
   contains
     procedure apply
     procedure init_from_config
  end type ic_constant

contains
  subroutine apply(this,eq,u)
    use m_state, only: nq, F_density
    use m_eos, only: eos, prim_to_cons
    use m_error, only: assert

    class(ic_constant) this
    class(eos), intent(in) :: eq
    real, intent(out) :: u(:,:,:)

    real rho, F(3,3), v(3), S, kappa
    
    real prim_state(nq), cons_state(nq)

    integer nx, ny

    rho = F_density(eq%rho0, this%F)
    v = this%v
    F = this%F
    S = this%S
    kappa = this%kappa

    call assert(size(u,1).eq.nq, 'size of first dimension of u in &
    &ic_constant%apply must be nq')
    
    nx = size(u,2)
    ny = size(u,3)

    prim_state = [rho, v, F, S, kappa]
    cons_state = prim_to_cons(eq,prim_state)

    u(:,:,:) = spread(spread(cons_state,2,nx),3,ny)
  end subroutine apply

  subroutine init_from_config(this, u)
    use m_matutil, only: identity
    class(ic_constant) this
    integer, intent(in) :: u ! unit number
    real v(3), F(3,3), S, kappa
    namelist/ic_constant/v,F,S,kappa

    F = identity(3)
    v = 0
    S = 0
    kappa = 0

    read (u,nml=ic_constant)

    this%v = v
    this%F = F
    this%S = S
    this%kappa = kappa
  end subroutine init_from_config
end module m_ic_constant
