! Riemann-problem initial conditions
module m_ic_RP
  use m_ic

  private
  public ic_RP

  type, extends(ic) :: ic_RP
     real a
     real v_left(3), F_left(3,3), S_left
     real v_right(3), F_right(3,3), S_right
   contains
     procedure apply
     procedure init_from_config
  end type ic_RP

contains
  subroutine apply(this,eq,u)
    use m_state, only: nq, F_density
    use m_eos, only: eos, prim_to_cons
    use m_error, only: assert
    use m_matutil, only: rot_z

    class(ic_RP) this
    class(eos), intent(in) :: eq
    real, intent(out) :: u(:,:,:)

    real a

    integer nx, ny, iy, iq, xoff
    real rho_left, v_left(3), F_left(3,3), S_left
    real rho_right, v_right(3), F_right(3,3), S_right
    real, dimension(nq) :: left_prim_state, right_prim_state
    real, dimension(nq) :: left_cons_state, right_cons_state

    real R(3,3), x, xfrac, xcorr

    a = this%a
    rho_left = F_density(eq%rho0,this%F_left)
    v_left  = this%v_left
    F_left  = this%F_left 
    S_left  = this%S_left 
    rho_right = F_density(eq%rho0,this%F_right)
    v_right = this%v_right
    F_right = this%F_right
    S_right = this%S_right

    call assert(size(u,1).eq.nq, 'size of first dimension of u in &
    &ic_RP%apply must be nq')
    
    nx = size(u,2)
    ny = size(u,3)

!   rotate everything: original 1d test problem, use identity
    R = rot_z(a)

    left_prim_state  = [rho_left, matmul(R,v_left),  matmul(R,F_left),  S_left]
    right_prim_state = [rho_right, matmul(R,v_right), matmul(R,F_right), S_right]

    left_cons_state  = prim_to_cons(eq,left_prim_state)
    right_cons_state = prim_to_cons(eq,right_prim_state)

    do iy=1,ny
       x = nx/2.0 + (0.5+ iy-(1.0*ny)/2.0)*tan(a)
       xoff = floor(x)
       xfrac = x - floor(x)
       do iq=1,nq
          u(iq,1:xoff-1,iy)  = left_cons_state(iq)
          u(iq,xoff+1:nx,iy) = right_cons_state(iq)
       end do

       if (xfrac+tan(a)/2.ge.1.0) then 
          xcorr = (1-(1-xfrac+0.5*tan(a))/tan(a)) * (xfrac + 0.5*tan(a) - 1)
          xfrac = xfrac - xcorr
          do iq=1,nq
             u(iq,xoff+1,iy) = xcorr*left_cons_state(iq) + (1-xcorr)*right_cons_state(iq)
          end do
!         u(:,xoff+1,iy) = prim_to_cons(eq, u(:,xoff+1,iy))
       else if (xfrac-tan(a)/2.le.0.0) then
          xcorr = (1-(xfrac+0.5*tan(a))/tan(a)) * (0.5*tan(a) - xfrac)
          xfrac = xfrac + xcorr
          do iq=1,nq
             u(iq,xoff-1,iy) = (1-xcorr)*left_cons_state(iq) + xcorr*right_cons_state(iq)
          end do
!         u(:,xoff-1,iy) = prim_to_cons(eq, u(:,xoff-1,iy))
       end if

       do iq=1,nq
          u(iq,xoff,iy) = xfrac*left_cons_state(iq) + (1-xfrac)*right_cons_state(iq)
       end do
!      u(:,xoff,iy) = prim_to_cons(eq, u(:,xoff,iy))
    end do

  end subroutine apply

  subroutine init_from_config(this, u)
    use m_matutil, only: identity
    class(ic_RP) this
    integer, intent(in) :: u ! unit number
    real a, v_left(3), F_left(3,3), S_left, v_right(3), F_right(3,3), S_right
    namelist/ic_RP/a,v_left,F_left,S_left,v_right,F_right,S_right

    a = 0
    F_left = identity(3)
    v_left = 0
    S_left = 0

    F_right = identity(3)
    v_right = 0
    S_right = 0

    read (u,nml=ic_RP)  

    this%a = a
    this%v_left  = v_left
    this%F_left  = F_left 
    this%S_left  = S_left 
    this%v_right = v_right
    this%F_right = F_right
    this%S_right = S_right
  end subroutine init_from_config
end module m_ic_RP
