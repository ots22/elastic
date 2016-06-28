module m_plastic
  use m_configurable

  private

  public plastic_relax, plastic_model

! NB 2016-06-17: ignoring hardening parameters for now

! plasticity using an associated flow rule
  type, abstract, extends(configurable) :: plastic_model
   contains
     procedure(yield_f), deferred :: yield_f
     procedure(df_dstress), deferred :: df_dstress
  end type plastic_model

  abstract interface
     function yield_f(this, stress)
       import plastic_model
       intent(in) stress
       class(plastic_model) this
       real stress(3,3), yield_f
     end function yield_f

     function df_dstress(this, stress)
       import plastic_model
       intent(in) stress
       class(plastic_model) this
       real stress(3,3), df_dstress(3,3)
     end function df_dstress
  end interface

contains
  function dge_dFp(g) ! holding g constant
    real g(3,3), dge_dFp(3,3,3,3)
    integer i,j,l
    dge_dFp = 0
    do i=1,3; do j=1,3; do l=1,3
       dge_dFp(i,j,i,l) = g(l,j)     ! k == i
    end do; end do; end do
  end function dge_dFp

  ! The flow rule.  Both F and g provided, since we have computed both when this function is used.
  function dFpdot_dzeta(df_dsigma_,F,ge)
    real df_dsigma_(3,3), F(3,3), ge(3,3), dFpdot_dzeta(3,3)
    dFpdot_dzeta = matmul(matmul(ge,df_dsigma_),F)
  end function dFpdot_dzeta

  ! If the state is determined to have yielded, return a state relaxed
  ! back onto the yield surface according to maximum dissipation;
  ! otherwise, return the state unchanged.  Fe is the trial state.
  ! Add parameters for accumulating the change in plastic deformation
  ! and an array of hardening parameters etc.
  subroutine plastic_relax(eq, plmodel, S, Fe)
    use m_state
    use m_eos
    use m_matutil, only: identity, inv3, det3
    use m_error, only: warn

    ! elastic and plastic models
    class(eos), intent(in) :: eq
    class(plastic_model), intent(in) :: plmodel

    ! entropy (for equation of state evaluation)
    real, intent(in) :: S

    ! IN: initial trial deformation gradient; OUT: relaxed deformation gradient
    real, intent(inout) :: Fe(3,3)
    real Ftot(3,3), gtot(3,3), ge(3,3), Fp(3,3), Fp_unscaled(3,3)

    real, parameter :: eps=1e-6 ! stopping criterion

    ! loop couters
    integer i,j,k,l,m,n

    real id3(3,3) ! identity matrix

    integer iiter
    integer, parameter :: maxiter=10

    real dsigma_dge_(3,3,3,3), dsigma_dFp_(3,3,3,3), dge_dFp_(3,3,3,3), yield_f
    real df_dzeta_, df_dsigma_(3,3), dFpdot_dzeta_(3,3), sigma(3,3), sigma0(3,3)
    real dzeta

    id3 = identity(3)

    ! does not need to be reevaluated in the loop
    sigma = eq%stress(S,Fe)
    sigma0 = sigma

    ge = inv3(Fe)
    Fp = id3
    Fp_unscaled = id3
    Ftot = Fe
    gtot = ge

    do iiter=1,maxiter
       ! get F and S arguments.  Sort out F, Fe, g, ge, Ftot etc.
       ! Fp can be assumed to be identity at the start

       yield_f = plmodel%yield_f(sigma)

       if (yield_f.le.eps) return

       if (iiter.eq.1) then
          dsigma_dge_ = eq%dstress_dg(S,Fe)
          dge_dFp_ = dge_dFp(gtot)

          dsigma_dFp_ = 0
          do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
             dsigma_dFp_(i,j,m,n) = dsigma_dFp_(i,j,m,n) + dsigma_dge_(i,j,k,l) * dge_dFp_(k,l,m,n)
          end do; end do; end do; end do; end do; end do
       end if

       df_dsigma_ = plmodel%df_dstress(sigma)

       dFpdot_dzeta_ = dFpdot_dzeta(df_dsigma_,Ftot,ge)

       df_dzeta_ = 0
       do i=1,3; do j=1,3; do k=1,3; do l=1,3
          df_dzeta_ = df_dzeta_ + df_dsigma_(i,j) * dsigma_dFp_(i,j,k,l) * dFpdot_dzeta_(k,l)
       end do; end do; end do; end do

       dzeta = -yield_f/df_dzeta_

       Fp_unscaled = Fp + dFpdot_dzeta_ * dzeta
       Fp = det3(Fp_unscaled)**(-1.0/3.0) * Fp_unscaled

       sigma = sigma0
       do i=1,3; do j=1,3; do k=1,3; do l=1,3
          sigma(i,j) = sigma(i,j) + dsigma_dFp_(i,j,k,l) * (Fp(k,l) - id3(k,l))
       end do; end do; end do; end do

       ge = matmul(Fp,gtot)
       Fe = inv3(ge)
    end do
!    call warn('reached maximum iterations in plastic_relax. Continuing.')
  end subroutine plastic_relax
end module m_plastic
