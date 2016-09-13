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
     procedure(df_dhardening), deferred :: df_dhardening
  end type plastic_model

  abstract interface
     function yield_f(this, stress, hardening)
       import plastic_model
       intent(in) stress, hardening
       class(plastic_model) this
       real stress(3,3), yield_f, hardening
     end function yield_f

     function df_dstress(this, stress, hardening)
       import plastic_model
       intent(in) stress, hardening
       class(plastic_model) this
       real stress(3,3), df_dstress(3,3), hardening
     end function df_dstress

     function df_dhardening(this, stress, hardening)
       import plastic_model
       intent(in) stress, hardening
       class(plastic_model) this
       real stress(3,3), df_dhardening, hardening
     end function df_dhardening
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
  subroutine plastic_relax(eq, plmodel, S, Fe, kappa, debug_output_p)
    use m_state
    use m_eos
    use m_matutil, only: identity, inv3, det3
    use m_error, only: warn

    logical debug_output_p

    ! elastic and plastic models
    class(eos), intent(in) :: eq
    class(plastic_model), intent(in) :: plmodel

    ! entropy (for equation of state evaluation)
    real, intent(in) :: S

    ! IN: initial trial deformation gradient; OUT: relaxed deformation gradient
    real, intent(inout) :: Fe(3,3), kappa
    real Ftot(3,3), gtot(3,3), ge(3,3), Fp(3,3), Fp_unscaled(3,3)

    real, parameter :: eps=1e-5 ! stopping criterion

    ! loop couters
    integer i,j,k,l,m,n

    real id3(3,3) ! identity matrix

    integer iiter
    integer, parameter :: maxiter=1000

    real dsigma_dge_(3,3,3,3), dsigma_dFp_(3,3,3,3), dge_dFp_(3,3,3,3), yield_f
    real dsigma_dkappa_(3,3), dkappadot_dzeta_, df_dhardening_, dhardening_dkappa_
    real df_dzeta_, df_dsigma_(3,3), dFpdot_dzeta_(3,3), sigma(3,3), sigma0(3,3)
    real kappa0
    real dzeta
    real hardening

    id3 = identity(3)

    ! does not need to be reevaluated in the loop
    sigma = eq%stress(S,Fe,kappa)
    sigma0 = sigma
    kappa0 = kappa

    ge = inv3(Fe)
    Fp = id3
    Fp_unscaled = id3
    Ftot = Fe
    gtot = ge

    if (debug_output_p) print *, sigma0

    hardening = eq%hardening(S,Fe,kappa)
    yield_f = plmodel%yield_f(sigma, hardening)
    ! Negative yield function (+eps) => unyielded, and the function
    ! returns without modifying Fe and kappa.
    if (yield_f.le.eps) return

    ! If we reach this point, the material has yielded, so perform an
    ! iterative relaxation to the yield surface.  The stopping
    ! criterion (tested at the end of the loop) is now whether the
    ! state is *on* the yield surface (+/- eps), not inside it, since
    ! this is where the relaxed state should go, and it is possible a
    ! step of the iterative solver sends the state significantly
    ! inside the yield surface before it has converged.

    dsigma_dge_ = eq%dstress_dg_E(S,Fe,kappa)
    dge_dFp_ = dge_dFp(gtot)

    dsigma_dFp_ = 0
    do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
       dsigma_dFp_(i,j,m,n) = dsigma_dFp_(i,j,m,n) + dsigma_dge_(i,j,k,l) * dge_dFp_(k,l,m,n)
    end do; end do; end do; end do; end do; end do

    dsigma_dkappa_ = eq%dstress_dkappa_E(S,Fe,kappa)

    do iiter=1,maxiter
       dhardening_dkappa_ = eq%dhardening_dkappa(S,Fe,kappa)

       df_dsigma_ = plmodel%df_dstress(sigma,hardening)

       dFpdot_dzeta_ = dFpdot_dzeta(df_dsigma_,Ftot,ge)

       df_dhardening_ = plmodel%df_dhardening(sigma,hardening)

       dkappadot_dzeta_ = -df_dhardening_

       df_dzeta_ = 0
       do i=1,3; do j=1,3
          do k=1,3; do l=1,3
             df_dzeta_ = df_dzeta_ + df_dsigma_(i,j) * dsigma_dFp_(i,j,k,l) * dFpdot_dzeta_(k,l)
          end do; end do
          df_dzeta_ = df_dzeta_ + df_dsigma_(i,j) * dsigma_dkappa_(i,j) * dkappadot_dzeta_
       end do; end do
       df_dzeta_ = df_dzeta_ + df_dhardening_ * dhardening_dkappa_ * dkappadot_dzeta_

       dzeta = -yield_f/df_dzeta_

       if (debug_output_p) then
          write (10,*) iiter, dzeta
          write (11,*) iiter, dFpdot_dzeta_ * dzeta
          write (12,*) iiter, sigma
          write (13,*) iiter, dsigma_dkappa_
          write (14,*) iiter, dkappadot_dzeta_
          write (15,*) iiter, yield_f
          write (16,'(I3,(3E16.7E3))') iiter, reshape(dsigma_dFp_,[3,3,3,3],order=[4,3,2,1])
          write (17,'(I3,(3E16.7E3))') iiter, reshape(dsigma_dge_,[3,3,3,3],order=[4,3,2,1])
          write (18,'(I3,(3E16.7E3))') iiter, reshape(dge_dFp_,[3,3,3,3],order=[4,3,2,1])
       end if

       Fp_unscaled = Fp + dFpdot_dzeta_ * dzeta
       Fp = det3(Fp_unscaled)**(-1.0/3.0) * Fp_unscaled

       kappa = kappa + dkappadot_dzeta_ * dzeta

       sigma = sigma0
       do i=1,3; do j=1,3;
          do k=1,3; do l=1,3
             sigma(i,j) = sigma(i,j) + dsigma_dFp_(i,j,k,l) * (Fp(k,l) - id3(k,l))
          end do; end do;
          sigma(i,j) = sigma(i,j) + dsigma_dkappa_(i,j) * (kappa - kappa0)
       end do; end do

       ge = matmul(Fp,gtot)
       Fe = inv3(ge)
       hardening = eq%hardening(S,Fe,kappa)
       yield_f = plmodel%yield_f(sigma, hardening)

       if (abs(yield_f).le.eps) return
    end do
  end subroutine plastic_relax
end module m_plastic
