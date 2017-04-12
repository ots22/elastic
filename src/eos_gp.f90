module m_eos_gp
  use m_eos
  use m_GP_dense
  use m_error
  
  public :: eos_gp, dE_dCinv

  private

  type, extends(eos) :: eos_gp
!    the underlying Gaussian process
     type(DenseGP) :: gp
!    the mininum and maxium entropy parameter in the database, used to
!    bound the solution when solving for S
     real Smin, Smax
!    convert the energy given from the Gaussian process to macroscopic
!    units used in the simulation:
!    eV/(simulation cell mass) ---> kJ/g
     real energy_unit
   contains
     procedure E
     procedure S
     procedure stress
     procedure init_from_config
  end type eos_gp

contains
  function E(this, S, F, kappa) result(result)
    use m_state, only: Cauchy_Green_right
    use m_matutil, only: matrix_to_voigt, inv3
    use m_error
    class(eos_gp) this
    real, intent(in) :: S, F(3,3), kappa
    real result
    real Cinv_voigt(6), optS
    Cinv_voigt = matrix_to_voigt(inv3(Cauchy_Green_right(F)))
    optS = S
!    optS = 0.0

!    write (17,*), "F", F
!    write (17,*), "Cinv", Cinv_voigt
!    write (17,*), this%energy_unit
!    write (17,*), "S", S 

    result = this%energy_unit * (optS + this%gp%predict([Cinv_voigt, S], 0))
    !result = 0.0
  end function E

  function S(this, E, F, kappa)
    use m_bounded_newton
    use m_state, only: Cauchy_Green_right
    use m_matutil, only: matrix_to_voigt, inv3
    class(eos_gp) this
    real, intent(in) :: E, F(3,3), kappa
    real S, bounds(2), input(7)
    integer err
    input(1:6) = matrix_to_voigt(inv3(Cauchy_Green_right(F)))
    bounds(:) = [this%Smin, this%Smax]
    S = 0.5*(this%Smin + this%Smax) ! initial guess
!    print *, E
    call bounded_newton(E0, dE0, bounds(1), bounds(2), S, maxit_arg=20, err=err)
    if (err.ne.NR_CONVERGED) then
       write (0,*) "  return code: ", err
       write (0,*) "  bounds: ", bounds
       write (0,*) "  Emin and Emax: ", E0(bounds(1)), E0(bounds(2))
       write (0,*) "  S: ", S
       write (0,*) "  E: ", E
       write (0,*) "  inverse left Cauchy-Green tensor (Voigt order): ", input(1:6)
       call panic("bounded_newton failed to find entropy!")
    end if
  contains
!   helper functions for bounded_newton
    real function E0(x)
      real, intent(in) :: x
      input(7) = x
      E0 = this%energy_unit * (x + this%gp%predict(input, 0)) - E
    end function E0
    real function dE0(x)
      real, intent(in) :: x
      input(7) = x
      dE0 = this%energy_unit * (1 + this%gp%predict(input, 7))
    end function dE0
  end function S

! helper routine for stress (see also dCinv_dF in m_state)
  function dE_dCinv(this,S,Cinv)
    use m_matutil, only: matrix_to_voigt, voigt_to_matrix
    real dE_dCinv(3,3), dE_dCinv_voigt(6), dE_dCinv_voigt_num(6)
    class(eos_gp) this
    real, intent(in) :: S, Cinv(3,3)
    real eps, h(6)
    integer i
    eps = 1e-5
    do i=1,6
       dE_dCinv_voigt(i) = this%energy_unit * this%gp%predict([matrix_to_voigt(Cinv),S], i)

 !      h = 0
 !      h(i) = eps
 !      dE_dCinv_voigt_num(i) = this%energy_unit * (this%gp%predict([matrix_to_voigt(Cinv)+h,S], 0) &
 !           - this%gp%predict([matrix_to_voigt(Cinv)-h,S], 0)) / (2*eps)
    end do

    ! correct for derivative being wrt symmetric components
    dE_dCinv_voigt(4:6) = dE_dCinv_voigt(4:6)/2
!    dE_dCinv_voigt_num(4:6) = dE_dCinv_voigt_num(4:6)/2

!    write (16,*) dE_dCinv_voigt_num - dE_dCinv_voigt

    dE_dCinv = voigt_to_matrix(dE_dCinv_voigt)
  end function dE_dCinv

  function stress(this, S, F, kappa)
    use m_state, only: Cauchy_Green_right, dCinv_dF
    use m_matutil, only: inv3, det3
    real stress(3,3), stress_alt(3,3)
    class(eos_gp) this
    real, intent(in) :: S, F(3,3), kappa
    real dEdCinv(3,3), dEdF_T(3,3), F_inv(3,3), dCinvdF(3,3)
    integer j,k
    dEdCinv = dE_dCinv(this, S, inv3(Cauchy_Green_right(F)))
    F_inv = inv3(F)

    dEdF_T = 0
    do j=1,3
       do k=1,3
          dCinvdF(:,:) = dCinv_dF(j,k,F_inv)
          dEdF_T(k,j) = dEdF_T(k,j) + dot_product(dEdCinv(1,:), dCinvdF(1,:))
          dEdF_T(k,j) = dEdF_T(k,j) + dot_product(dEdCinv(2,:), dCinvdF(2,:))
          dEdF_T(k,j) = dEdF_T(k,j) + dot_product(dEdCinv(3,:), dCinvdF(3,:))
       end do
    end do
    stress = (this%rho0/det3(F)) * matmul(F,dEdF_T)

!    stress_alt = -2.0 * (this%rho0/det3(F)) * matmul(transpose(F_inv),matmul(dEdCinv,F_inv))

!    write (18,"(3F12.6)") stress - stress_alt

  end function stress

  subroutine init_from_config(this,u)
    integer, parameter :: MAX_FILENAME_LEN=1000
    class(eos_gp) this
    integer, intent(in) :: u ! unit number
    character(len=MAX_FILENAME_LEN) gp_file

    real rho0,Smin,Smax,energy_unit

    namelist/EOS_GP/rho0,Smin,Smax,energy_unit,gp_file

    rho0 = 0
    Smin = 1
    Smax = 0
    energy_unit = 0
    gp_file = ''

    read (u,nml=EOS_GP)

    call assert(rho0.gt.0, "eos_gp: density must be positive")
    call assert(Smin.le.Smax, "eos_gp: Smin and Smax must both be specified, and Smin must be smaller than Smax")
    call assert(energy_unit.ne.0, "eos_gp: energy_unit must be set and non-zero")  

    this%rho0 = rho0
    this%Smin = Smin
    this%Smax = Smax
    this%energy_unit = energy_unit
  
    if (gp_file.eq.'') call panic("eos_gp: gp_file must be specified")

    this%gp = DenseGP(gp_file)
  end subroutine init_from_config
end module m_eos_gp
