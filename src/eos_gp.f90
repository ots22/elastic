module m_eos_gp
  use m_eos
  use m_SparseGP
  use m_error
  
  public :: eos_gp

  private

  type, extends(eos) :: eos_gp
!    the underlying Gaussian process
     type(SparseGP) :: gp
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
  function E(this, S, F) result(result)
    use m_state, only: Cauchy_C
    use m_matutil, only: matrix_to_voigt
    use m_error
    class(eos_gp) this
    real, intent(in) :: S, F(3,3)
    real result
    real C_voigt(6)
    C_voigt = matrix_to_voigt(Cauchy_C(F))
    result = this%energy_unit * (S + SparseGP_predict(this%gp, [C_voigt, S], 0))
  end function E

  function S(this, E, F)
    use m_bounded_newton
    use m_state, only: Cauchy_C
    use m_matutil, only: matrix_to_voigt
    class(eos_gp) this
    real, intent(in) :: E, F(3,3)
    real S, bounds(2), input(7)
    integer err
    input(1:6) = matrix_to_voigt(Cauchy_C(F))
    bounds(:) = [this%Smin, this%Smax]
    S = 0.5*(this%Smin + this%Smax) ! initial guess
    call bounded_newton(E0, dE0, bounds(1), bounds(2), S, err=err)
    if (err.ne.NR_CONVERGED) then
       write (0,*) "  return code: ", err
       write (0,*) "  bounds: ", bounds
       write (0,*) "  Emin and Emax: ", E0(bounds(1)), E0(bounds(2))
       write (0,*) "  S: ", S
       write (0,*) "  Cauchy C tensor (Voigt order): ", input(1:6)
       call panic("bounded_newton failed to find entropy!")
    end if
  contains
!   helper functions for bounded_newton
    real function E0(x)
      real, intent(in) :: x
      input(7) = x
      E0 = this%energy_unit * (x + SparseGP_predict(this%gp, input, 0))
    end function E0
    real function dE0(x)
      real, intent(in) :: x
      input(7) = x
      dE0 = this%energy_unit * SparseGP_predict(this%gp, input, 7)
    end function dE0
  end function S

! helper routine for stress (see also dC_dF in m_state)
  function dE_dC(this,S,C)
    use m_matutil, only: matrix_to_voigt, voigt_to_matrix
    real dE_dC(3,3), dE_dCv(6)
    class(eos_gp) this
    real, intent(in) :: S, C(3,3)
    integer i
    forall (i=1:6) dE_dCv(i) = this%energy_unit * &
         & SparseGP_predict(this%gp, [matrix_to_voigt(C),S], i)
    forall (i=4:6) dE_dCv(i) = dE_dCv(i)/2
    dE_dC = voigt_to_matrix(dE_dCv)
  end function dE_dC

  function stress(this, S, F)
    use m_state, only: Cauchy_C, dC_dF
    use m_matutil, only: inv3, det3
    real stress(3,3)
    class(eos_gp) this
    real, intent(in) :: S, F(3,3)
    real dEdC(3,3), dEdF_T(3,3), F_inv(3,3), dCdF(3,3)
    integer j,k
    dEdC = dE_dC(this, S, Cauchy_C(F))
    F_inv = inv3(F)

    dEdF_T = 0
    do j=1,3
       do k=1,3
          dCdF(:,:) = dC_dF(j,k,F_inv)
          dEdF_T(k,j) = dEdF_T(k,j) + dot_product(dEdC(1,:), dCdF(1,:))
          dEdF_T(k,j) = dEdF_T(k,j) + dot_product(dEdC(2,:), dCdF(2,:))
          dEdF_T(k,j) = dEdF_T(k,j) + dot_product(dEdC(3,:), dCdF(3,:))
       end do
    end do
    stress = (this%rho0/det3(F)) * matmul(F,dEdF_T)
  end function stress

  subroutine init_from_config(this,u)
    use m_config, only: MAX_FILENAME_LEN
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

    this%gp = SparseGP(gp_file)
  end subroutine init_from_config
end module m_eos_gp
