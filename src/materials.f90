! defines the equation of state parameters for specific materials
! 'init_materials' must be called before using any of them

module m_materials
  use m_eos_romenski

  private
  public copper_romenski

  real, parameter  :: c0=4.6, b0=2.1
! alpha, beta and gamma are not exact integers to avoid a bug in glibc
  type(eos_romenski), parameter :: copper_romenski = &
       eos_romenski(rho0  = 8.93, &
                    K0    = c0**2 - (4.0/3.0)*b0**2, &
                    B0    = b0**2, &
                    alpha = 1.00001, &
                    beta  = 3.00001, &
                    gamma = 2.00001, &
                    cv    = 3.9E-4, &
                    T0    = 300.0)
end module m_materials
