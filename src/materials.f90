! defines the equation of state parameters for specific materials
! 'init_materials' must be called before using any of them

module m_materials
  use m_eos_romenski
  use m_eos_mooney_rivlin

  private
  public copper_romenski, example_material_mr

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

  type(eos_mooney_rivlin), parameter :: example_material_mr = &
       eos_mooney_rivlin(rho0 = 1.0, &
                         mu_0 = 0.6, &
                         lambda_0 = 0.6, &
                         mu_s = 0.01, &
                         lambda_s = 0.01, &
                         theta_0 = 0.1, &
                         theta_1 = 10.0)
end module m_materials
