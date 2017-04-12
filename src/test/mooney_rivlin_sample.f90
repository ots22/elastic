program mooney_rivlin_sample
  use m_eos_mooney_rivlin
  use m_matutil
  use m_materials, only: example_material_mr

  real C(3,3), rndn(7), S, F(3,3)

  ! infinite stream of samples from the eos
  do
     call random_number(rndn)
     rndn(1:3) = 0.95 + 0.1 * rndn(1:3)
     rndn(4:6) = 0.1 * rndn(4:6)
     rndn(7) = 100 * rndn(7)
     C = voigt_to_matrix(rndn(1:6))
     S = rndn(7)
     
     F = cholesky(C)

     write (*,'(I3,8E18.11)'), 0, rndn, example_material_mr%E(S,F,0.0)

  end do


end program mooney_rivlin_sample
