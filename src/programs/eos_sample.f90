program eos_sample
  use m_matutil
  use m_eos
  use m_eos_mooney_rivlin
  use m_eos_romenski
  use m_materials

!  type(eos_mooney_rivlin) eq
  type(eos_romenski) eq
  real :: Cinv(3,3), Cinv_voigt(6), F(3,3), S, E, dEdCinv(3,3), dEdCinv_voigt(6)
  integer i, j

  eq = copper_romenski

  do
     ! generate random input for the equation of state
     call random_number(S)
     S = S * 1E-4
     call random_number(Cinv_voigt)
     Cinv_voigt(1:3) = 1.0 + 0.2 * (Cinv_voigt(1:3) - 0.5)
     Cinv_voigt(4:6) = 0.1 * (Cinv_voigt(4:6) - 0.5)
     Cinv = voigt_to_matrix(Cinv_voigt)
     F = inv3(cholesky(Cinv))
     
     ! internal energy
     E = eq%E(S, F, 0.0)
     
     ! partial derivative of internal energy w.r.t components of Cinv
     stress = eq%stress(S, F, 0.0)

     ! TODO: need to multiply this by the density factor
     dEdCinv = matmul(transpose(F), matmul(stress, F))
     dEdCinv_voigt = matrix_to_voigt(dEdCinv)
     
     ! multiply off diagonal elts by two (derivative is wrt the
     ! symmetric component)
     dEdCinv_voigt(4:6) = 2 * dEdCinv_voigt(4:6)
     
     write (*,'(7E16.7,I3,E16.7)') Cinv_voigt, S, 0, E
     do i=1,6
        write (*,'(7E16.7,I3,E16.7)') Cinv_voigt, S, i, dEdCinv_voigt(i)
     end do
  end do

end program eos_sample
