! Sample the energy and gradient (dE/dG, where G is the inverse
! right(?) Cauchy-Green tensor), at points (G, potT) read from STDIN

program eos_gp_gradient

  use m_eos_gp
  use m_matutil

  integer, parameter :: STDERR=0, STDIN=5, STDOUT=6
  integer u, read_status
  type(eos_gp) eq

  real Gvoigt(6), G(3,3), S, grad_voigt(6), stress_voigt(6), F(3,3)

  open(newunit=u, file='gp.eos')

  call eq%init_from_config(u)

  print *, "G11 G22 G33 G23 G31 G12 potT E &
       &dEdG11 dEdG22 dEdG33 dEdG23 dEdG31 dEdG12 &
       &sigma11 sigma22 sigma33 sigma23 sigma31 sigma12"

  do

     read (*,*, iostat=read_status) Gvoigt, S
     
     if (read_status<0) then
        EXIT
     else if (read_status>0) then
        STOP 'BAD INPUT'
     end if

     G = voigt_to_matrix(Gvoigt)
     F = cholesky(inv3(G))
     grad_voigt = matrix_to_voigt(dE_dCinv(eq,S,G))
     stress_voigt = matrix_to_voigt(eq%stress(S,F,0.0))
     write(*,'(20E16.7)') Gvoigt, S, eq%E(S, F, 0.0), grad_voigt, stress_voigt

  end do

end program eos_gp_gradient
