e_internal = (1.0d0/2.0d0)*B0*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + &
      Binv13**2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + &
      Binv22 + Binv33)**2) + (1.0d0/2.0d0)*K0*((Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)**2/alpha**2 + T0*cv &
      *(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*gamma)