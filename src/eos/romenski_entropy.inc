entropy = cv*log(0.166666666666667d0*(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2.0d0*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**(-0.5d0*gamma)*(6.0d0*E*alpha**2 + 6.0d0*K0*( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2.0d0*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(0.5d0*alpha) - &
      3.0d0*K0*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2.0d0*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**alpha - &
      3.0d0*K0 + 6.0d0*T0*alpha**2*cv*(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2.0d0*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**(0.5d0*gamma) + alpha**2*(-B0*(Binv11 + Binv22 &
      + Binv33)**2 + 3.0d0*B0*(Binv11*Binv22 + Binv11*Binv33 - Binv12** &
      2 - Binv13**2 + Binv22*Binv33 - Binv23**2))*(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2.0d0*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**(0.5d0*beta))/(T0*alpha**2*cv))