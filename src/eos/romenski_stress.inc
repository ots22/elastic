stress(1,1) = -2*Binv11*rho0*((1.0d0/4.0d0)*B0*beta*(Binv22*Binv33 - &
      Binv23**2)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*B0*((2.0d0/3.0d0)*Binv11 - 1.0d0/3.0d0*Binv22 - 1.0d0/ &
      3.0d0*Binv33)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/2.0d0)*K0*(Binv22*Binv33 - Binv23**2)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(Binv22*Binv33 - Binv23**2)*(exp(s/cv) - 1)*( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22) - 2*Binv12*rho0*((1.0d0/2.0d0)*B0* &
      Binv12*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*(-Binv12*Binv33 + Binv13* &
      Binv23)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*K0*(-Binv12*Binv33 + Binv13*Binv23)*((Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      -Binv12*Binv33 + Binv13*Binv23)*(exp(s/cv) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) - 2*Binv13*rho0*((1.0d0/2.0d0)*B0*Binv13*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(Binv12*Binv23 - Binv13*Binv22)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(Binv12*Binv23 - Binv13* &
      Binv22)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(Binv12*Binv23 - Binv13* &
      Binv22)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)

stress(1,2) = -2*Binv11*rho0*((1.0d0/2.0d0)*B0*Binv12*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(-Binv12*Binv33 + Binv13*Binv23)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(-Binv12*Binv33 + Binv13* &
      Binv23)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(-Binv12*Binv33 + Binv13* &
      Binv23)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2*Binv12*rho0*(( &
      1.0d0/4.0d0)*B0*beta*(Binv11*Binv33 - Binv13**2)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 &
      - Binv11*Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23 &
      **2 + (1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22 &
      *Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22) + (1.0d0/2.0d0)*B0*(-1.0d0/3.0d0* &
      Binv11 + (2.0d0/3.0d0)*Binv22 - 1.0d0/3.0d0*Binv33)*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/ &
      2.0d0)*K0*(Binv11*Binv33 - Binv13**2)*((Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      Binv11*Binv33 - Binv13**2)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2 &
      *Binv13*rho0*((1.0d0/2.0d0)*B0*Binv23*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*( &
      -Binv11*Binv23 + Binv12*Binv13)*(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + &
      Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0) &
      *(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) + (1.0d0/2.0d0)*K0*(-Binv11*Binv23 + Binv12*Binv13)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(-Binv11*Binv23 + Binv12*Binv13)*(exp(s/cv) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma) &
      /(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)

stress(1,3) = -2*Binv11*rho0*((1.0d0/2.0d0)*B0*Binv13*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(Binv12*Binv23 - Binv13*Binv22)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(Binv12*Binv23 - Binv13* &
      Binv22)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(Binv12*Binv23 - Binv13* &
      Binv22)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2*Binv12*rho0*(( &
      1.0d0/2.0d0)*B0*Binv23*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*(-Binv11*Binv23 + &
      Binv12*Binv13)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12 &
      **2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*K0*(-Binv11*Binv23 + Binv12*Binv13)*((Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      -Binv11*Binv23 + Binv12*Binv13)*(exp(s/cv) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) - 2*Binv13*rho0*((1.0d0/4.0d0)*B0*beta*(Binv11*Binv22 - &
      Binv12**2)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*B0*(-1.0d0/3.0d0*Binv11 - 1.0d0/3.0d0*Binv22 + (2.0d0/ &
      3.0d0)*Binv33)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12 &
      **2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/2.0d0)*K0*(Binv11*Binv22 - Binv12**2)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(Binv11*Binv22 - Binv12**2)*(exp(s/cv) - 1)*( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)

stress(2,1) = -2*Binv12*rho0*((1.0d0/4.0d0)*B0*beta*(Binv22*Binv33 - &
      Binv23**2)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*B0*((2.0d0/3.0d0)*Binv11 - 1.0d0/3.0d0*Binv22 - 1.0d0/ &
      3.0d0*Binv33)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/2.0d0)*K0*(Binv22*Binv33 - Binv23**2)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(Binv22*Binv33 - Binv23**2)*(exp(s/cv) - 1)*( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22) - 2*Binv22*rho0*((1.0d0/2.0d0)*B0* &
      Binv12*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*(-Binv12*Binv33 + Binv13* &
      Binv23)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*K0*(-Binv12*Binv33 + Binv13*Binv23)*((Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      -Binv12*Binv33 + Binv13*Binv23)*(exp(s/cv) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) - 2*Binv23*rho0*((1.0d0/2.0d0)*B0*Binv13*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(Binv12*Binv23 - Binv13*Binv22)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(Binv12*Binv23 - Binv13* &
      Binv22)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(Binv12*Binv23 - Binv13* &
      Binv22)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)

stress(2,2) = -2*Binv12*rho0*((1.0d0/2.0d0)*B0*Binv12*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(-Binv12*Binv33 + Binv13*Binv23)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(-Binv12*Binv33 + Binv13* &
      Binv23)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(-Binv12*Binv33 + Binv13* &
      Binv23)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2*Binv22*rho0*(( &
      1.0d0/4.0d0)*B0*beta*(Binv11*Binv33 - Binv13**2)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 &
      - Binv11*Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23 &
      **2 + (1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22 &
      *Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22) + (1.0d0/2.0d0)*B0*(-1.0d0/3.0d0* &
      Binv11 + (2.0d0/3.0d0)*Binv22 - 1.0d0/3.0d0*Binv33)*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/ &
      2.0d0)*K0*(Binv11*Binv33 - Binv13**2)*((Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      Binv11*Binv33 - Binv13**2)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2 &
      *Binv23*rho0*((1.0d0/2.0d0)*B0*Binv23*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*( &
      -Binv11*Binv23 + Binv12*Binv13)*(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + &
      Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0) &
      *(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) + (1.0d0/2.0d0)*K0*(-Binv11*Binv23 + Binv12*Binv13)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(-Binv11*Binv23 + Binv12*Binv13)*(exp(s/cv) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma) &
      /(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)

stress(2,3) = -2*Binv12*rho0*((1.0d0/2.0d0)*B0*Binv13*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(Binv12*Binv23 - Binv13*Binv22)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(Binv12*Binv23 - Binv13* &
      Binv22)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(Binv12*Binv23 - Binv13* &
      Binv22)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2*Binv22*rho0*(( &
      1.0d0/2.0d0)*B0*Binv23*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*(-Binv11*Binv23 + &
      Binv12*Binv13)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12 &
      **2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*K0*(-Binv11*Binv23 + Binv12*Binv13)*((Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      -Binv11*Binv23 + Binv12*Binv13)*(exp(s/cv) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) - 2*Binv23*rho0*((1.0d0/4.0d0)*B0*beta*(Binv11*Binv22 - &
      Binv12**2)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*B0*(-1.0d0/3.0d0*Binv11 - 1.0d0/3.0d0*Binv22 + (2.0d0/ &
      3.0d0)*Binv33)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12 &
      **2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/2.0d0)*K0*(Binv11*Binv22 - Binv12**2)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(Binv11*Binv22 - Binv12**2)*(exp(s/cv) - 1)*( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)

stress(3,1) = -2*Binv13*rho0*((1.0d0/4.0d0)*B0*beta*(Binv22*Binv33 - &
      Binv23**2)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*B0*((2.0d0/3.0d0)*Binv11 - 1.0d0/3.0d0*Binv22 - 1.0d0/ &
      3.0d0*Binv33)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/2.0d0)*K0*(Binv22*Binv33 - Binv23**2)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(Binv22*Binv33 - Binv23**2)*(exp(s/cv) - 1)*( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22) - 2*Binv23*rho0*((1.0d0/2.0d0)*B0* &
      Binv12*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*(-Binv12*Binv33 + Binv13* &
      Binv23)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*K0*(-Binv12*Binv33 + Binv13*Binv23)*((Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      -Binv12*Binv33 + Binv13*Binv23)*(exp(s/cv) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) - 2*Binv33*rho0*((1.0d0/2.0d0)*B0*Binv13*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(Binv12*Binv23 - Binv13*Binv22)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(Binv12*Binv23 - Binv13* &
      Binv22)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(Binv12*Binv23 - Binv13* &
      Binv22)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)

stress(3,2) = -2*Binv13*rho0*((1.0d0/2.0d0)*B0*Binv12*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(-Binv12*Binv33 + Binv13*Binv23)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(-Binv12*Binv33 + Binv13* &
      Binv23)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(-Binv12*Binv33 + Binv13* &
      Binv23)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2*Binv23*rho0*(( &
      1.0d0/4.0d0)*B0*beta*(Binv11*Binv33 - Binv13**2)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 &
      - Binv11*Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23 &
      **2 + (1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22 &
      *Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22) + (1.0d0/2.0d0)*B0*(-1.0d0/3.0d0* &
      Binv11 + (2.0d0/3.0d0)*Binv22 - 1.0d0/3.0d0*Binv33)*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/ &
      2.0d0)*K0*(Binv11*Binv33 - Binv13**2)*((Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      Binv11*Binv33 - Binv13**2)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2 &
      *Binv33*rho0*((1.0d0/2.0d0)*B0*Binv23*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*( &
      -Binv11*Binv23 + Binv12*Binv13)*(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + &
      Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0) &
      *(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) + (1.0d0/2.0d0)*K0*(-Binv11*Binv23 + Binv12*Binv13)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(-Binv11*Binv23 + Binv12*Binv13)*(exp(s/cv) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma) &
      /(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)

stress(3,3) = -2*Binv13*rho0*((1.0d0/2.0d0)*B0*Binv13*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)* &
      B0*beta*(Binv12*Binv23 - Binv13*Binv22)*(Binv11*Binv22*Binv33 - &
      Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*beta)*(-Binv11*Binv22 - Binv11* &
      Binv33 + Binv12**2 + Binv13**2 - Binv22*Binv33 + Binv23**2 + ( &
      1.0d0/3.0d0)*(Binv11 + Binv22 + Binv33)**2)/(Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22) + (1.0d0/2.0d0)*K0*(Binv12*Binv23 - Binv13* &
      Binv22)*((Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*alpha) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*alpha)/(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23 &
      **2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2* &
      Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*(Binv12*Binv23 - Binv13* &
      Binv22)*(exp(s/cv) - 1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 &
      - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)** &
      ((1.0d0/2.0d0)*gamma)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22))* &
      sqrt(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) - 2*Binv23*rho0*(( &
      1.0d0/2.0d0)*B0*Binv23*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - &
      Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**(( &
      1.0d0/2.0d0)*beta) + (1.0d0/4.0d0)*B0*beta*(-Binv11*Binv23 + &
      Binv12*Binv13)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12 &
      **2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*K0*(-Binv11*Binv23 + Binv12*Binv13)*((Binv11*Binv22*Binv33 &
      - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - &
      Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha)/(alpha*(Binv11* &
      Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12* &
      Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/2.0d0)*T0*cv*gamma*( &
      -Binv11*Binv23 + Binv12*Binv13)*(exp(s/cv) - 1)*(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22*Binv33 - Binv11* &
      Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2 &
      *Binv22) - 2*Binv33*rho0*((1.0d0/4.0d0)*B0*beta*(Binv11*Binv22 - &
      Binv12**2)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta)*(-Binv11*Binv22 - Binv11*Binv33 + Binv12**2 + Binv13 &
      **2 - Binv22*Binv33 + Binv23**2 + (1.0d0/3.0d0)*(Binv11 + Binv22 &
      + Binv33)**2)/(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12** &
      2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22) + (1.0d0/ &
      2.0d0)*B0*(-1.0d0/3.0d0*Binv11 - 1.0d0/3.0d0*Binv22 + (2.0d0/ &
      3.0d0)*Binv33)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12 &
      **2*Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/ &
      2.0d0)*beta) + (1.0d0/2.0d0)*K0*(Binv11*Binv22 - Binv12**2)*(( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) - &
      1)*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + &
      2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*alpha) &
      /(alpha*(Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2* &
      Binv33 + 2*Binv12*Binv13*Binv23 - Binv13**2*Binv22)) + (1.0d0/ &
      2.0d0)*T0*cv*gamma*(Binv11*Binv22 - Binv12**2)*(exp(s/cv) - 1)*( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22)**((1.0d0/2.0d0)*gamma)/( &
      Binv11*Binv22*Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2* &
      Binv12*Binv13*Binv23 - Binv13**2*Binv22))*sqrt(Binv11*Binv22* &
      Binv33 - Binv11*Binv23**2 - Binv12**2*Binv33 + 2*Binv12*Binv13* &
      Binv23 - Binv13**2*Binv22)

