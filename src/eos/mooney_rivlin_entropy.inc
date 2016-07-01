entropy = (-C11*mu_0*theta_1*exp(kappa*theta_1) - C22*mu_0*theta_1*exp( &
      kappa*theta_1) - C33*mu_0*theta_1*exp(kappa*theta_1) + 2.0d0*E* &
      theta_1*exp(kappa*theta_1) - 2.0d0*kappa*theta_0*theta_1*sqrt(C11 &
      *C22*C33 - C11*C23**2 - C12**2*C33 + 2*C12*C13*C23 - C13**2*C22)* &
      exp(kappa*theta_1) - lambda_0*theta_1*exp(kappa*theta_1)*log(sqrt &
      (C11*C22*C33 - C11*C23**2 - C12**2*C33 + 2.0d0*C12*C13*C23 - C13 &
      **2*C22))**2 + mu_0*theta_1*exp(kappa*theta_1)*log(C11*C22*C33 - &
      C11*C23**2 - C12**2*C33 + 2.0d0*C12*C13*C23 - C13**2*C22) - 2.0d0 &
      *theta_0*sqrt(C11*C22*C33 - C11*C23**2 - C12**2*C33 + 2*C12*C13* &
      C23 - C13**2*C22))*exp(-kappa*theta_1)/(theta_1*(C11*mu_s + C22* &
      mu_s + C33*mu_s + lambda_s*log(sqrt(C11*C22*C33 - C11*C23**2 - &
      C12**2*C33 + 2.0d0*C12*C13*C23 - C13**2*C22))**2 - mu_s*log(C11* &
      C22*C33 - C11*C23**2 - C12**2*C33 + 2.0d0*C12*C13*C23 - C13**2* &
      C22)))