e_internal = theta_0*(kappa + exp(-kappa*theta_1)/theta_1)*sqrt(C11*C22* &
      C33 - C11*C23**2 - C12**2*C33 + 2*C12*C13*C23 - C13**2*C22) + ( &
      0.5d0*lambda_0 + 0.5d0*lambda_s*s)*log(sqrt(C11*C22*C33 - C11*C23 &
      **2 - C12**2*C33 + 2.0d0*C12*C13*C23 - C13**2*C22))**2 + (0.5d0* &
      mu_0 + 0.5d0*mu_s*s)*(C11 + C22 + C33) - (0.5d0*mu_0 + 0.5d0*mu_s &
      *s)*log(C11*C22*C33 - C11*C23**2 - C12**2*C33 + 2.0d0*C12*C13*C23 &
      - C13**2*C22)