from syms import *
from sympy import *
import eos

params = ['lambda_0', 'lambda_s', 'mu_0', 'mu_s', 'theta_0', 'theta_1', 'kappa']
list_to_symbol_dict(params, globals())

e_internal = 0.5 * (lambda_0 + lambda_s * s) * (log(sqrt(I3(c))))**2      \
    + 0.5 * (mu_0 + mu_s * s) * (I1(c))                                   \
    - 0.5 * (mu_0 + mu_s * s) * (log(I3(c)))                              \
    + (theta_0 * sqrt(I3(c))) * (kappa + exp(-theta_1 * kappa)/theta_1)

eos.eos_to_f90('eos/mooney_rivlin', e_internal)
