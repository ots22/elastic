from syms import *
from sympy import *
import eos

params = ['K0', 'B0', 'alpha', 'beta', 'gamma', 'cv', 'T0']
list_to_symbol_dict(params, globals())

e_internal = K0/(2*alpha**2) * (I3(Binv)**(alpha/2) - 1)**2                \
    + cv*T0*I3(Binv)**(gamma/2) * (exp(s/cv) - 1)                          \
    + (B0/2)*I3(Binv)**(beta/2) * (I1(Binv)**2/3 - I2(Binv))

E = symbols('E') # (sympy uses E for base of natural logs, and the import clobbers it)
# from axiom (sympy had trouble with the solve)
entropy = cv*log((6*T0*alpha*alpha*cv*I3(Binv)**(gamma/2)+(3*B0*I2(Binv)+(-B0*I1(Binv)*I1(Binv)))*alpha*alpha*I3(Binv)**(beta/2)+(-3*K0*I3(Binv)**alpha)+6*K0*I3(Binv)**(alpha/2)+6*E*alpha*alpha+(-3*K0))/(6*T0*alpha*alpha*cv*I3(Binv)**(gamma/2)))

eos.eos_to_f90('eos/romenski', e_internal, entropy)
