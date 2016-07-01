import numpy as np
from sympy import *
from sympy.utilities.codegen import codegen
from syms import *
from matrix_deriv import *

def eos_to_f90(file_prefix, internal_energy, entropy=None):
    # if entropy is not supplied, it will attempt to call solve on energy
    internal_energy_subs = internal_energy.subs(symmetric_pairs)

    if entropy is None:
        print "eos_to_f90: Attempting to solve the internal energy expression for entropy.  If this appears to hang, consider passing an expression for entropy to this function."
        entropy = sympy.solve(internal_energy - E, s)[0]
        
    entropy_subs = entropy.subs(symmetric_pairs)

    stress = -2 * rho0 * sqrt(I3(Binv)) * Binv * diff_by_matrix(internal_energy, Binv)
    stress += 2 * (rho0 / sqrt(I3(B))) * B * diff_by_matrix(internal_energy, B)
    stress += 2 * (rho0 / sqrt(I3(c))) * F * diff_by_matrix(internal_energy, c) * F
###  this one isn't correct I think:
###    stress += -2 * rho0 * sqrt(I3(Cinv)) * F * diff_by_matrix(internal_energy, Cinv) * F
    stress = stress.subs(symmetric_pairs)

    with open(file_prefix + '_energy.inc', 'w') as f:
        f.write(printing.fcode(internal_energy_subs, assign_to='e_internal', source_format='free'))

    with open(file_prefix + '_entropy.inc','w') as f:
        f.write(printing.fcode(entropy_subs, assign_to='entropy', source_format='free'))

    with open(file_prefix + '_stress.inc','w') as f:
        for i in range(1,4):
            for j in range(1,4):
                f.write(printing.fcode(stress[i-1,j-1], assign_to="stress(%d,%d)" % (i,j), source_format='free'))
                f.write("\n\n")
