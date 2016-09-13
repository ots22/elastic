import numpy as np
from sympy import *
from sympy.utilities.codegen import codegen
from syms import *
from matrix_deriv import *
from itertools import product
from cse import array_to_f90_with_cse

## abstract out an output_expression function

## sometimes end up with ndarrays of Matrix objects, which makes this tricky

# def output_expression(file_prefix, name, expression):
#     with open(file_prefix + '_' + name + '.inc', 'w') as f:
#         # something with an iterator here, and conitional on the expression being a matrix or scalar...
#         indexed_name = name + '(%d,%d)'
#         f.write(printing.fcode(expression.subs(symmetric_pairs), assign_to=indexed_name, source_format='free'))
#         f.write("\n\n")

def eos_to_f90(file_prefix, internal_energy, entropy=None):
    # if entropy is not supplied, it will attempt to call solve on energy
    internal_energy_subs = internal_energy.subs(symmetric_pairs)

    if entropy is None:
        print "eos_to_f90: Attempting to solve the internal energy expression for entropy.  If this appears to hang, consider passing an expression for entropy to this function."
        entropy = sympy.solve(internal_energy - E, s)[0]
        
    entropy_subs = entropy.subs(symmetric_pairs)

    stress = -2 * rho0 * sqrt(I3(Binv)) * Binv * diff_by_matrix(internal_energy, Binv)
    stress += 2 * (rho0 / sqrt(I3(B))) * B * diff_by_matrix(internal_energy, B)
    stress += 2 * (rho0 / sqrt(I3(c))) * F * diff_by_matrix(internal_energy, c) * sympy.transpose(F)
###  this one isn't correct I think:
###    stress += -2 * rho0 * sqrt(I3(Cinv)) * F * diff_by_matrix(internal_energy, Cinv) * F

    # stress as a function of energy, and its kappa-derivative (at constant internal energy)
    #stress_e = stress.subs(s,entropy)
    #dstress_dkappa_e = diff_of_matrix(stress_e, kappa)

    # Substitute all of the strain measures with their expressions in
    # terms of F, so we can take the derivative w.r.t. g
    #stress_g_e = stress_e.subs(strain_g_pairs)

    #dstress_dg_e = np.empty((3,3), dtype=object)
    #for i,j in product(range(3),range(3)):
    #    dstress_dg_e[i,j] = diff_by_matrix(stress_g_e[i,j],g)
    
    #dstress_dkappa_e = dstress_dkappa_e.subs(symmetric_pairs)
    stress = stress.subs(symmetric_pairs)
    #stress_e = stress_e.subs(symmetric_pairs)
    #for i,j in product(range(3),range(3)):
    #    dstress_dg_e[i,j] = dstress_dg_e[i,j].subs(symmetric_pairs)

    with open(file_prefix + '_energy.inc', 'w') as f:
        f.write(printing.fcode(internal_energy_subs, assign_to='e_internal', source_format='free'))

    with open(file_prefix + '_entropy.inc','w') as f:
        f.write(printing.fcode(entropy_subs, assign_to='entropy', source_format='free'))

    with open(file_prefix + '_stress.inc','w') as f:
        for i in range(1,4):
            for j in range(1,4):
                f.write(printing.fcode(stress[i-1,j-1], assign_to="stress(%d,%d)" % (i,j), source_format='free'))
                f.write("\n\n")
                
    #with open(file_prefix + '_dstress_dkappa_e.inc','w') as f:
    #    for i in range(1,4):
    #        for j in range(1,4):
    #            f.write(printing.fcode(dstress_dkappa_e[i-1,j-1], assign_to="dstress_dkappa_e(%d,%d)" % (i,j), source_format='free'))
    #            f.write("\n\n")

    #dstress_dg_e_ndarray = np.empty((3,3,3,3), dtype=object)
    #for i,j,k,l in product(range(3),range(3),range(3),range(3)):
    #    dstress_dg_e_ndarray[i,j,k,l] = dstress_dg_e[i,j][k,l]

#    dstress_dg_e_indexed = [(idx, dstress_dg_e_ijkl) for idx, dstress_dg_e_ijkl in np.ndenumerate(dstress_dg_e_ndarray)]

#    indices = [(idx, dstress_dg_e_ijkl) for idx, dstress_dg_e_ijkl in np.ndenumerate(dstress_dg_e_ndarray)]
    
    #header, body = array_to_f90_with_cse(dstress_dg_e_ndarray, 'dstress_dg_e')
    #with open(file_prefix + '_dstress_dg_e-vars.inc','w') as f:
    #    f.write(header)
    #with open(file_prefix + '_dstress_dg_e.inc','w') as f:
    #    f.write(body)

    # with open(file_prefix + '_dstress_dg_e.inc', 'w') as f:
    #     for i,j,k,l in product(range(3),range(3),range(3),range(3)):
    #         f.write(printing.fcode(dstress_dg_e[i,j][k,l], assign_to="dstress_dg_e(%d,%d,%d,%d)" % (i+1,j+1,k+1,l+1), source_format='free'))
    #         f.write("\n\n")

