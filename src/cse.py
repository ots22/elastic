import StringIO
import numpy as np
import sympy
from syms import *
from matrix_deriv import *
from itertools import product

# returns a string of temporary variable declarations, and the array
def array_to_f90_with_cse(xs, arr_name):
    indices, values = zip(*[(idx, x) for idx, x in np.ndenumerate(xs)])
    indices = list(indices)
    values = list(values)
    
    print 'about to do cse'
    
    symbols, expressions = sympy.cse(values)

    print 'finished cse'

    header = StringIO.StringIO()
    body = StringIO.StringIO()

    for symb, assignment in symbols:
        print >>header, 'real ', symb

    for symb, assignment in symbols:
        print >>body, symb, '=', assignment

    for idx, expr in zip(indices,expressions):
        fort_idx = tuple(map(lambda x:x+1,idx))
        print >>body, sympy.printing.fcode(expr, assign_to=arr_name+str(fort_idx), source_format='free')

    header_str = header.getvalue()
    header.close()

    body_str = body.getvalue()
    body.close()

    return (header_str, body_str)

