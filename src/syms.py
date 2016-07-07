import numpy as np
import sympy

def string_to_symbol_dict(params,d):
    for p in params.split(): d[p] = sympy.symbols(p)

def list_to_symbol_dict(params,d):
    for p in params: d[p] = sympy.symbols(p)

##### energy, entropy and initial density

#E, s, rho0 = sympy.symbols('E, s, rho0')
string_to_symbol_dict('E s rho0 kappa',globals())

##### F

F_names_full = 'F11 F12 F13 F21 F22 F23 F31 F32 F33'
string_to_symbol_dict(F_names_full,globals())

F = sympy.Matrix([[F11, F12, F13],
                 [F21, F22, F23],
                 [F31, F32, F33]])

##### B

B_names_full = 'B11 B12 B13 B21 B22 B23 B31 B32 B33'
string_to_symbol_dict(B_names_full,globals())

B = sympy.Matrix([[B11, B12, B13],
                 [B21, B22, B23],
                 [B31, B32, B33]])

B_symmetric_pairs = [(B21, B12), (B31, B13), (B32, B23)]

##### C

c_names_full = 'C11 C12 C13 C21 C22 C23 C31 C32 C33'
string_to_symbol_dict(c_names_full,globals())

c = sympy.Matrix([[C11, C12, C13],
                  [C21, C22, C23],
                  [C31, C32, C33]])

c_symmetric_pairs = [(C21, C12), (C31, C13), (C32, C23)]

##### Binv

Binv_names_full = 'Binv11 Binv12 Binv13 Binv21 Binv22 Binv23 Binv31 Binv32 Binv33'
string_to_symbol_dict(Binv_names_full,globals())

Binv = sympy.Matrix([[Binv11, Binv12, Binv13],
                     [Binv21, Binv22, Binv23],
                     [Binv31, Binv32, Binv33]])

Binv_symmetric_pairs = [(Binv21, Binv12), (Binv31, Binv13), (Binv32, Binv23)]

##### Cinv

Cinv_names_full = 'Cinv11 Cinv12 Cinv13 Cinv21 Cinv22 Cinv23 Cinv31 Cinv32 Cinv33'
string_to_symbol_dict(Cinv_names_full,globals())

Cinv = sympy.Matrix([[Cinv11, Cinv12, Cinv13],
                     [Cinv21, Cinv22, Cinv23],
                     [Cinv31, Cinv32, Cinv33]])

Cinv_symmetric_pairs = [(Cinv21, Cinv12), (Cinv31, Cinv13), (Cinv32, Cinv23)]

##### all symmetric pairs
symmetric_pairs = B_symmetric_pairs + c_symmetric_pairs + Binv_symmetric_pairs + Cinv_symmetric_pairs

##### invariants

def I1(M):
    return sympy.trace(M)

def I2(M):
    return sympy.simplify((sympy.trace(M)**2 - sympy.trace(M**2))/2)

def I3(M):
    return sympy.det(M)

