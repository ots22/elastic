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

##### g (= F^{-1})

g_names_full = 'g11 g12 g13 g21 g22 g23 g31 g32 g33'
string_to_symbol_dict(g_names_full,globals())

g = sympy.Matrix([[g11, g12, g13],
                 [g21, g22, g23],
                 [g31, g32, g33]])


##### express these in terms of each other

g_F = F.inv()

F_g = g.inv()

F_g_pairs = [(Fij, F_g[idx]) for idx, Fij in np.ndenumerate(F)]

##### B = F * F^T

B_names_full = 'B11 B12 B13 B21 B22 B23 B31 B32 B33'
string_to_symbol_dict(B_names_full,globals())

B = sympy.Matrix([[B11, B12, B13],
                 [B21, B22, B23],
                 [B31, B32, B33]])

B_symmetric_pairs = [(B21, B12), (B31, B13), (B32, B23)]

B_F = sympy.simplify(F * sympy.transpose(F))

B_g = sympy.simplify(g.inv() * sympy.transpose(g.inv()))

B_g_pairs = [(Bij, B_g[idx]) for idx, Bij in np.ndenumerate(B)]

##### C = F^T * F

c_names_full = 'C11 C12 C13 C21 C22 C23 C31 C32 C33'
string_to_symbol_dict(c_names_full,globals())

c = sympy.Matrix([[C11, C12, C13],
                  [C21, C22, C23],
                  [C31, C32, C33]])

c_symmetric_pairs = [(C21, C12), (C31, C13), (C32, C23)]

c_F = sympy.simplify(sympy.transpose(F), F)

c_g = sympy.simplify(sympy.transpose(g.inv()), g.inv())

c_g_pairs = [(cij, c_g[idx]) for idx, cij in np.ndenumerate(c)]

##### Binv

Binv_names_full = 'Binv11 Binv12 Binv13 Binv21 Binv22 Binv23 Binv31 Binv32 Binv33'
string_to_symbol_dict(Binv_names_full,globals())

Binv = sympy.Matrix([[Binv11, Binv12, Binv13],
                     [Binv21, Binv22, Binv23],
                     [Binv31, Binv32, Binv33]])

Binv_symmetric_pairs = [(Binv21, Binv12), (Binv31, Binv13), (Binv32, Binv23)]

Binv_F = sympy.simplify(B_F.inv())

Binv_g = sympy.simplify(sympy.transpose(g) * g)

Binv_g_pairs = [(Binv_ij, Binv_g[idx]) for idx, Binv_ij in np.ndenumerate(Binv)]

##### Cinv

Cinv_names_full = 'Cinv11 Cinv12 Cinv13 Cinv21 Cinv22 Cinv23 Cinv31 Cinv32 Cinv33'
string_to_symbol_dict(Cinv_names_full,globals())

Cinv = sympy.Matrix([[Cinv11, Cinv12, Cinv13],
                     [Cinv21, Cinv22, Cinv23],
                     [Cinv31, Cinv32, Cinv33]])

Cinv_symmetric_pairs = [(Cinv21, Cinv12), (Cinv31, Cinv13), (Cinv32, Cinv23)]

Cinv_F = sympy.simplify(c_F.inv())

Cinv_g = sympy.simplify(g * sympy.transpose(g))

Cinv_g_pairs = [(Cinv_ij, Cinv_g[idx]) for idx, Cinv_ij in np.ndenumerate(Cinv)]

##### all symmetric pairs
symmetric_pairs = B_symmetric_pairs + c_symmetric_pairs + Binv_symmetric_pairs + Cinv_symmetric_pairs

##### all strains --> g
strain_g_pairs = F_g_pairs + B_g_pairs + c_g_pairs + Binv_g_pairs + Cinv_g_pairs

##### invariants

def I1(M):
    return sympy.trace(M)

def I2(M):
    return sympy.simplify((sympy.trace(M)**2 - sympy.trace(M**2))/2)

def I3(M):
    return sympy.det(M)

