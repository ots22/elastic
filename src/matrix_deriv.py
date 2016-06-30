import numpy as np
import sympy

def diff_of_matrix(M, x):
    return sympy.Matrix(M.rows, M.cols, lambda i,j: sympy.diff(M[i,j], x))

def diff_by_matrix(expr, M):
    return sympy.Matrix(M.rows, M.cols, lambda i,j: sympy.diff(expr, M[i,j]))
