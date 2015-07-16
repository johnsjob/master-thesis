from __future__ import division
#--------------------------#
import sys
#--------------------------#
import sympy as s
from sympy import cos, sin, pi
from sympy.matrices import Matrix as mat
#=====================================================#
from helperfunctions import matmul_series
sys.path.append("../int/misc-tools/")
import parsingtools as  parse
#--------------------------------------------------------------------#
def diff_mat(M,param):
    diff = lambda y: s.diff(y,param)
    sh = M.shape
    return mat(map(diff, M)).reshape(*sh)
#--------------------------------------------------------------------#
def matmul_series( *args ):
    return reduce(s.Matrix.multiply, args)
#--------------------------------------------------------------------#
def _spin_tensor_diff_mat_ang(R,ang):
    return matmul_series(diff_mat(mat_rot_z(ang),ang), mat_rot_z(ang).T)
#--------------------------------------------------------------------#
def spin_tensor(R,ang):
    dang_dt = s.sympify('d'+str(ang)+'/dt')
    return _spin_tensor_diff_mat_ang(R,ang) * dang_dt
#--------------------------------------------------------------------#
def mat_trans_x( tx ):
    return mat([[1,        0,      0,     tx],
                [0,        1,      0,      0],
                [0,        0,      1,      0],
                [0,        0,      0,      1]])
#--------------------------------------------------------------------#
def mat_trans_z( tz ):
    return mat([[1,     0,     0,      0],
                [0,     1,     0,      0],
                [0,     0,     1,      tz],
                [0,     0,     0,      1]])
#--------------------------------------------------------------------#
def mat_rot_x( ang ):
    #convert to radians
    c = cos(ang)
    s = sin(ang)
    return mat([[1,     0,      0,      0],
                [0,     c,     -s,      0],
                [0,     s,      c,      0],
                [0,     0,      0,      1]])
#--------------------------------------------------------------------#
def mat_rot_z( ang ):
    #convert to radians
    c = cos(ang)
    s = sin(ang)
    return mat([[c,     -s,     0,     0],
                [s,      c,     0,     0],
                [0,      0,     1,     0],
                [0,      0,     0,     1]])
#--------------------------------------------------------------------#
def transform_to_next(A, alpha, D, theta):
    Rz_J = mat_rot_z(theta)
    Tz_J = mat_trans_z(D)
    Tx_I = mat_trans_x(A)
    Rx_I = mat_rot_x(alpha)
    return matmul_series(Rz_J, Tz_J, Tx_I, Rx_I)
#--------------------------------------------------------------------#
def DH_params( *params ):
    nbr_of_sections = int(len(params) / 4)
    if len(params) == 1 and type(params[0]) in [list, tuple]:
        raise ArithmeticError("Function does not use lists or tuples, please unpack using *.")
    elif not (len(params) % 4 == 0):
        raise ArithmeticError("Invalid number of Denavit-Hartenberg parameters.")

    matrices = []
    for k in xrange(0, nbr_of_sections):
        A, alpha, D, theta = params[4*k:4*k+4]
        matrices.append( transform_to_next(A, alpha, D, theta) )
    return matmul_series(*matrices)
#--------------------------------------------------------------------#
def calc_tool_IRB120(a=None,b=None,c=None,d=None,e=None,f=None):
    if a is None:
        a = s.sympify('a')
    if b is None:
        b = s.sympify('b')
    if c is None:
        c = s.sympify('c')
    if d is None:
        d = s.sympify('d')
    if e is None:
        e = s.sympify('e')
    if f is None:
        f = s.sympify('f')
    flange = DH_params(
                        0, 90,0.290,180+a,
                       0.270,0,0,90+b,
                      -0.070, 90, 0, 180+c,
                      0, 90, 0.302, 180+d,
                      0, 90, 0, 180+e,
                      0, 0, 0.072, 0+f
                        )
    return flange
#--------------------------------------------------------------------#
def custom_round(v, prec = 1e-4):
    coef = 1 / prec
    return n.round(v * coef) / coef
#--------------------------------------------------------------------#
if __name__ == '__main__':
    a = s.sympify('a')
    
