from __future__ import division
#--------------------------#
import numpy as N
import sympy as S
from sympy import init_printing
from numpy import pi
#--------------------------#
def is_sympy_matrix(a):
    return 'sympy.matrices' in str(type(a))
def is_numpy_matrix(a):
    return type(a) == N.ndarray
def mul(A, B, mode='sympy'):
    mode = mode.lower()
    if mode == 'sympy':
        if is_sympy_matrix(A) and is_sympy_matrix(B):
            op = A.multiply(B)
        else:
            raise Exception("Conflicting types!")
    elif mode == 'numpy':
        if is_numpy_matrix(A) and is_numpy_matrix(B):
            op = A.dot(B)
        else:
            raise Exception("Conflicting types!")
    else:
        raise Exception("Supported Modes: 'sympy', 'numpy'.")
    return op
#--------------------------#
#j = i-1
#--------------------------#
def matRotZj(oi = None, mode = 'sympy'):
    mode = mode.lower()
    if  mode == 'sympy':
        mat = S.Matrix
        if oi == None:
            oi = S.sympify('\theta_i')
        c = S.cos
        s = S.sin
    elif mode == 'numpy':
        mat = N.array
        if oi == None:
            oi = 0;
        c = N.cos
        s = S.sin
    else:
        raise Exception("Supported Modes: 'sympy', 'numpy'.")

    M = [[ c(oi), -s(oi), 0, 0],
         [ s(oi),  c(oi), 0, 0],
         [ 0,          0, 1, 0],
         [ 0,          0, 0, 1]]
    return mat(M)
#--------------------------#
def matTransZj(di = None, mode = 'sympy'):
    mode = mode.lower()
    if  mode == 'sympy':
        mat = S.Matrix
        if di == None:
            di = S.sympify('d_i')
    elif mode == 'numpy':
        mat = N.array
        if di == None:
            di = 0;
    else:
        raise Exception("Supported Modes: 'sympy', 'numpy'.")

    M = [[ 1, 0, 0, 0],
         [ 0, 1, 0, 0],
         [ 0, 0, 1, di],
         [ 0, 0, 0, 1]]
    return mat(M)
#--------------------------#
def matTransXi(ai = None, mode = 'sympy'):
    mode = mode.lower()
    if  mode == 'sympy':
        mat = S.Matrix
        if ai == None:
            ai = S.sympify('a_i')
    elif mode == 'numpy':
        mat = N.array
        if ai == None:
            ai = 0;
    else:
        raise Exception("Supported Modes: 'sympy', 'numpy'.")

    M = [[ 1, 0, 0, ai],
         [ 0, 1, 0, 0],
         [ 0, 0, 1, 0],
         [ 0, 0, 0, 1]]
    return mat(M)
#--------------------------#
def matRotXi(ai = None, mode = 'sympy'):
    mode = mode.lower()
    if  mode == 'sympy':
        mat = S.Matrix
        if ai == None:
            ai = S.sympify('\alpha_i')
        c = S.cos
        s = S.sin
    elif mode == 'numpy':
        mat = N.array
        if ai == None:
            ai = 0;
        c = N.cos
        s = S.sin
    else:
        raise Exception("Supported Modes: 'sympy', 'numpy'.")

    M = [[  1,          0,      0, 0],
         [  0,      c(ai), -s(ai), 0],
         [ 0,      s(ai),  c(ai), 0],
         [ 0,          0,      0, 1]]
    return mat(M)
#--------------------------#
def DHji(ai= None, alphi = None, di = None, theti = None, mode = 'sympy'):
    if mode.lower() == 'sympy':
        ai = S.sympify(ai)
        alphi = S.sympify(alphi)
        di = S.sympify(di)
        theti = S.sympify(theti)
    left = mul(matRotZj(theti, mode),matTransZj(di,mode), mode)
    right = mul(matTransXi(ai, mode), matRotXi(alphi, mode), mode)
    Ai = mul(left,  right, mode)
    return Ai
#--------------------------#
def DH(a, b, g, z, e, f, _mode = 'sympy'):
    A01 = DHji(ai = 0,      alphi = 0,          di = 0.290, theti = 0       , mode = _mode)
    A12 = DHji(ai = 0,      alphi = pi/2,       di = 0,     theti = pi      , mode = _mode)
    A23 = DHji(ai = 340,    alphi = 0,          di = 0,     theti = pi/2    , mode = _mode)
    A34 = DHji(ai = 0,      alphi = pi/2,       di = 0,     theti = pi    , mode = _mode)
    A45 = DHji(ai = 0,      alphi = pi/2,       di = 302,   theti = pi    , mode = _mode)
    A56 = DHji(ai = 0,      alphi = pi/2,       di = 0,     theti = pi       , mode = _mode)
    A67 = DHji(ai = 0,      alphi = 0,          di = 72,    theti = 0       , mode = _mode)

    #maps from tool-space to world-space
    if _mode.lower() == 'sympy':
        return A01 * A12 * A23 * A34 * A45 * A56
    else:
        return A01.dot(A12).dot(A23).dot(A34).dot(A45).dot(A56)
        
#--------------------------#
def round(a, prec=1e-3):
    prec = 1/prec
    return N.float64(N.int32(a * prec + 0.5) / prec)
#--------------------------#
if __name__ == '__main__':
    mat = N.array
    init_printing()
#============== Initial State =====================#
    T44 = N.array([[0.000648032180188635, 0.0088508747110115, -0.99996062026019, 0], [7.93095000264074E-06, 0.999960830147902, 0.00885088170848774, 0], [0.999999789995675, -1.36662938540362E-05, 0.00064793660101917, 0], [0.374465331936172, -9.83973157376271E-07, 0.629797169829582, 1]])
    T44 = T44.T
    res = DH(0,0,0,0,0,0, 'numpy')
    print "#Initial State#"
    S.pprint(round(res[:,3]))
    #print T44
#============== Random State =====================#
##    print "#Random State#"
##    T44 = mat([[0.364179993608069, -0.88395872629214, -0.293240349317084, 0],
##           [0.594536079422304, 0.463016221588474, -0.657375713584773, 0],
##           [0.716868037033862, 0.0650611155599963, 0.694166600119385, 0],
##           [0.0836863220332209, -0.0181462742616406, 0.648626596909796, 1]])
##    T44 = T44.T
##    res = DH(-35.4455740562938,
##             -56.7487189266047,
##             38.5964427732583,
##             -43.7627466350371,
##             -42.6628949333683,
##             -0.0244652973799493)
##    
