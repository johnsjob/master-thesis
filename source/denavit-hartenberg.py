from __future__ import division
#--------------------------#
import numpy as n
from numpy import cos, sin, pi
from numpy import array as mat
import operator
#--------------------------#
def prod(factors):
    return reduce(operator.mul, factors, 1)
#--------------------------#
def matmul(matrix_factors):
	return reduce(n.dot, matrix_factors, 1)

def round(v, prec=1e-4):
    factor = 1/prec
    if v > 0:
        return n.int32(v * factor + 0.5)/factor
    elif v < 0:
        return n.int32(v * factor - 0.5)/factor
    return 0
#--------------------------#
def list_round(L):
    sh = L.shape
    L = L.reshape(prod(sh))
    res = mat(map(round, L)).reshape(sh)
    return res
#--------------------------#
def rotZj_theta_i(theta_i):
    c = lambda: cos(theta_i)
    s = lambda: sin(theta_i)
    return mat([[c(), -s(), 0, 0],
                [s(),  c(), 0 , 0],
                [0,      0, 1,  0],
                [0,      0, 0,  1]])    

def rotXi_alpha_i(alpha_i):
    c = lambda: cos(alpha_i)
    s = lambda: sin(alpha_i)
    return mat([[1, 0, 0, 0],
                [0, c(), -s(), 0],
                [0, s(),  c(), 0],
                [0,   0,    0, 1]])

def transZj_di(d_i):
    return mat([[1, 0, 0,    0],
                [0, 1, 0,    0],
                [0, 0, 1,  d_i],
                [0, 0, 0,    1]])    

def transXi_ai(a_i):
    return mat([[1, 0, 0,  a_i],
                [0, 1, 0,    0],
                [0, 0, 1,    0],
                [0, 0, 0,    1]])

def DHji(ai, alphai, di, thetai):
    res = rotZj_theta_i(thetai)
    res = res.dot( transZj_di(di) )
    res = res.dot( transXi_ai(ai) )
    res = res.dot( rotXi_alpha_i(alphai) )
    return res

def DH(a,b,c,d,e,f):
    A01 = DHji(0,   pi/2,   0.290,    pi+a)
    A12 = DHji(0.340, 0,      0,      pi/2+b)
    A23 = DHji(0,   pi/2,   0,      pi+c)
    A34 = DHji(0,   pi/2,   0.302,    pi+d)
    A45 = DHji(0,   pi/2,   0,      pi+e)
    A56 = DHji(0,   0,      0.72,     f)
    matrices = (A01, A12, A23, A34, A45, A56)
    return matmul(matrices)

def pr(s):
    print list_round(s)

if __name__ == '__main__':
    T44 = mat([[0.364179993608069, -0.88395872629214, -0.293240349317084, 0],
                [0.594536079422304, 0.463016221588474, -0.657375713584773, 0],
                [0.716868037033862, 0.0650611155599963, 0.694166600119385, 0],
                [0.0836863220332209, -0.0181462742616406, 0.648626596909796,1]]).T
    j0 = -35.4455740562938
    j1 = -56.7487189266047
    j2 =  38.5964427732583
    j3 = -43.7627466350371
    j4 = -42.6628949333683
    j5 = -0.0244652973799493
    T = DH(j0,j1,j2,j3,j4,j5)
    
