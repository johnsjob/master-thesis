from __future__ import division
#--------------------------#
import numpy as n
from numpy import cos, sin, pi
from numpy import array as mat
#--------------------------#
import operator
#--------------------------#
def prod(factors):
    return reduce(operator.mul, factors, 1)
#--------------------------#
def matmul(matrix_factors):
	return reduce(n.dot, matrix_factors, 1)
#=====================================================#
import sys
sys.path.append("../int/misc-tools/")
import parsingtools as  parse
#----------------------------------------------------------------------------------------------------------#
def mat_trans_x( tx ):
    return mat([[1,        0,      0,     tx],
                     [0,        1,      0,      0],
                     [0,        0,      1,      0],
                     [0,        0,      0,      1]])
#----------------------------------------------------------------------------------------------------------#
def mat_trans_z( tz ):
    return mat([[1,     0,     0,      0],
                     [0,     1,     0,      0],
                     [0,     0,     1,      tz],
                     [0,     0,     0,      1]])
def mat_rot_x( ang ):
    #convert to radians
    c = cos(ang * pi / 180)
    s = sin(ang * pi / 180)
    return mat([[1,     0,      0,      0],
                     [0,     c,     -s,     0],
                     [0,     s,     c,      0],
                     [0,     0,     0,      1]])
#----------------------------------------------------------------------------------------------------------#
def mat_rot_z( ang ):
    #convert to radians
    c = cos(ang * pi / 180)
    s = sin(ang * pi / 180)
    return mat([[c,     -s,    0,     0],
                     [s,     c,     0,     0],
                     [0,     0,     1,     0],
                     [0,     0,     0,     1]])
#----------------------------------------------------------------------------------------------------------#
def Aji(A, alpha, D, theta):
    Rz_J = mat_rot_z(theta)
    Tz_J = mat_trans_z(D)
    Tx_I = mat_trans_x(A)
    Rx_I = mat_rot_x(alpha)
    return n.round((Rz_J.dot(Tz_J).dot(Tx_I).dot(Rx_I))*1000)/1000.0
#----------------------------------------------------------------------------------------------------------#
def DH(a,b,c,d,e,f):
    a1 = Aji(0,90,290,180+a) #ok
    a2 = Aji(340,0,0,90+b) #ok
    a3 = Aji(0, 90, 0, 180+c) #ok
    a4 = Aji(0, 90, 302, 180+d) #ok
    a5 = Aji(0, 90, 0, 180+e) #ok
    a6 = Aji(0, 0, 72, 0+f) #ok
    a = a1.dot(a2).dot(a3).dot(a4).dot(a5).dot(a6) #ok
    return a
#----------------------------------------------------------------------------------------------------------#
def vec_to_mat44(v):
    return v.reshape((4,4))
#----------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    data = parse.parse_file("C:\\robot-studio-output.txt")
