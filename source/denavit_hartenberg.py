from __future__ import division
#--------------------------#
import sys
#--------------------------#
import numpy as n
from numpy import cos, sin, pi
from numpy import array as mat
#=====================================================#
from helperfunctions import matmul_series
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
def transform_to_next(A, alpha, D, theta):
    Rz_J = mat_rot_z(theta)
    Tz_J = mat_trans_z(D)
    Tx_I = mat_trans_x(A)
    Rx_I = mat_rot_x(alpha)
    return matmul_series([Rz_J, Tz_J, Tx_I, Rx_I])
#----------------------------------------------------------------------------------------------------------#
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
    return matmul_series(matrices)
#----------------------------------------------------------------------------------------------------------#
def calc_tool_IRB120(a,b,c,d,e,f):
    flange = DH_params(
                        0, 90,0.290,180+a,
                       0.270,0,0,90+b,
                      -0.070, 90, 0, 180+c,
                      0, 90, 0.302, 180+d,
                      0, 90, 0, 180+e,
                      0, 0, 0.072, 0+f
                        )
    return flange
#----------------------------------------------------------------------------------------------------------#
def custom_round(v, prec = 1e-4):
    coef = 1 / prec
    return n.round(v * coef) / coef
#----------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    index = -1
    data = parse.parse_file("C:\\robot-studio-output.txt")
    params = zip(data['A'], data['alpha'], data['D'], data['theta'])
    T44 = mat(data['T44'][index]).reshape((4,4))
    a,b,c,d,e,f = data['Joint_1'][index], data['Joint_2'][index], data['Joint_3'][index], data['Joint_4'][index], data['Joint_5'][index], data['Joint_6'][index],

    A= calc_tool_IRB120(a,b,c,d,e,f)

    print n.abs(T44 - A)
    print n.linalg.norm( n.abs(T44 - A) )

###############################################################
##
##    Robot Studio notes for IRB120:
##    DH-parameters (NOT modified Denavit-Hartenbergaccording to the SDK):
##---------------------------------------------------------------------------------------
##    a1 = transform_to_next(0,0,0,0) #ok
##    a2 = transform_to_next(0,90,0,90) #ok
##    a3 = transform_to_next(0.270, 0, 0, 0) #ok
##    a4 = transform_to_next(0.070, -90, 0.302, 0) #ok
##    a5 = transform_to_next(0, 90, 0, -180) #ok
##    a6 = transform_to_next(0, 90, 0, 0) #ok
##    a = a1.dot(a2).dot(a3).dot(a4).dot(a5).dot(a6) #ok
##
##
##
##
##
##
##
##
##
##
##
