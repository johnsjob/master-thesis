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

matmul = matmul_series
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
#----------------------------------------------------------------------------------------------------------#
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
    return matmul_series(Rz_J, Tz_J, Tx_I, Rx_I)
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
    return matmul_series(*matrices), mat(matrices)
#----------------------------------------------------------------------------------------------------------#
def calc_tool_IRB120(a,b,c,d,e,f):
    flange = DH_params(
                    0,      90, 0.290,  180+a,
                    0.270,   0, 0,      90+b,
                   -0.070,  90, 0,      180+c,
                    0,      90, 0.302,  180+d,
                    0,      90, 0,      180+e,
                    0,      0, 0.072,   0+f
                    )
    return flange
#----------------------------------------------------------------------------------------------------------#
def custom_round(v, prec = 1e-4):
    coef = 1 / prec
    return n.round(v * coef) / coef
#----------------------------------------------------------------------------------------------------------#
def clear():
    for i in xrange(0,100):
        print ''
#----------------------------------------------------------------------------------------------------------#
rad = lambda x: x * pi / 180.0
deg = lambda x: x * 180.0 / pi
up_to = lambda i: matmul_series(*[debug[x] for x in range(i)])
cos_sats = lambda a,b,th: a**2 + b**2 - 2*a*b*cos(rad(th)); #ok
ang_sats = lambda c,a,b: deg(acos((c**2 - a**2 - b**2)/(-2*a*b))); #ok
#----------------------------------------------------------------------------------------------------------#
from numpy.linalg import norm
from numpy import arctan2 as atan2, arccos as acos, sqrt, arctan as atan
#----------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    """
    GENERAL COMMENTS:
    ################
        data['Joint_(i)_T'] = cumulative multiplication of up to i-1 of joint T-matrices

        'debug' contains i:th matrices of python implementation of DH-algorithm

        lambda function up_to(i) performs cumulative multiplication up to i:th joint
    """
    index = 0

    data = parse.parse_file("C:\\robot-studio-output.txt")
    num_joint_conf = len(data['Joint_1_T'])

    for key in data:
        if 'T' in key:
            data[key] = custom_round(mat(data[key]).reshape(num_joint_conf,4,4))
            
            
    params = zip(data['A'], data['alpha'], data['D'], data['theta'])
    T44 = mat(data['T44'][index]).reshape((4,4))

    print "\nNumber of configurations: " + str(len(data['Joint_1'])) + "\n"
    a,b,c,d,e,f = data['Joint_1'][index], data['Joint_2'][index], data['Joint_3'][index], data['Joint_4'][index], data['Joint_5'][index], data['Joint_6'][index],

    A, debug = calc_tool_IRB120(a,b,c,d,e,f)
    debug = custom_round(debug)

    print "T44 error compared to robot-studio: \n" + str(n.abs(T44 - A))
    print "Error norm: " + str(n.linalg.norm( n.abs(T44 - A) ))

    #INVERSE KINEMATICS STARTS HERE
    wcp = (A[:,3]-A[:,2]*0.072)[0:3]
    print "wcp-norm: "+str(norm(data['Joint_6_T'][index][0:3,3] - wcp))
    x0 = norm((wcp[0],wcp[1]))
    
    h1 = norm(mat([0,0,290e-3]))
    x0p = norm(mat([0,0,h1]) - mat([wcp[0],wcp[1],0]))
    h2 = wcp[2]
    x1 = norm(mat([0,0,h1]) - wcp[0:3])

    gamma0 = (atan2(wcp[1],wcp[0])*180/pi)
    print 'a-norm: ' + str( a - gamma0 )

    beta = sqrt(70e-3**2 + 302e-3**2)
    alpha = 270e-3
    th1 = ang_sats(x0, h1, x0p)
    
    th2 = ang_sats(h2, x0p, x1)
    
    th3 = ang_sats(beta, alpha, x1)
    
    gamma1 = 180-(th1+th2+th3)
    print 'b-norm: ' + str( b - gamma1 )

    th41 = deg(atan(0.070 / 0.302))
    s = h2 - h1
    th4 = acos((x1**2 - alpha**2 - beta**2)/(-2*alpha*beta)) * 180/pi
    th4 = th4 - th41
    #th4 = acos((x1**2 - (270e-3 + 70e-3)**2 - 302e-3**2)/(-2*(270e-3+70e-3)*302e-3)) * 180/pi
    print 'c-norm: ' + str( th4-c )
    

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
