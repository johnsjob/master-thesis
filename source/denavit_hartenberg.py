from __future__ import division
#--------------------------#
import sys
#--------------------------#
import numpy as n
from numpy import cos, sin, pi
from numpy import array as mat

from sympy import init_printing, pprint
#from Ipython.display import display

init_printing(use_latex='mathjax')
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
def calc_tool_IRB120_sub(a,b,c):
    flange = DH_params(
                    0,      90, 0.290,  180+a,
                    0.270,   0, 0,      90+b,
                   -0.070,  90, 0,      180+c
                    )
    return flange
#----------------------------------------------------------------------------------------------------------#
def calc_inv_kin_IRB120(T44):
    if type(T44) is list:
        T44 = mat(T44)
    dim = T44.shape
    if len(dim) != 2:
        raise ArithmeticError('Forward-kinematics must be a 4x4 matrix!')
    if dim[0] != dim[1]:
        raise ArithmeticError('Forward-kinematics must be square!')
    if dim[0] != 4:
        raise ArithmeticError('Forward-kinematics must have dimension of 4!')
#----------------------------------------------------------------------------------------------------------#
def custom_round(v, prec = 1e-8):
    coef = 1 / prec
    return n.round(v * coef) / coef
#----------------------------------------------------------------------------------------------------------#
def clear():
    for i in xrange(0,100):
        print ''
#----------------------------------------------------------------------------------------------------------#
from numpy.linalg import norm
from numpy import arctan2 as atan2, arccos as acos, arcsin as asin, sqrt, arctan as atan
#----------------------------------------------------------------------------------------------------------#
rad = lambda x: x * pi / 180.0
deg = lambda x: x * 180.0 / pi
cos2 = lambda x: n.cos(rad(x))
sin2 = lambda x: n.sin(rad(x))

cos_sats = lambda a,b,th: a**2 + b**2 - 2*a*b*cos(rad(th)); #ok
ang_sats = lambda c,a,b: deg(acos((c**2 - a**2 - b**2)/(-2*a*b))); #ok
ang_sats2 = lambda c,a,b: deg(acos((c**2 - a**2 - b**2)/(2*a*b))); #ok
round = lambda x: custom_round(x)
atan = lambda x: deg(n.arctan(x))
atan2 = lambda y,x: deg(n.arctan2(y,x))

up_to = lambda i: custom_round(matmul_series(*[debug[x] for x in range(i)]))
#----------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    """
    GENERAL COMMENTS:
    ################
        data['Joint_(i)_T'] = cumulative multiplication of up to i-1 of joint T-matrices

        'debug' contains i:th matrices of python implementation of DH-algorithm

        lambda function up_to(i) performs cumulative multiplication up to i:th joint

        wcp = data['Joint_6_T'][index] = up_to(5)
    """
    clear()
    index = -1

    data = parse.parse_file("C:\\robot-studio-output.txt")
    num_joint_conf = len(data['Joint_1_T'])
            
            
    params = zip(data['A'], data['alpha'], data['D'], data['theta'])
    T44 = mat(data['T44'][index]).reshape((4,4))

    for key in data:
        if 'T' in key:
            data[key] = custom_round(mat(data[key]).reshape(num_joint_conf,4,4))

    print "\nNumber of configurations: " + str(len(data['Joint_1'])) + "\n"
    a,b,c,d,e,f = data['Joint_1'][index], data['Joint_2'][index], data['Joint_3'][index], data['Joint_4'][index], data['Joint_5'][index], data['Joint_6'][index],

    A, debug = calc_tool_IRB120(a,b,c,d,e,f)

    print "T44 sanity check-norm: " + str(norm(T44 - A))

    #INVERSE KINEMATICS STARTS HERE
    p_end = T44[0:3,3]
    wcp = (T44[:,3] - T44[:,2]*0.072)[0:3]
        
    print "wcp-norm: "+str(norm(data['Joint_6_T'][index][0:3,3] - wcp))
    x0 = norm((wcp[0],wcp[1]))
    
    h1 = norm(mat([0,0,290e-3]))
    x0p = norm(mat([0,0,h1]) - mat([wcp[0],wcp[1],0]))
    h2 = wcp[2]
    s = abs(h2 - h1)
    x1 = norm(mat([0,0,h1]) - wcp[0:3])

    #First angle - j1
    gamma0 = atan2(wcp[1],wcp[0])
    print 'a-norm: ' + str( norm(a - gamma0 ))

    beta = sqrt(0.070**2 + 0.302**2)
    alpha = 270e-3
    m = atan(0.070 / 0.302)


##    ### elbow-up ###
##    # Third angle - j3
##    th3 = ang_sats2(x1, alpha, beta)
##    gamma2 = th3 + m - 90
##    print 'c-norm: ' + str(norm( gamma2-c ))
##    # Second angle - j2
##    th21 = atan2(s, x0)
##    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
##    th2 = th21 + th22
##    gamma1 = 90 - th2
##    print 'b-norm: ' + str(norm( b - gamma1 ))


##    ### elbow-up (actually inverse, upside-down) ###
##    # Third angle - j3
##    if gamma0 > 0:
##        gamma0 = gamma0 - 180
##    elif gamma0 < 0:
##        gamma0 = gamma0 + 180
##    th3 = ang_sats2(x1, alpha, beta)
##    gamma2 = -(90 - th3 - m)
##    print 'c-norm: ' + str(norm( gamma2-c ))
##    # Second angle - j2
##    th21 = atan2(s, x0)
##    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
##    th2 = th21 - th22
##    gamma1 = -(90-th2)
##    print 'b-norm: ' + str(norm( b - gamma1 ))

##    ### elbow-down ###
##    # Third angle - j3
##    th3 = ang_sats2(x1, alpha, beta)
##    gamma2 = -(th3 - m + 90)
##    print 'c-norm: ' + str(norm( gamma2-c ))
##    # Second angle - j2
##    th21 = atan2(s, x0)
##    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
##    th2 =  th21 - th22
##    gamma1 = 90 - th2
##    print 'b-norm: ' + str(norm( b - gamma1 ))

##    ### elbow-down (actually inverse, upside-down) ###
##    # Third angle - j3
##    if gamma0 > 0:
##        gamma0 = gamma0 - 180
##    elif gamma0 < 0:
##        gamma0 = gamma0 + 180
##    th3 = ang_sats2(x1, alpha, beta)
##    gamma2 = -(th3 - m + 90)
##    # Second angle - j2
##    print 'c-norm: ' + str(norm( gamma2-c ))
##    th21 = atan2(s, x0)
##    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
##    th2 =  th21 + th22
##    gamma1 = -(90 - th2)
##    print 'b-norm: ' + str(norm( b - gamma1 ))

    #Calculate last angles
    R = T44[0:3,0:3]    
    H3,_ = calc_tool_IRB120_sub(gamma0, gamma1, gamma2)
    R3 = H3[0:3, 0:3]
    R36 = R3.T.dot(R)
    X = R36[:,0]
    Y = R36[:,1]
    Z = R36[:,2]
    # for order of parameters check numpy.info(numpy.arctan2)
    gamma3 = atan2(Z[1],Z[0])
    print 'd-norm: ' + str(norm( gamma3-d ))

    # for order of parameters check numpy.info(numpy.arctan2)
    gamma4 = atan2(norm(Z[0:2]), Z[2])
    print 'e-norm: ' + str(norm( gamma4-e ))

    gamma5 = atan2(X[2], Y[2]) + 90
    print 'f-norm: ' + str(norm( gamma5-f ))

    A, debug = calc_tool_IRB120(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5)
    p0 = debug[0][:,3]
    p1 = matmul(debug[0],debug[1])[:,3]
    p2 = matmul(debug[0],debug[1],debug[2])[:,3]
    p3 = matmul(debug[0],debug[1],debug[2], debug[3])[:,3]
    p4 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4])[:,3]
    p5 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4], debug[5])[:,3]

    print "FK-norm: " + str(norm(A - T44))

    #Plotting
    from pylab import plot, show
    M = mat(zip([0,0,0],p0,p1,p2,p3,p4,p5)).T
    plot(wcp[0], wcp[2], 'ro')
    plot(p_end[0], p_end[2], 'ko')
    plot(M[:,0],M[:,2])
    plot([-1,-1,1,1],[0,1,1,0],'w')
    show()

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
