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
    """
    Calculats transform from frame J = I-1 to I
    in order AI = Rzj( theta )Tzj( D )Txi( A )Rxi( alpha )
    and parameters are given in order A, alpha, D, theta.
    """
    Rz_J = mat_rot_z(theta)
    Tz_J = mat_trans_z(D)
    Tx_I = mat_trans_x(A)
    Rx_I = mat_rot_x(alpha)
    return matmul_series(Rz_J, Tz_J, Tx_I, Rx_I)
#----------------------------------------------------------------------------------------------------------#
def DH_params(*params, **kwargs):
    """
    Performs the denavit-hartenberg algorithm
    and calculates T44 = A0 * A1 * A2 * ... * An multiplication.

    The parameters are entered in the order A, alpha, D, theta repeatedly.
    """

    # handle which unit is being used
    if not kwargs.has_key('unit'):
        kwargs['unit'] = 'm'
    unit = kwargs['unit']

    nbr_of_sections = int(len(params) / 4)
    if len(params) == 1 and type(params[0]) in [list, tuple]:
        raise ArithmeticError("Function does not use lists or tuples, please unpack using *.")
    elif not (len(params) % 4 == 0):
        raise ArithmeticError("Invalid number of Denavit-Hartenberg parameters.")

    matrices = []
    for k in xrange(0, nbr_of_sections):
        A, alpha, D, theta = params[4*k:4*k+4]
        if unit == 'mm':
            A = A * 1e-3
            D = D * 1e-3
        elif unit == 'm':
            pass
        else:
            raise ArithmeticError("Unknown unit of length, only meters(\'m\') or millimeters meters(\'mm\') allowed.")
            
        matrices.append( transform_to_next(A, alpha, D, theta) )
    return matmul_series(*matrices), mat(matrices)
#----------------------------------------------------------------------------------------------------------#
def calc_tool_IRB140(a,b,c,d,e,f):
    tool0, Ai = DH_params(
                            -70, 90, 352, 180 + a,
                            360, 0, 0, 90 + b,
                            0, 90, 0, 180+c,
                            0, 90, 380, 180+d,
                            0, 90, 0, 180+e,
                            0,0,65,f,
                            unit='mm')
    return tool0, Ai
#----------------------------------------------------------------------------------------------------------#
def calc_tool_IRB120(a,b,c,d,e,f):
    tool0, Ai = DH_params(
                    0,      90, 0.290,  180+a,
                    0.270,   0, 0,      90+b,
                   -0.070,  90, 0,      180+c,
                    0,      90, 0.302,  180+d,
                    0,      90, 0,      180+e,
                    0,      0, 0.072,   0+f,
                    unit = 'm')
    return tool0, Ai
#----------------------------------------------------------------------------------------------------------#
def calc_tool_IRB120_sub(a,b,c):
    flange = DH_params(
                    0,      90, 0.290,  180+a,
                    0.270,   0, 0,      90+b,
                   -0.070,  90, 0,      180+c,
                    unit = 'm')
    return flange
#----------------------------------------------------------------------------------------------------------#
def __IK_irb120__orientation(j1, j2, j3, T44):
    #Calculate last angles
    R = T44[0:3,0:3]    
    H3,_ = calc_tool_IRB120_sub(j1, j2, j3)
    R3 = H3[0:3, 0:3]
    R36 = R3.T.dot(R)
    X = R36[:,0]
    Y = R36[:,1]
    Z = R36[:,2]
    # for order of parameters check numpy.info(numpy.arctan2)
    j4 = atan2(Z[1],Z[0])
    j5 = atan2(norm(Z[0:2]), Z[2])
    j6 = atan2(X[2], Y[2]) + 90

    R36 = R36.T
    X = R36[:,0]
    Y = R36[:,1]
    Z = R36[:,2]
    # for order of parameters check numpy.info(numpy.arctan2)
    j41 = -(atan2(X[2], Y[2]) + 90)
    j51 = -atan2(norm(Z[0:2]), Z[2])
    j61 = -(atan2(Z[1],Z[0]))
    return j4, j5, j6, j41, j51, j61

def __IK_irb120_position_elbow_up(T44, flipped=False):
    #Geometrical paramaters
    wcp = calc_wcp(T44)
    x0 = norm((wcp[0],wcp[1]))
    h1 = norm(mat([0,0,290e-3]))
    x0p = norm(mat([0,0,h1]) - mat([wcp[0],wcp[1],0]))
    h2 = wcp[2]
    s = abs(h2 - h1)
    x1 = norm(mat([0,0,h1]) - wcp[0:3])
    beta = sqrt(0.070**2 + 0.302**2)
    alpha = 270e-3
    m = atan(0.070 / 0.302)
    #First angle - j1
    j1 = atan2(wcp[1],wcp[0])
    if not flipped:
        ### elbow-up ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = th3 + m - 90

        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 = th21 + th22
        j2 = 90 - th2
    else:
        ### elbow-up (actually inverse, upside-down) ###
        if j1 > 0:
            j1 = j1 - 180
        elif j1 < 0:
            j1 = j1 + 180
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(90 - (th3 + m))
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 = th21 - th22
        j2 = -(90-th2)

    j4, j5, j6,\
        j41,j51,j61 = __IK_irb120__orientation(j1, j2, j3, T44)
    return (j1, j2, j3, j4, j5, j6), (j1, j2, j3, j41, j51, j61)

def __IK_irb120_position_elbow_down(T44, flipped=False):
    #Geometrical paramaters
    wcp = calc_wcp(T44)
    x0 = norm((wcp[0],wcp[1]))
    h1 = norm(mat([0,0,290e-3]))
    x0p = norm(mat([0,0,h1]) - mat([wcp[0],wcp[1],0]))
    h2 = wcp[2]
    s = abs(h2 - h1)
    x1 = norm(mat([0,0,h1]) - wcp[0:3])
    beta = sqrt(0.070**2 + 0.302**2)
    alpha = 270e-3
    m = atan(0.070 / 0.302)
    #First angle - j1
    j1 = atan2(wcp[1],wcp[0])
    if not flipped:
        ### elbow-down ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(th3 - m + 90)
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 =  th21 - th22
        j2 = 90 - th2
    else:
        ### elbow-down (actually inverse, upside-down) ###
        if j1 > 0:
            j1 = j1 - 180
        elif j1 < 0:
            j1 = j1 + 180
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(th3 - m + 90)
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 =  th21 + th22
        j2 = -(90 - th2)

    j4, j5, j6,\
        j41,j51,j61 = __IK_irb120__orientation(j1, j2, j3, T44)
    return (j1, j2, j3, j4, j5, j6), (j1, j2, j3, j41, j51, j61)
    
def calc_wcp(T44):
    return (T44[:,3] - T44[:,2]*0.072)[0:3]

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

    sol1, sol11 = __IK_irb120_position_elbow_up(T44)
    sol2, sol21 = __IK_irb120_position_elbow_down(T44)
    sol3, sol31 = __IK_irb120_position_elbow_up(T44, flipped = True)
    sol4, sol41 = __IK_irb120_position_elbow_down(T44, flipped = True)

    #first columnt is first solution and so forth
    return zip(sol1, sol2, sol3, sol4, sol11, sol21, sol31, sol41)
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

    A, debug = calc_tool_IRB140(a,b,c,d,e,f)

    print "T44 sanity check-norm: " + str(norm(T44 - A))

##    #INVERSE KINEMATICS STARTS HERE
##    p_end = T44[0:3,3]
##    wcp = calc_wcp(T44)
##    
##    sol = mat( calc_inv_kin_IRB120(T44) )
##    s0 = mat([a,b,c,d,e,f])
##    for i in xrange(0, 8):
##        s = sol[:,i]
##        gamma0,gamma1,gamma2,gamma3,gamma4,gamma5 = s
##        A, debug = calc_tool_IRB120(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5)
##        p0 = debug[0][:,3]
##        p1 = matmul(debug[0],debug[1])[:,3]
##        p2 = matmul(debug[0],debug[1],debug[2])[:,3]
##        p3 = matmul(debug[0],debug[1],debug[2], debug[3])[:,3]
##        p4 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4])[:,3]
##        p5 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4], debug[5])[:,3]
##        print "\n[ Solution %s ]" % str(i)
##        print "FK-norm: " + str( norm(A - T44) )
##        print "angle-norm: %0.2f" % norm(s - s0)
##
##        #Plotting
##        from pylab import plot, show, legend
##        M = mat(zip([0,0,0],p0,p1,p2,p3,p4,p5)).T
##        if (i % 4) == 0:
##            col = 'b-'
##            lw = 3
##        if (i % 4) == 1:
##            col = 'r-'
##            lw = 3
##        if (i % 4) == 2:
##            col = 'b-.'
##            lw = 2
##        if (i % 4) == 3:
##            col = 'r-.'
##            lw = 2
##        plot(M[:,0], M[:,2], col, linewidth = lw)
##        legend(['elbow-up', 'elbow-down', 'elbow-up-flipped', 'elbow-down-flipped',
##                'elbow-up-2', 'elbow-down-2', 'elbow-up-flipped-2', 'elbow-down-flipped-2'])
##    plot([-1,-1,1,1],[0,1,1,0],'w')
##    plot(wcp[0], wcp[2], 'ro')
##    plot(p_end[0], p_end[2], 'ko')
##    show()

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
