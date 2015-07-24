from __future__ import division
#--------------------------#
import sys
#--------------------------#
import numpy as n
#--------------------------#
from helperfunctions_math import *
from denavit_hartenberg import *
#=====================================================#
sys.path.append("../int/misc-tools/")
import parsingtools as  parse
#----------------------------------------------------------------------------------------------------------#
def calc_tool_IRB140(a,b,c,d,e,f):
#    tool0, Ai, DH = DH_params(
#                            -70, 90, 352, 180 + a,'R',
#                            360, 0, 0, 90 + b,'R',
#                            0, 90, 0, 180 + c,'R',
#                            0, 90, 380, 180 + d,'R',
#                            0, 90, 0, 180 + e,'R',
#                            0,0,65,0 + f,'R',
#                            unit='mm')

    tool0, Ai, DH = DH_params(
                            -70, 90, 352, 180,'R',
                            360, 0, 0, 90,'R',
                            0, 90, 0, 180,'R',
                            0, 90, 380, 180,'R',
                            0, 90, 0, 180,'R',
                            0,0,65,0,'R',
                            unit='mm')
    tool0, Ai, DH = calc_tool(DH, a,b,c,d,e,f)
    return tool0, Ai, DH
#----------------------------------------------------------------------------------------------------------#
def calc_tool(DH, *joint_values):
    joint_types = DH['table'][4 :: 5]
    print str(joint_types)
    revolute_index = DH['order'].index('theta')
    prismatic_index = DH['order'].index('D')
    table = list( DH['table'] )

    for i,k in enumerate(joint_types):
        if k == 'R':
            table[i*5+3] += joint_values[i]

    tool0, Ai, DH_table = DH_params(*table, unit='mm')
    return tool0, Ai, DH_table
#----------------------------------------------------------------------------------------------------------#
def calc_tool_IRB140_sub(a,b,c):
    tool0, Ai, DH_table = DH_params(
                            -70, 90, 352, 180 + a, 'R',
                            360, 0, 0, 90 + b, 'R',
                            0, 90, 0, 180 + c, 'R',
                            unit='mm')
    return tool0, Ai, DH_table
#----------------------------------------------------------------------------------------------------------#
def __IK_irb140__orientation(j1, j2, j3, T44):
    #Calculate last angles
    R = T44[0:3,0:3]    
    H3, _, _ = calc_tool_IRB140_sub(j1, j2, j3)
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

def __IK_irb140_position_elbow_up(T44, flipped=False):
    #Geometrical paramaters
    wcp = calc_wcp(T44, 0.065)

    #First angle - j1, used to adjust a point-position
    j1 = atan2(wcp[1],wcp[0])
    if flipped:
        if j1 > 0:
            j1 = j1 - 180
        elif j1 < 0:
            j1 = j1 + 180

    p0 = mat([70e-3, 0, 352e-3])
    p0 = homogenous_rotation_z(j1)[0:3,0:3].dot(p0)

    x0 = norm(wcp[0:2] - p0[0:2])
    h1 = 0.352
    h2 = wcp[2]
    s = abs(h2 - h1)
    x1 = norm(p0 - wcp)
    beta = 380e-3
    alpha = 360e-3
    m = atan(70e-3/352e-3)

    if not flipped:
        ### elbow-up ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = th3 - 90

        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 = th21 + th22
        j2 = 90 - th2
    else:
        ### elbow-up (actually inverse, upside-down) ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(90 - th3)
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 = th21 - th22
        j2 = -(90-th2)

    j4, j5, j6,\
        j41,j51,j61 = __IK_irb140__orientation(j1, j2, j3, T44)
    return (j1, j2, j3, j4, j5, j6), (j1, j2, j3, j41, j51, j61)

def __IK_irb140_position_elbow_down(T44, flipped=False):
    #Geometrical paramaters
    wcp = calc_wcp(T44, 0.065)

    #First angle - j1, used to adjust a point-position
    j1 = atan2(wcp[1],wcp[0])
    if flipped is True:
        ### elbow-down (actually inverse, upside-down) ###
        if j1 > 0:
            j1 = j1 - 180
        elif j1 < 0:
            j1 = j1 + 180

    p0 = mat([70e-3, 0, 352e-3])
    p0 = homogenous_rotation_z(j1)[0:3,0:3].dot(p0)

    x0 = norm(wcp[0:2] - p0[0:2])
    h1 = 0.352
    h2 = wcp[2]
    s = abs(h2 - h1)
    x1 = norm(p0 - wcp)
    beta = 380e-3
    alpha = 360e-3
    m = atan(70e-3/352e-3)
    if not flipped:
        ### elbow-down ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(th3 + 90)
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 =  th21 - th22
        j2 = 90 - th2
    else:
        ### elbow-down (actually inverse, upside-down) ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(th3 + 90)
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        th2 =  th21 + th22
        j2 = -(90 - th2)

    j4, j5, j6,\
        j41,j51,j61 = __IK_irb140__orientation(j1, j2, j3, T44)
    return (j1, j2, j3, j4, j5, j6), (j1, j2, j3, j41, j51, j61)
    
def check_range(x, _min, _max, inclusive=True):

    #swap if needed
    if _max < _min:
        _max, _min = _min, _max

    if inclusive == True:
        return _min <= x <= _max
    else:
        return _min < x < _max

def check_solution(j1,j2,j3,j4,j5,j6, inclusive=True):
    sol  = check_range(j1, -180, 180, inclusive)
    sol &= check_range(j2, -90,  110, inclusive)
    sol &= check_range(j3, -230, 50,  inclusive)
    sol &= check_range(j4, -200, 200, inclusive)
    sol &= check_range(j5, -115, 115, inclusive)
    sol &= check_range(j6, -400, 400, inclusive)
    return sol

def calc_inv_kin_IRB140(T44):
    if type(T44) is list:
        T44 = mat(T44)
    dim = T44.shape
    if len(dim) != 2:
        raise ArithmeticError('Forward-kinematics must be a 4x4 matrix!')
    if dim[0] != dim[1]:
        raise ArithmeticError('Forward-kinematics must be square!')
    if dim[0] != 4:
        raise ArithmeticError('Forward-kinematics must have dimension of 4!')

    sol1, sol11 = __IK_irb140_position_elbow_up(T44)
    sol2, sol21 = __IK_irb140_position_elbow_down(T44)
    sol3, sol31 = __IK_irb140_position_elbow_up(T44, flipped = True)
    sol4, sol41 = __IK_irb140_position_elbow_down(T44, flipped = True)

    #first columnt is first solution and so forth
    return mat(zip(sol1, sol2, sol3, sol4, sol11, sol21, sol31, sol41))

def filter_solutions(solutions, filter_function = check_solution):
    result = []
    for s in solutions.T:
        if filter_function(*s) == True:
            result.append( s )
    return mat(zip(*result))

def calc_valid_inv_kin_IRB140(T44):
    return filter_solutions( calc_inv_kin_IRB140(T44) )

def create_T44(pos, orientation):
    T44 = n.zeros((4,4))
    T44[0:3,0:3] = orientation
    T44[0:3,3] = pos
    T44[3, :] = [0,0,0,1]
    return T44
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

up_to = lambda i: custom_round(matmul(*[debug[x] for x in range(i)]))
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

    A, debug, DH = calc_tool_IRB140(a,b,c,d,e,f)
    calc_tool(DH, a,b,c,d,e,f)

    print "T44 sanity check-norm: " + str(norm(T44 - A))

    p_end = T44[0:3,3]
    wcp = calc_wcp(T44, 0.065)
    sol = mat( calc_inv_kin_IRB140(T44) )
    s0 = mat([a,b,c,d,e,f])
    all_norms = 0
    for i in xrange(0, 8):
        s = sol[:,i]
        gamma0,gamma1,gamma2,gamma3,gamma4,gamma5 = s
        A, debug, DH_table = calc_tool_IRB140(gamma0, gamma1, gamma2,
                                              gamma3, gamma4, gamma5)
        p0 = debug[0][:,3]
        p1 = matmul(debug[0],debug[1])[:,3]
        p2 = matmul(debug[0],debug[1],debug[2])[:,3]
        p3 = matmul(debug[0],debug[1],debug[2], debug[3])[:,3]
        p4 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4])[:,3]
        p5 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4], debug[5])[:,3]
        print "\n[ Solution %s ]" % str(i)
        print "FK-norm: " + str( norm(A - T44) )
        all_norms = all_norms + norm(A - T44)
        print "angle-norm: %0.2f" % norm(s - s0)

        #Plotting
        from pylab import plot, show, legend
        M = mat(zip([0,0,0],p0,p1,p2,p3,p4,p5)).T
        if (i % 4) == 0:
            col = 'b-'
            lw = 3
        if (i % 4) == 1:
            col = 'r-'
            lw = 3
        if (i % 4) == 2:
            col = 'b-.'
            lw = 2
        if (i % 4) == 3:
            col = 'r-.'
            lw = 2
        plot(M[:,0], M[:,2], col, linewidth = lw)
        legend(['elbow-up', 'elbow-down', 'elbow-up-flipped', 'elbow-down-flipped',
                'elbow-up-2', 'elbow-down-2', 'elbow-up-flipped-2', 'elbow-down-flipped-2'])

    print "FK-norm-summary: " + str( all_norms )
    plot([-1,-1,1,1],[0,1,1,0],'w')
    plot(wcp[0], wcp[2], 'ro')
    plot(p_end[0], p_end[2], 'ko')
    show()
