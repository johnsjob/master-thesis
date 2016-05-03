from __future__ import division
#--------------------------#
import sys
#--------------------------#
import numpy as n
#--------------------------#
import unittest
#--------------------------#
import numpy
from numpy import pi, arctan2 as atan2, arccos as acos,\
     arcsin as asin, sqrt, arctan as atan
from numpy.linalg import norm
from helperfunctions_math import rand_range, mat, homogenous_rotation_z
from denavit_hartenberg import inverse_kinematics_spherical_wrist as inv_wrist,\
     forward_kinematics, calc_wcp, calc_j1, pack_elbow_and_wrists
#=====================================================#
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
DH_TABLE = {  'table':[-70, 90, 352, 180, 'R',
                       360,  0,   0,  90, 'R',
                         0, 90,   0, 180, 'R',
                         0, 90, 380, 180, 'R',
                         0, 90,   0, 180, 'R',
                         0,  0,  65,   0, 'R'],
             'unit': 'mm',
             'order': ['A','alpha','D','theta'],
             'convention': 'standard'
            }
#----------------------------------------------------------------------------------------------------------#
def backward_facing_elbow_down(dh_table, T44):
    '''
    Calculates backward-facing elbow-down.

    Note: This is in reality an elbow-down solution,
          the name merely implies it is "flipped" up-side-down
          due to the base have turned 180 degrees.
    '''
    #Geometrical paramaters
    wcp = calc_wcp(T44, 0.065)

    #First angle - j1, used to adjust a point-position
    j1 = calc_j1(wcp, flipped=True)

    p0 = mat([70e-3, 0, 352e-3])
    p0 = homogenous_rotation_z(j1)[0:3,0:3].dot(p0)

    x0 = norm(wcp[0:2] - p0[0:2])
    h1 = p0[2]
    h2 = wcp[2]
    s = h2 - h1
    x1 = norm(p0 - wcp)
    beta = 380e-3
    alpha = 360e-3

    th3 = ang_sats2(x1, alpha, beta)
    j3 = -90 + th3

    th21 = atan2(s, x0)
    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
    j2 = -90 + th21 - th22
            
    #packing the solutions in a dynamic way
##    j4, j5, j6,\
##    j41,j51,j61, \
##    j42,j52,j62 = inverse_kinematics_spherical_wrist(dh_table, j1, j2, j3, T44)
##
##    return (j1, j2, j3, j4, j5, j6),\
##           (j1, j2, j3, j41, j51, j61), \
##           (j1, j2, j3, j42, j52, j62)
    result = pack_elbow_and_wrists(dh_table, j1, j2, j3, T44)
    return result
#----------------------------------------------------------------------------------------------------------#
def backward_facing_elbow_up(dh_table, T44):
    '''
    Calculates backward-facing elbow-up.

    Note: This is in reality an elbow-up solution,
          the name merely implies it is "flipped" up-side-down
          due to the base have turned 180 degrees.
    '''
    #Geometrical paramaters
    wcp = calc_wcp(T44, 0.065)

    #First angle - j1, used to adjust a point-position
    j1 = calc_j1(wcp, flipped=True)
    
    p0 = mat([70e-3, 0, 352e-3])
    p0 = homogenous_rotation_z(j1)[0:3,0:3].dot(p0)

    x0 = norm(wcp[0:2] - p0[0:2])
    h1 = p0[2]
    h2 = wcp[2]
    s = h2 - h1
    x1 = norm(p0 - wcp)
    beta = 380e-3
    alpha = 360e-3

    th3 = ang_sats2(x1, alpha, beta)
    j3 = -90 - th3

    th21 = atan2(s, x0)
    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
    j2 = -90 + (th21 + th22)

##    j4, j5, j6,\
##    j41,j51,j61, \
##    j42,j52,j62 = inverse_kinematics_spherical_wrist(dh_table, j1, j2, j3, T44)
##
##    return (j1, j2, j3, j4, j5, j6),\
##           (j1, j2, j3, j41, j51, j61), \
##           (j1, j2, j3, j42, j52, j62)
    result = pack_elbow_and_wrists(dh_table, j1, j2, j3, T44)
    return result
#----------------------------------------------------------------------------------------------------------#
def elbow_up(dh_table, T44):
    '''
    Calculates forward-facing elbow-up.
    '''
    #Geometrical paramaters
    wcp = calc_wcp(T44, 0.065)

    #First angle - j1, used to adjust a point-position
    j1 = calc_j1(wcp, flipped=False)

    p0 = mat([70e-3, 0, 352e-3])
    p0 = homogenous_rotation_z(j1)[0:3,0:3].dot(p0)

    x0 = norm(wcp[0:2] - p0[0:2])
    h1 = p0[2]
    h2 = wcp[2]
    s = h2 - h1
    x1 = norm(p0 - wcp)
    beta = 380e-3
    alpha = 360e-3

    th3 = ang_sats2(x1, alpha, beta)
    j3 = th3 - 90

    th21 = atan2(s, x0)
    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
    j2 = 90 - (th21 + th22)
    if norm(wcp[:2])-norm(p0[:2]) < 0:
        j2 = -90 + (th21 - th22)
##        j3 = -90-th3
            
##    j4, j5, j6,\
##    j41,j51,j61, \
##    j42,j52,j62 = inverse_kinematics_spherical_wrist(dh_table, j1, j2, j3, T44)
####    import pdb; pdb.set_trace()
##    return (j1, j2, j3, j4, j5, j6),\
##           (j1, j2, j3, j41, j51, j61), \
##           (j1, j2, j3, j42, j52, j62)
    result = pack_elbow_and_wrists(dh_table, j1, j2, j3, T44)
    return result
#----------------------------------------------------------------------------------------------------------#
def elbow_down(dh_table, T44):
    '''
    Calculates forward-facing elbow-down.
    '''
    #Geometrical paramaters
    wcp = calc_wcp(T44, 0.065)

    #First angle - j1, used to adjust a point-position
    j1 = calc_j1(wcp, flipped=False)

    p0 = mat([70e-3, 0, 352e-3])
    p0 = homogenous_rotation_z(j1)[0:3,0:3].dot(p0)

    x0 = norm(wcp[0:2] - p0[0:2])
    h1 = p0[2]
    h2 = wcp[2]
    s = h2 - h1
    x1 = norm(p0 - wcp)
    beta = 380e-3
    alpha = 360e-3

    th3 = ang_sats2(x1, alpha, beta)
    j3 = -th3 - 90

    th21 = atan2(s, x0)
    th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
    j2 = 90 - (th21 - th22)
    if norm(wcp[:2])-norm(p0[:2]) < 0:
        j2 = -90 + (th21+th22)
##        j3 = -90-th3
    
##    j4, j5, j6,\
##    j41,j51,j61, \
##    j42,j52,j62 = inverse_kinematics_spherical_wrist(dh_table, j1, j2, j3, T44)
##
##    return (j1, j2, j3, j4, j5, j6),\
##           (j1, j2, j3, j41, j51, j61), \
##           (j1, j2, j3, j42, j52, j62)
    result = pack_elbow_and_wrists(dh_table, j1, j2, j3, T44)
    return result
#----------------------------------------------------------------------------------------------------------#
def inverse_kinematics_elbow_up(dh_table, T44, flipped = False):
    '''
    Wrapper function for forward-facing elbow-up,
    and backward-facing elbow-dup
    '''
    if not flipped:
        return elbow_up(dh_table, T44)
    else:
        return backward_facing_elbow_up(dh_table, T44)

def inverse_kinematics_elbow_down(dh_table, T44, flipped = False):
    '''
    Wrapper function for forward-facing elbow-down,
    and backward-facing elbow-down
    '''
    if not flipped:
        return elbow_down(dh_table, T44)
    else:
        return backward_facing_elbow_down(dh_table, T44)
#----------------------------------------------------------------------------------------------------------#
# INVERSE KINEMATICS - WRAPPERS
#----------------------------------------------------------------------------------------------------------#
def inverse_kinematics_irb140(dh_table, T44):
    '''
    Wrapper function that calculates the 20 analytical solutions
    to the IRB140 robot.
    '''
    if type(T44) is list:
        T44 = mat(T44)
    dim = T44.shape
    if len(dim) != 2:
        raise ArithmeticError('Forward-kinematics must be a 4x4 matrix!')
    if dim[0] != dim[1]:
        raise ArithmeticError('Forward-kinematics must be square!')
    if dim[0] != 4:
        raise ArithmeticError('Forward-kinematics must have dimension of 4!')

    # x5 for each elb_x
    sol_elbup      = inverse_kinematics_elbow_up(dh_table, T44)
    sol_elbdown    = inverse_kinematics_elbow_down(dh_table, T44)
    sol_elbup_backward   = inverse_kinematics_elbow_up(dh_table, T44, flipped = True)
    sol_elbdown_backward = inverse_kinematics_elbow_down(dh_table, T44, flipped = True)

    #first columnt is first solution and so forth
##    ret = mat(zip(sol_elbup1, sol_elbdown1, sol_elbup1_fl, sol_elbdown1_fl,
##                  sol_elbup2, sol_elbdown2, sol_elbup2_fl, sol_elbdown2_fl,
##                  sol_elbup3, sol_elbdown3, sol_elbup3_fl, sol_elbdown3_fl))

    #first columnt is first solution and so forth
    ret = mat( zip( sol_elbup, sol_elbdown, sol_elbup_backward, sol_elbdown_backward ) )
    k,m,n = ret.shape
    ret = ret.reshape(k*m,n)
    return ret.T

def calc_valid_raw_invkin_irb140(T44):
    return __inverse_kinematics_pose(T44, filtering=True, raw_solutions=True)

def calc_valid_invkin_irb140(T44):
    return __inverse_kinematics_pose(T44, filtering=True, raw_solutions=False)

def calc_invkin_irb140(T44, filtering=False, raw_solutions=False):
    return __inverse_kinematics_pose(T44, filtering, raw_solutions)
#----------------------------------------------------------------------------------------------------------#
# INVERSE KINEMATICS - SOLUTION HANDLING
#----------------------------------------------------------------------------------------------------------#
def check_solution(j1,j2,j3,j4,j5,j6, inclusive=True):
    sol  = check_range(j1, -180, 180, inclusive)
    sol &= check_range(j2, -90,  110, inclusive)
    sol &= check_range(j3, -230, 50,  inclusive)
    sol &= check_range(j4, -200, 200, inclusive)
    sol &= check_range(j5, -115, 115, inclusive)
    sol &= check_range(j6, -400, 400, inclusive)
    return sol

def filter_solutions(solutions, filter_function = check_solution):
    result = []
    for s in solutions.T:
        if filter_function(*s) == True:
            result.append( s )
    # returns non-flipped, flipped
    return mat(zip(*result))

def merge_solutions(*args):
    result = []
    for m in args:
        result += zip(*m)
    return mat(zip(*result))

def __modulo_solutions(solution_matrix, index, modulo=360.0):
    for s in solution_matrix.T:
        result = s.copy()
        value = result[index]
        result[index] = value + modulo
        yield result

def generate_modulo_solutions(solution_matrix, index, modulo=360.0):
    return mat(zip(*__modulo_solutions(solution_matrix, index, modulo)))
#----------------------------------------------------------------------------------------------------------#
# INVERSE KINEMATICS - SOLUTION CURVE
#----------------------------------------------------------------------------------------------------------#
def inverse_kinematics_curve(trans_frames):
    # perform inverse kinematics over a curve and collect all solutions
    all_solutions = []
    for point_frame in trans_frames:
        all_solutions.append(_inverse_kinematics_pose(point_frame))
    return mat(all_solutions)

def __inverse_kinematics_pose(T44, filtering=False, raw_solutions=False):
    # perform inverse kinematics on a single frame
    angle_solutions = inverse_kinematics_irb140(DH_TABLE, T44)
    if not raw_solutions:
        extra = [angle_solutions]
        for index in xrange(6):
            extra.append( generate_modulo_solutions(angle_solutions, index, 360.0))
            extra.append( generate_modulo_solutions(angle_solutions, index, -360.0))
            pass
        angle_solutions = merge_solutions(*extra)
    if filtering:
        angle_solutions = filter_solutions(angle_solutions)
    return mat(angle_solutions.T)
#----------------------------------------------------------------------------------------------------------#
# MISC
#----------------------------------------------------------------------------------------------------------#
def check_range(x, _min, _max, inclusive=True):
    #swap if needed
    if _max < _min:
        _max, _min = _min, _max

    if inclusive == True:
        return _min <= x <= _max
    else:
        return _min < x < _max

def custom_round(v, prec = 1e-8):
    coef = 1 / prec
    return n.round(v * coef) / coef
#----------------------------------------------------------------------------------------------------------#
def clear():
    for i in xrange(0,100):
        print ''
#----------------------------------------------------------------------------------------------------------#
def iterdim(a, axis=0):
  """
  Relevant Stackoverflow:
        http://stackoverflow.com/questions/1589706/iterating-over-arbitrary-dimension-of-numpy-array
  """
  a = numpy.asarray(a)
  leading_indices = (slice(None),)*axis
  for i in xrange(a.shape[axis]) :
    yield a[leading_indices+(i,)]
#----------------------------------------------------------------------------------------------------------#
# TEST-CASE IRB140
#----------------------------------------------------------------------------------------------------------#
class TestIRB140(unittest.TestCase):
        
    def test_elbow_down(self):
        print 'test_elbow_down'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = 40
            j3 = -100
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['T44'], robot_info['robot_geometry_local']

            s = inverse_kinematics_elbow_down(DH_TABLE, A)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_up(DH_TABLE, A)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_fl, delbow_down_fl
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(1, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)


    def test_elbow_up(self):
        print '\ntest_elbow_up'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = 40
            j3 = -30
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['T44'], robot_info['robot_geometry_local']

            s = inverse_kinematics_elbow_up(DH_TABLE, A)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_down(DH_TABLE, A)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_fl, delbow_down_fl
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(0, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)


    def test_elbow_down_backward_facing(self):
        print '\ntest_elbow_down_backward'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = -90
            j3 = -30
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['T44'], robot_info['robot_geometry_local']

            s = inverse_kinematics_elbow_down(DH_TABLE, A, flipped = True)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_up(DH_TABLE, A, flipped = True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_fl, delbow_down_fl
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(3, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)

    def test_elbow_up_backward_facing(self):
        print '\ntest_elbow_up_backward'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = -40
            j3 = -100
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['T44'], robot_info['robot_geometry_local']

            sol = mat([j1,j2,j3,j4,j5,j6])
            s = inverse_kinematics_elbow_up(DH_TABLE, A, flipped = True)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_down(DH_TABLE, A, flipped = True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_backward, delbow_down_backward
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(2, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)

    def test_non_reach_config(self):
            print '\ntest_non_reach_configs'
            for _ in xrange(0,100):
                j1 = rand_range(-180, 180)
                j2 = 90
                j3 = -89
                j4 = rand_range(-200, 200)
                j5 = rand_range(-115, 115)
                j6 = rand_range(-400, 400)

                s0 = j1,j2,j3,j4,j5,j6

                robot_info = forward_kinematics(*s0, **DH_TABLE)
                T44, debug1  = robot_info['T44'], robot_info['robot_geometry_local'] 
                sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
                sol = sol.T

                for i,s in enumerate(sol):
                    robot_info = forward_kinematics(*s, **DH_TABLE)
                    A, debug2  = robot_info['T44'], robot_info['robot_geometry_local'] 
                    if i in [l+m*8 for m,_ in enumerate(range(0, len(sol), 8)) for l in [0,1,4,5]]: #all non-flipped solutions only
                        self.assertAlmostEqual(norm(A-T44), 0)
                    else:
                        self.assertTrue(n.isnan(norm(A-T44)))

    def test_just_barely_reach_flipped_configs(self):
            print '\ntest_just_barely_reach_flipped_configs'
            for _ in xrange(0,100):
                j1 = rand_range(-180, 180)
                j2 = -90
                j3 = -89
                j4 = rand_range(-200, 200)
                j5 = rand_range(-115, 115)
                j6 = rand_range(-400, 400)

                s0 = j1,j2,j3,j4,j5,j6

                robot_info = forward_kinematics(*s0, **DH_TABLE)
                T44, debug1  = robot_info['T44'], robot_info['robot_geometry_local'] 
                sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
                for s in sol.T:
                    robot_info = forward_kinematics(*s, **DH_TABLE)
                    A, debug2  = robot_info['T44'], robot_info['robot_geometry_local'] 
                    self.assertAlmostEqual(norm(A-T44), 0)
                    
    def test_forward_kinematics_general(self):
        print '\ntest_forward_kinematics_general'

        for counter in xrange(10000):
            fcounter = (counter / 10000.0)*100
            if fcounter % 1.0 == 0.0:
                print fcounter
            j1 = rand_range(-180, 180)
            j2 = rand_range(-90, 110)
            j3 = rand_range(-230, 50)
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            # makes sure we never end up at a singular point                

            while (abs(j3) - 90) < 1e-7:
                j3 = rand_range(-230, 50)

            s0 = j1,j2,j3,j4,j5,j6
            robot_info = forward_kinematics(j1, j2, j3, j4, j5, j6, **DH_TABLE)
            T44, debug1 = robot_info['T44'], robot_info['robot_geometry_local']
            
            while norm(calc_wcp(T44,L=0.065)[:2]) < 1e-7:
                j2 = rand_range(-90, 110)
                T44, debug1  = forward_kinematics(j1, j2, j3, j4, j5, j6, **DH_TABLE)
            
            sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
            num_valid_solutions = 0
            for s in sol.T:
                robot_info = forward_kinematics(*s, **DH_TABLE)
                A, debug2  = robot_info['T44'], robot_info['robot_geometry_local']
                num_valid_solutions += check_solution(*s)
                error = norm(A - T44)
                if not n.isnan(error):
                    self.assertAlmostEqual(error, 0)
            self.assertGreaterEqual(num_valid_solutions, 1)
            self.assertEqual(num_valid_solutions, calc_valid_raw_invkin_irb140(T44).shape[0])

            L = []
            for s in iterdim(sol,1):
                if check_solution(*s) == True:
                    L.append(s)
            L = mat(L)
            self.assertTrue(norm(calc_valid_raw_invkin_irb140(T44) - L) == 0.0)
#----------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    unittest.main()
