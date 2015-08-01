from __future__ import division
#--------------------------#
import sys
#--------------------------#
import numpy as n
#--------------------------#
import unittest
#--------------------------#
from helperfunctions_math import *
from denavit_hartenberg import *
from denavit_hartenberg import inverse_kinematics_spherical_wrist as inv_wrist
#=====================================================#
sys.path.append("../int/misc-tools/")
import parsingtools as  parse
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
DH_TABLE = {  'table':[-70, 90, 352, 180,'R',
                     360, 0, 0, 90,'R',
                       0, 90, 0, 180,'R',
                       0, 90, 380, 180,'R',
                       0, 90, 0, 180,'R',
                       0,  0,65,  0,'R'],
             'unit': 'mm',
             'order': ['A','alpha','D','theta'],
             'convention': 'standard'
            }
#----------------------------------------------------------------------------------------------------------#
def inverse_kinematics_elbow_up(dh_table, T44, flipped=False):
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
        if h2-h1 < 0:
            th2 = -(th21-th22)
        else:
            th2 = th21 + th22
        j2 = 90 - th2
        #import pdb; pdb.set_trace()
    else:
        ### elbow-up (actually inverse, upside-down) ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(90 - th3)
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        if h2-h1 < 0:
            th2 = -(th21 + th22)
        else:
            th2 = th21 - th22
        j2 = -(90-th2)
    j4, j5, j6,\
        j41,j51,j61 = inverse_kinematics_spherical_wrist(dh_table, j1, j2, j3, T44)
    return (j1, j2, j3, j4, j5, j6), (j1, j2, j3, j41, j51, j61)

def inverse_kinematics_elbow_down(dh_table, T44, flipped=False):
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
        if h2-h1 < 0:
            th2 = -(th21+th22)
        else:
            th2 =  th21 - th22
        j2 = 90 - th2
        #import pdb; pdb.set_trace()
    else:
        ### elbow-down (actually inverse, upside-down) ###
        # Third angle - j3
        th3 = ang_sats2(x1, alpha, beta)
        j3 = -(th3 + 90)
        # Second angle - j2
        th21 = atan2(s, x0)
        th22 = atan2(beta*sin2(th3), alpha + beta*cos2(th3))
        if h2-h1 < 0:
            th2 =  -(th21 - th22)
        else:
            th2 =  th21 + th22
        j2 = -(90 - th2)

    j4, j5, j6,\
        j41,j51,j61 = inverse_kinematics_spherical_wrist(dh_table, j1, j2, j3, T44)
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

def inverse_kinematics_irb140(dh_table, T44):
    if type(T44) is list:
        T44 = mat(T44)
    dim = T44.shape
    if len(dim) != 2:
        raise ArithmeticError('Forward-kinematics must be a 4x4 matrix!')
    if dim[0] != dim[1]:
        raise ArithmeticError('Forward-kinematics must be square!')
    if dim[0] != 4:
        raise ArithmeticError('Forward-kinematics must have dimension of 4!')

    sol_elbup1,      sol_elbup2      = inverse_kinematics_elbow_up(dh_table, T44)
    sol_elbdown1,    sol_elbdown2    = inverse_kinematics_elbow_down(dh_table, T44)
    sol_elbup1_fl,   sol_elbup2_fl   = inverse_kinematics_elbow_up(dh_table, T44, flipped = True)
    sol_elbdown1_fl, sol_elbdown2_fl = inverse_kinematics_elbow_down(dh_table, T44, flipped = True)

    #first columnt is first solution and so forth
    return mat(zip(sol_elbup1, sol_elbdown1, sol_elbup1_fl, sol_elbdown1_fl,
                   sol_elbup2, sol_elbdown2, sol_elbup2_fl, sol_elbdown2_fl))

def filter_solutions(solutions, filter_function = check_solution):
    result = []
    for s in solutions.T:
        if filter_function(*s) == True:
            result.append( s )
    return mat(zip(*result))

def calc_valid_inv_kin_IRB140(T44):
    return filter_solutions( inverse_kinematics_irb140(T44) )

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
class TestIRB140(unittest.TestCase):
        
    def test_forward_kinematics_from_file(self):
        data = parse.parse_file("C:\\robot-studio-output.txt")
        num_joint_conf = len(data['Joint_1_T'])

        self.data = data
        self.num_joint_conf = num_joint_conf
        self.T44 = mat(data['T44'][:]).reshape((num_joint_conf,4,4))

        a,b,c,d,e,f = data['Joint_1'][:], data['Joint_2'][:],\
                      data['Joint_3'][:], data['Joint_4'][:],\
                      data['Joint_5'][:], data['Joint_6'][:]
        self.joint_values = zip(a,b,c,d,e,f)

        data = self.data
        num_joint_conf = self.num_joint_conf
        for index in xrange(num_joint_conf):
            #T44 = mat(data['T44'][index]).reshape((4,4))
            T44 = self.T44[index]

##            for key in data:
##                if 'T' in key:
##                    data[key] = mat(data[key]).reshape(num_joint_conf,4,4)

            a,b,c,d,e,f = self.joint_values[index]
            A, debug  = forward_kinematics(a,b,c,d,e,f, **DH_TABLE)

            print "T44 sanity check-norm: " + str(norm(T44 - A))
            self.assertLess(norm(T44 - A), 1e-7)
        
    def test_inverse_kinematics_from_file(self):
        print '\ntest_inverse_kinematics'
        data = parse.parse_file("C:\\robot-studio-output.txt")
        num_joint_conf = len(data['Joint_1_T'])

        self.data = data
        self.num_joint_conf = num_joint_conf
        self.T44 = mat(data['T44'][:]).reshape((num_joint_conf,4,4))
        print "\nNumber of configurations: " + str(self.num_joint_conf) + "\n"

        a,b,c,d,e,f = data['Joint_1'][:], data['Joint_2'][:],\
                      data['Joint_3'][:], data['Joint_4'][:],\
                      data['Joint_5'][:], data['Joint_6'][:]
        self.joint_values = zip(a,b,c,d,e,f)

        print '### TEST-inverse-kinematics-from-file'
        nans = 0
        total_valid_sols = 0
        num_zeroed_angles = 0
        for index in xrange(self.num_joint_conf):
            T44 = self.T44[index]
            sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )

            a,b,c,d,e,f = self.joint_values[index]
            s0 = mat([a,b,c,d,e,f])
            print "\n[ IK %s ]" % str(index)

            num_zeroed_angle_norms = 0
            num_valid_solutions = 0
            for i in xrange(0, 8):
                solution_names = ['sol_elbup1', 'sol_elbdown1', 'sol_elbup1_fl', 'sol_elbdown1_fl',
                                  'sol_elbup2', 'sol_elbdown2', 'sol_elbup2_fl', 'sol_elbdown2_fl']
                s = sol[:,i]
                is_valid_solution = check_solution(*s)
                num_valid_solutions += is_valid_solution
                
                gamma0,gamma1,gamma2,gamma3,gamma4,gamma5 = s
                A, debug = forward_kinematics(gamma0, gamma1, gamma2,
                                                 gamma3, gamma4, gamma5, **DH_TABLE)
                p0 = debug[0][:,3]
                p1 = matmul(debug[0],debug[1])[:,3]
                p2 = matmul(debug[0],debug[1],debug[2])[:,3]
                p3 = matmul(debug[0],debug[1],debug[2], debug[3])[:,3]
                p4 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4])[:,3]
                p5 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4], debug[5])[:,3]
                fk_norm = norm(A - T44)
                angle_norm = norm(s - s0)

                print "\n\tSolution: %s" % solution_names[i]
                print "\tFK-norm: " + str(fk_norm)
                if not n.isnan(fk_norm):
                    self.assertLess(fk_norm, 1e-14)
                else:
                    nans += 1
                    
                print "\tangle-norm: %s" % str(angle_norm)
                if not n.isnan(fk_norm):
                    if angle_norm < 1e-8:
                        num_zeroed_angle_norms += 1
                    elif int(angle_norm*1000) == 360000:
                        num_zeroed_angle_norms += 1
                    elif angle_norm < 1:
                        print str(s - s0)
            self.assertGreater(num_zeroed_angle_norms, 0)

            print "\n\tnum valid solutions: %s" % str(num_valid_solutions)
            print "\tnum angle norms ~0: %s" % str(num_zeroed_angle_norms)
            num_zeroed_angles += num_zeroed_angle_norms
            total_valid_sols += num_valid_solutions
            self.assertGreater(num_valid_solutions, 0)
        print 'nans = '+str(nans)
        print 'total_valid = '+str(total_valid_sols)
        print 'total_possible = '+str(self.num_joint_conf*8)
        print 'total_zeroed_angles = '+str(num_zeroed_angles)

    def test_elbow_down(self):
        print 'test_elbow_down'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = 40
            j3 = -100
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)
            A, debug = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)

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


    def test_elbow_up(self):
        print '\ntest_elbow_up'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = 40
            j3 = -30
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)
            A, debug = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)

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

    def test_elbow_up_flipped(self):
        print '\ntest_elbow_up_flipped'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = -90
            j3 = -30
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)
            A, debug = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)

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

    def test_elbow_down_flipped(self):
        print '\ntest_elbow_down_flipped'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = -40
            j3 = -100
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)
            A, debug = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)

            sol = mat([j1,j2,j3,j4,j5,j6])
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

                T44, debug1  = forward_kinematics(*s0, **DH_TABLE)
                sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
                for i,s in enumerate(sol.T):
                    A, debug2  = forward_kinematics(*s, **DH_TABLE)
                    if i in [0,1,4,5]: #all non-flipped solutions only
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

                T44, debug1  = forward_kinematics(*s0, **DH_TABLE)
                sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
                for s in sol.T:
                    A, debug2  = forward_kinematics(*s, **DH_TABLE)
                    self.assertAlmostEqual(norm(A-T44), 0)

    def test_forward_kinematics_general(self):
        print '\ntest_forward_kinematics_general'

        for _ in xrange(100):
            j1 = rand_range(-180, 180)
            j2 = rand_range(-90, 110)
            j3 = rand_range(-230, 50)
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            # makes sure we never end up at a singular point
            if abs(90 - j2) <= 1e-7:
                if abs(j2 - j3) <= 1e-7:
                    j3 = j3 + n.sign(rand_range(-1, 1)) * rand_range(0.5, 1)

            if abs(j2) <= 1e-7:
                if abs(j2 - j3) <= 1e-7:
                    j3 = j3 + n.sign(rand_range(-1, 1)) * rand_range(0.5, 1)

            s0 = j1,j2,j3,j4,j5,j6
            T44, debug1  = forward_kinematics(j1, j2, j3, j4, j5, j6, **DH_TABLE)
            
            sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )

            num_valid_solutions = 0
            for s in sol.T:
                A, debug2  = forward_kinematics(*s, **DH_TABLE)
                num_valid_solutions += check_solution(*s)
                error = norm(A - T44)
                if not n.isnan(error):
                    try:
                        self.assertAlmostEqual(error, 0)
                    except Exception:
                        import pdb; pdb.set_trace()
            self.assertGreaterEqual(num_valid_solutions, 1)
#----------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
##    unittest.main()
##    """
##    GENERAL COMMENTS:
##    ################
##        data['Joint_(i)_T'] = cumulative multiplication of up to i-1 of joint T-matrices
##
##        'debug' contains i:th matrices of python implementation of DH-algorithm
##
##        lambda function up_to(i) performs cumulative multiplication up to i:th joint
##
##        wcp = data['Joint_6_T'][index] = up_to(5)
##    """
##    clear()
##    index = -1
##
##    data = parse.parse_file("C:\\robot-studio-output.txt")
##    num_joint_conf = len(data['Joint_1_T'])
##            
##            
##    params = zip(data['A'], data['alpha'], data['D'], data['theta'])
##    T44 = mat(data['T44'][index]).reshape((4,4))
##
##    for key in data:
##        if 'T' in key:
##            data[key] = custom_round(mat(data[key]).reshape(num_joint_conf,4,4))
##
##    print "\nNumber of configurations: " + str(len(data['Joint_1'])) + "\n"
##    a,b,c,d,e,f = data['Joint_1'][index], data['Joint_2'][index], data['Joint_3'][index], data['Joint_4'][index], data['Joint_5'][index], data['Joint_6'][index],

    j1 = -101.19523566942526
    j2 =  67.81576827839055
    j3 = -222.64836963410147
    j4 =  157.78092915244753
    j5 =  109.07390106077256
    j6 = -95.98652095283182
    a,b,c,d,e,f = j1,j2,j3,j4,j5,j6


    T44, debug  = forward_kinematics(a,b,c,d,e,f, **DH_TABLE)

    #print "T44 sanity check-norm: " + str(norm(T44 - A))

    p_end = T44[0:3,3]
    wcp = calc_wcp(T44, 0.065)
    sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
    s0 = mat([a,b,c,d,e,f])
    all_norms = 0

    solution_names = ['sol_elbup1', 'sol_elbdown1', 'sol_elbup1_fl', 'sol_elbdown1_fl',
                      'sol_elbup2', 'sol_elbdown2', 'sol_elbup2_fl', 'sol_elbdown2_fl']

    num_valid_solutions = 0
    for i in xrange(0, 8):
##        if not ((i == 1)):
##            continue 
        s = sol[:,i]
        num_valid_solutions += check_solution(*s)
        print check_solution(*s)
        print s
        gamma0,gamma1,gamma2,gamma3,gamma4,gamma5 = s
        A, debug = forward_kinematics(gamma0, gamma1, gamma2,
                                         gamma3, gamma4, gamma5, **DH_TABLE)
        A = A
        p0 = debug[0][:,3]
        p1 = matmul(debug[0],debug[1])[:,3]
        p2 = matmul(debug[0],debug[1],debug[2])[:,3]
        p3 = matmul(debug[0],debug[1],debug[2], debug[3])[:,3]
        p4 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4])[:,3]
        p5 = matmul(debug[0],debug[1],debug[2], debug[3], debug[4], debug[5])[:,3]

        print "\n\t[ Solution %s ]" % solution_names[i]
        print "\tFK-norm: " + str( norm(A - T44) )
        all_norms = all_norms + norm(A - T44)
        print "\tangle-norm: %0.2f" % norm(s - s0)

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
    print "Num valid solutions: " + str( num_valid_solutions )
    plot([-1,-1,1,1],[0,1,1,0],'w')
    plot(wcp[0], wcp[2], 'ro')
    plot(p_end[0], p_end[2], 'ko')
    show()
