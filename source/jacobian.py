# imports - core
import sys, os
import numpy as n
n.set_printoptions(suppress=True)

from numpy.linalg import norm, inv, det
from numpy import array as mat, linspace, diff, cross, inf

# unit-testing
import unittest

# imports - thesis
from helperfunctions_math import homogenous_matrix as hom,\
                                 homogenous_matrices as homs,\
                                 rotation_matrix_rot_tilt_skew as mat_rts,\
                                 rot_tilt_skew as rts

from helperfunctions_math import matmul, nmap, homogenous_matrix as hom

from denavit_hartenberg140 import forward_kinematics,\
     inverse_kinematics_irb140, DH_TABLE,\
     inverse_kinematics_curve, calc_valid_invkin_irb140

# functions
def jacobian_from_joints(*joints):
    fk = forward_kinematics(*joints, **DH_TABLE)['robot_geometry_global']
    if not (DH_TABLE['tool'] is None):
        fk[-1] = fk[-1].dot(DH_TABLE['tool'])
    fk = mat(fk)
    zi = fk[:-1,:3,2]
    Jw = zi.T
    
    on = fk[-1,:3,3]
    oi = fk[:-1,:3,3]
    odiff = on - oi
    Jv = nmap(lambda x: reduce(cross, x), zip(zi, odiff)).T
    J = mat([Jv, Jw])
    J = J.reshape((6,6))
    return J

def __tcp_from_joints(*joints, **dh_table):
    return forward_kinematics(*joints, **dh_table)['tcp']

def __flange_from_joints(*joints, **dh_table):
    return forward_kinematics(*joints, **dh_table)['flange']

def __robot_frames(*joints, **dh_table):
    return forward_kinematics(*joints, **dh_table)['robot_geometry_global']
    
def show_func_name(a_func):
    def wrapper(*args, **kwargs):
        print '\n{} :: {}'.format(args[0].__class__.__name__, a_func.__name__)
        return a_func(*args, **kwargs)
    return wrapper

# tests
class JacobianTestWithTool(unittest.TestCase):
    def setUp(self):
        self.dh_table = DH_TABLE
        self.dh_table['tool'] = hom(0,0,0,[0.1,0.12,0.13])
        self.joints = (0,0,0,0,0,0)
        self.tcp = forward_kinematics(*self.joints,**self.dh_table)['tcp']

    @show_func_name
    def test_jacobian_q1_angular_velocity(self):
    	angular_velocities = (1,0,0,0,0,0)
        expected_vw = mat((-0.12,0.645,0,0,0,1))
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q2_angular_velocity(self):
    	angular_velocities = (0,1,0,0,0,0)
        expected_vw = mat([0.26, 0, -0.575, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)
        
    @show_func_name
    def test_jacobian_q3_angular_velocity(self):
    	angular_velocities = (0,0,1,0,0,0)
        expected_vw = mat([-0.1, 0, -0.575, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_q4_angular_velocity(self):
    	angular_velocities = (0,0,0,1,0,0)
        expected_vw = mat([0, 0.1, 0.12, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)

        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_q5_angular_velocity(self):
    	angular_velocities = (0,0,0,0,1,0)
        expected_vw = mat([-0.1, 0, -0.195, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q6_angular_velocity(self):
    	angular_velocities = (0,0,0,0,0,1)
        expected_vw = mat([0, 0.1, 0.12, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_w1_angular_velocity_keep_flange_orientation_constant(self):
    	angular_velocities = (1,0,0,0,0,1)
        expected_vw = mat([ 0, 0.45,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(0,0,0,0,90,0).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_w2_pure_translation_positive_x(self):
    	angular_velocities = (0,1,-1,0,0,0)
        expected_vw = mat([ 0.36, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format(angular_velocities))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_w3_keep_flange_orientation_and_position_constant(self):
    	angular_velocities = (0,0,0,-1,0,1)
        expected_vw = mat([ 0, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format(angular_velocities))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

class JacobianTestWithoutTool(unittest.TestCase):
    def setUp(self):
        self.dh_table = DH_TABLE
        self.dh_table['tool'] = None
        self.joints = (0,0,0,0,0,0)
        self.tcp = forward_kinematics(*self.joints,**self.dh_table)['tcp']

    @show_func_name
    def test_jacobian_q1_angular_velocity(self):
    	angular_velocities = (1,0,0,0,0,0)
        expected_vw = mat((0,0.515,0,0,0,1))
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q2_angular_velocity(self):
    	angular_velocities = (0,1,0,0,0,0)
        expected_vw = mat([0.36, 0, -0.445, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_q3_angular_velocity(self):
    	angular_velocities = (0,0,1,0,0,0)
        expected_vw = mat([0, 0, -0.445, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_q4_angular_velocity(self):
    	angular_velocities = (0,0,0,1,0,0)
        expected_vw = mat([0, 0, 0, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_q5_angular_velocity(self):
    	angular_velocities = (0,0,0,0,1,0)
        expected_vw = mat([0, 0, -0.065, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_q6_angular_velocity(self):
    	angular_velocities = (0,0,0,0,0,1)
        expected_vw = mat([0, 0, 0, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format((angular_velocities)))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_pure_translation_positive_x(self):
    	angular_velocities = (0,1,-1,0,0,0)
        expected_vw = mat([ 0.36, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format(angular_velocities))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_keep_flange_orientation_constant(self):
    	angular_velocities = (1,0,0,0,0,1)
        expected_vw = mat([ 0, 0.45,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(0,0,0,0,90,0).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format(angular_velocities))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

    @show_func_name
    def test_jacobian_keep_flange_orientation_and_position_constant(self):
    	angular_velocities = (0,0,0,-1,0,1)
        expected_vw = mat([ 0, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot(angular_velocities)
        max_error = norm(expected_vw - result_vw, inf)
        print('Joint values: {}'.format(angular_velocities))
        print('Expected: {}\nResult:   {}\nMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)        

class JacobianTestGenericWithoutTool(unittest.TestCase):
    def setUp(self):
        self.dh_table = DH_TABLE
        self.dh_table['tool'] = None
        self.joints = (0,0,0,0,0,0)
        self.tcp = forward_kinematics(*self.joints,**self.dh_table)['tcp']
        self.flange = forward_kinematics(*self.joints,**self.dh_table)['flange']
        
    @show_func_name
    def test_jacobian_determinant_should_be_non_zero(self):
    	joint_values = (10,20,30,40,50,60)
    	print "Joint values: {}".format(joint_values)
        det_val = abs(det(jacobian_from_joints(*joint_values)))
        print det_val
        self.assertNotAlmostEqual(det_val, 0)

    @show_func_name
    def test_jacobian_determinant_should_be_zero(self):
    	joint_values = (0,0,0,0,0,0)
    	print "Joint values: {}".format(joint_values)
        det_val = abs(det(jacobian_from_joints(*joint_values)))
        print det_val
        self.assertAlmostEqual(det_val, 0)
        
    @show_func_name
    def test_jacobian_inversion(self):
        #q1', q2', q3', ...
        joint_values = (1,2,3,4,5,6)
    	print "Joint values: {}".format(joint_values)

        jacobian = jacobian_from_joints(*joint_values)
        vw = jacobian.dot(joint_values)
        print joint_values
        # vx, vy, vz, wx, wy, wz
        q = inv(jacobian).dot(vw)

        for a,b in zip(joint_values, q):
            self.assertAlmostEqual(a,b)

    @show_func_name
    def test_jacobian_calculation_for_non_zero_joints(self):    	
        #q1', q2', q3', ...
        joint_values = (10,20,30,40,50,60)
    	print "Joint values: {}".format(joint_values)
        jacobian = jacobian_from_joints(*joint_values)
        
        res = mat([[-0.10706101, -0.00919023, -0.3423402 ,  0.01752216, -0.0603293 ,
        -0.        ],
       [ 0.4228565 , -0.00162049, -0.06036381,  0.04182162,  0.01663305,
         0.        ],
       [ 0.        , -0.36502331, -0.24189606,  0.02057322,  0.01757034,
         0.        ],
       [ 0.        , -0.17364818, -0.17364818,  0.63302222,  0.35190093,
        -0.12131011],
       [ 0.        ,  0.98480775,  0.98480775,  0.1116189 ,  0.83991154,
         0.47860976],
       [ 1.        ,  0.        ,  0.        , -0.76604444,  0.41317591,
        -0.86960713]])
        self.assertAlmostEqual(n.max(n.abs(res - jacobian)), 0)
        print res

    @show_func_name
    def test_jacobian_calculation_for_zero_joints(self):    	
        #q1', q2', q3', ...
        joint_values = (0,0,0,0,0,0)
    	print "Joint values: {}".format(joint_values)
        jacobian = jacobian_from_joints(*joint_values)
        
        res = mat([[ 0.   ,  0.36 ,  0.   , -0.   , -0.   ,  0.   ],
               [ 0.515, -0.   ,  0.   , -0.   , -0.   ,  0.   ],
               [-0.   , -0.445, -0.445, -0.   , -0.065, -0.   ],
               [ 0.   ,  0.   ,  0.   ,  1.   ,  0.   ,  1.   ],
               [ 0.   ,  1.   ,  1.   , -0.   ,  1.   , -0.   ],
               [ 1.   ,  0.   ,  0.   ,  0.   , -0.   ,  0.   ]])
        self.assertAlmostEqual(n.max(n.abs(res - jacobian)), 0)
        print res



if __name__ == '__main__':
    q = (10,20,30,40,50,60); f = __robot_frames(*q, **DH_TABLE); tcp = __tcp_from_joints(*q, **DH_TABLE)
    unittest.main()
