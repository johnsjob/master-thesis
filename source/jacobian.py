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

def tcp_from_joints(*joints, **dh_table):
    return forward_kinematics(*joints, **dh_table)['tcp']

def flange_from_joints(*joints, **dh_table):
    return forward_kinematics(*joints, **dh_table)['flange']

def robot_frames(*joints, **dh_table):
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
    def test_forward_kinematics_setup(self):
        tcp = n.round(self.tcp, 12)
        test_val = ( tcp == mat([[ 0. ,  0. ,  1. ,  0.645],
                                 [-0. ,  1. , -0. ,  0.12 ],
                                 [-1. , -0. ,  0. ,  0.612],
                                 [ 0. ,  0. ,  0. ,  1.   ]]).reshape(4,4)).all()
        self.assertTrue(test_val)
        print(tcp)

    @show_func_name
    def test_jacobian_q1(self):
        expected_vw = mat((-0.12,0.645,0,0,0,1))
        result_vw  = jacobian_from_joints(*self.joints).dot((1,0,0,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q2(self):
        expected_vw = mat([0.26, 0, -0.575, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,1,0,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q3(self):
        expected_vw = mat([-0.1, 0, -0.575, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,1,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q4(self):
        expected_vw = mat([0, 0.1, 0.12, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,1,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q5(self):
        expected_vw = mat([-0.1, 0, -0.195, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,0,1,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q6(self):
        expected_vw = mat([0, 0.1, 0.12, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,0,0,1))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_w1_keep_flange_orientation_constant(self):
        expected_vw = mat([ 0, 0.45,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(0,0,0,0,90,0).dot((1,0,0,0,0,1))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_w2_pure_translation_positive_x(self):
        expected_vw = mat([ 0.36, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,1,-1,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_w3_keep_flange_orientation_and_position_constant(self):
        expected_vw = mat([ 0, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,-1,0,1))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
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
    def test_forward_kinematics_setup(self):
        tcp = n.round(self.tcp, 12)
        test_val = ( tcp == mat([[ 0. ,  0. ,  1. ,  0.515],
                                 [-0. ,  1. , -0. , -0.   ],
                                 [-1. , -0. ,  0. ,  0.712],
                                 [ 0. ,  0. ,  0. ,  1.   ]]).reshape(4,4)).all()
        self.assertTrue(test_val)
        print(tcp)

    @show_func_name
    def test_jacobian_q1(self):
        expected_vw = mat((0,0.515,0,0,0,1))
        result_vw  = jacobian_from_joints(*self.joints).dot((1,0,0,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q2(self):
        expected_vw = mat([0.36, 0, -0.445, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,1,0,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q3(self):
        expected_vw = mat([0, 0, -0.445, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,1,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q4(self):
        expected_vw = mat([0, 0, 0, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,1,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q5(self):
        expected_vw = mat([0, 0, -0.065, 0, 1, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,0,1,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_q6(self):
        expected_vw = mat([0, 0, 0, 1, 0, 0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,0,0,1))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_pure_translation_positive_x(self):
        expected_vw = mat([ 0.36, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,1,-1,0,0,0))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_keep_flange_orientation_constant(self):
        expected_vw = mat([ 0, 0.45,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(0,0,0,0,90,0).dot((1,0,0,0,0,1))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
                                                               result_vw,
                                                               max_error))
        self.assertAlmostEqual(0, max_error)

    @show_func_name
    def test_jacobian_keep_flange_orientation_and_position_constant(self):
        expected_vw = mat([ 0, 0,  0,  0,  0,  0])
        result_vw  = jacobian_from_joints(*self.joints).dot((0,0,0,-1,0,1))
        max_error = norm(expected_vw - result_vw, inf)
        print('\tExpected: {}\n\tResult:   {}\n\tMax Error: {}\n'.format(expected_vw,
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
    def test_fordward_kinematics_setup(self):
        tcp = n.round(self.tcp,12)
        test_val = ( tcp == mat([ 0, 0, 1, 0.515,
                                  0, 1, 0, 0,
                                 -1, 0, 0, 0.712,
                                  0, 0, 0, 1]).reshape(4,4)).all()
        self.assertTrue(test_val)
        print 'tcp = \n{}\n'.format(tcp)

    @show_func_name
    def test_jacobian_determinant_should_be_non_zero(self):
        det_val = abs(det(jacobian_from_joints(10,20,30,40,50,60)))
        print det_val
        self.assertNotAlmostEqual(det_val, 0)

    @show_func_name
    def test_jacobian_determinant_should_be_zero(self):
        det_val = abs(det(jacobian_from_joints(0,0,0,0,0,0)))
        print det_val
        self.assertAlmostEqual(det_val, 0)

    @show_func_name
    def test_jacobian_inversion(self):
        #q1', q2', q3', ...
        test_joints = (1,2,3,4,5,6)
        jacobian = jacobian_from_joints(*test_joints)
        vw = jacobian.dot(test_joints)
        print test_joints
        # vx, vy, vz, wx, wy, wz
        q = inv(jacobian).dot(vw)
        print q
        for a,b in zip(test_joints, q):
            self.assertAlmostEqual(a,b)

if __name__ == '__main__':
    q = (0,0,0,0,0,0); f = robot_frames(*q, **DH_TABLE); tcp = tcp_from_joints(*q, **DH_TABLE)
    unittest.main()
