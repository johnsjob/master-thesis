# imports - core
import sys, os
import numpy as n
n.set_printoptions(suppress=True)

from numpy.linalg import norm, inv, det
from numpy import array as mat, linspace, diff, cross

# unit-testing
import unittest
import inspect

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

# tests
class JacobianTest(unittest.TestCase):
    def setUp(self):
        self.dh_table = DH_TABLE
        self.dh_table['tool'] = hom(0,0,0,[-0.2,0.0,0])
        self.joints = (0,0,0,
                      0,90,0)
        self.tcp = forward_kinematics(*self.joints,**self.dh_table)['tcp']
        self.flange = forward_kinematics(*self.joints,**self.dh_table)['flange']
        self.invkin_joints = calc_valid_invkin_irb140(self.tcp)
        self.jacobian = jacobian_from_joints(*self.invkin_joints[0])
        

    def test_fordward_kinematics_setup(self):
        print inspect.stack()[0][3]
        tcp = n.round(self.tcp,12)
        test_val = ( tcp == mat([-1, 0, 0, 0.65,
                                  0, 1, 0, 0,
                                  0, 0,-1, 0.647,
                                  0, 0, 0, 1]).reshape(4,4)).all()
        self.assertTrue(test_val)
        print tcp

        flange = n.round(self.flange,12)
        test_val = ( flange == mat([-1, 0, 0, 0.45,
                                     0, 1, 0, 0,
                                     0, 0,-1, 0.647,
                                     0, 0, 0, 1]).reshape(4,4)).all()
        self.assertTrue(test_val)
        print flange

    def test_jacobian_determinant(self):
        print inspect.stack()[0][3]
        det_val = abs(det(self.jacobian))
        print det_val
        self.assertNotAlmostEqual(det_val, 0)

    def test_jacobian_inversion(self):
        print inspect.stack()[0][3]
        #j1', j2', j3', ...
        test_joints = (1,2,3,4,5,6)
        vw = self.jacobian.dot(test_joints)
        print test_joints
        # vx, vy, vz, wx, wy, wz
        q = inv(self.jacobian).dot(vw)
        print q
        for a,b in zip(test_joints, q):
            self.assertAlmostEqual(a,b)


    def test_rotation_with_non_rigid_link(self):
        print inspect.stack()[0][3]
        #j1', j2', j3', ...
        vw = self.jacobian.dot((1,0,0,
                                0,0,1))
        print vw
        for i in [0,2,3,4,5]:
            self.assertAlmostEqual(vw[i],0)
        self.assertAlmostEqual(vw[1], self.flange[0,3])

    def test_rotation_with_rigid_link(self):
        print inspect.stack()[0][3]
        #j1', j2', j3', ...
        vw = self.jacobian.dot((1,0,0,
                                0,0,0))
        print vw
        for i in [0,2,3,4]:
            self.assertAlmostEqual(vw[i],0)
        self.assertAlmostEqual(vw[-1], 1)
        self.assertAlmostEqual(vw[1], self.tcp[0,3])


if __name__ == '__main__':
    unittest.main()
