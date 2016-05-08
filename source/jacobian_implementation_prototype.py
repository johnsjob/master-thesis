# imports - core
import sys, os

import numpy as n
from numpy.linalg import norm, inv, det
from numpy import array as mat, linspace, diff, cross
n.set_printoptions(suppress=True)

from random import random as rand
# imports - anoto

# imports - thesis
sys.path.append('../int/master_thesis_code/source/')
from helperfunctions_math import homogenous_matrix as hom,\
                                 homogenous_matrices as homs,\
                                 rotation_matrix_rot_tilt_skew as mat_rts,\
                                 rot_tilt_skew as rts
from helperfunctions_math import matmul, nmap, homogenous_matrix as hom

from denavit_hartenberg140 import forward_kinematics,\
     inverse_kinematics_irb140, DH_TABLE,\
     inverse_kinematics_curve, calc_valid_invkin_irb140
import denavit_hartenberg140 as robot

from calibration_algorithm import _solve_orientation as SVD
from plane_relative import generate_curve

# globals
coord_conv = mat([[0,1,0], [1,0,0], [0,0,-1]]).T
tool0 = hom(0,0,0,
            0,0,0)

def jacobian_from_joints(*joints):
    fk = forward_kinematics(*joints, **DH_TABLE)['robot_geometry_global']
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

if __name__ == '__main__':
    joints = (0,0,0,
              0,90,0)
    DH_TABLE['tool'] = hom(0,0,0,[-0.2,0.0,0])
    pose = forward_kinematics(*joints,**DH_TABLE)['tcp']
    flange = forward_kinematics(*joints,**DH_TABLE)['flange']
    #pose = hom(-90,90,0, [0.5,0,0.5])
    res = calc_valid_invkin_irb140(pose)
    j = jacobian_from_joints(*res[0])
    #j1', j2', j3', ...
    vw = j.dot((1,0,0,
                0,0,1))
    print vw
    # vx, vy, vz, wx, wy, wz
    q = inv(j).dot((0,1,0,
                    0,0,1))
    print q
    
    print 'FINISHED'     
