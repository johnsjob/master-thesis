import random

import numpy, numpy as n
from numpy import pi, linspace
from numpy.linalg import norm, det, inv
from helperfunctions_plot import *
from helperfunctions_math import rand_range,\
                                 rotation_matrix_rot_tilt_skew as ori,\
                                 homogenous_matrices, nzip, nmap,\
                                 quat_slerp
from pylab import axhline
from plane_relative import generate_symmetric_curve,\
                           generate_curve,\
                           get_transformed_points, attach_to_base_frame,\
                           create_circle_curve, place_curve, attach_frames,\
                           orientation_frames

from denavit_hartenberg140 import forward_kinematics, DH_TABLE as dh_table,\
     inverse_kinematics_curve, find_single_path, \
     calc_valid_invkin_irb140 as invkin

from denavit_hartenberg import homogenous_matrix as hom
from jacobian import jacobian_from_joints

#from pyqtplot import QtPlot
from standardplot import StPlot
import pylab as plt

import itertools as it

import sys
import time

#test

import utils

##sys.path.append('../int/djikstra/')
##from graph import shortestPath as shortest_path


numpy.set_printoptions(precision=4)
numpy.set_printoptions(suppress=True)

# generate angles
num_points = 50
rot  = numpy.linspace(0,0,num_points)
tilt = numpy.linspace(0,0,num_points)
skew = numpy.linspace(0,0,num_points)
angles = nzip(rot, tilt, skew)

normalize = lambda x: x / norm(x)

def calc_robot_tcp(*joints):
    return forward_kinematics(*joints, **dh_table)['tcp']

def calc_robot_curve_j1(*jvals):
    return mat([calc_robot_tcp(v,0,0,0,0,0) for v in jvals])


if __name__ == '__main__':
    dh_table['tool'] = hom(0,0,0,[0.0,0.0,0.1])
    wobj = hom(-90,180,0,[0.6,0,0.5])
    # robot movement
    robot_movement_j1 = calc_robot_curve_j1(*linspace(0,90,50))

    trajectory = robot_movement_j1
    
    # create velocities
    # traverse curve in 10 seconds
    velocity = n.diff(trajectory[:,:3,3], axis=0) / (10.0/len(trajectory)) #m/s
    velocity = mat([[0,0,0]] + list(velocity))
    angular_velocity = velocity*0
    vw = nzip(velocity, angular_velocity).reshape(50,6)

    #inverse kinematics over curve
    result = inverse_kinematics_curve(trajectory)
    path = find_single_path(result)
##    joint_vel = nmap(lambda x: calc_joint_velocities(*x), zip(path, vw))
    J = nmap(lambda x: jacobian_from_joints(*x), path)
    Jinv = nmap(inv, J)
    joint_angular_vel = nmap(lambda x: reduce(n.dot, x), zip(Jinv, vw))

## 
##    print 'Path found: \n {}'.format(path)
##    print 'Joint angular velocities: \n {}'.format(joint_vel)
##
##    for count in xrange(1):
##        pl = StPlot()
####        joints = path[count]
####                
####        robot_info = forward_kinematics(*joints, **dh_table)
####        pl.draw_robot(robot_info['robot_geometry_global'])
####        pl.draw_trajectory(trajectory)
####        pl.draw_tool(robot_info['flange'],
####                           dh_table['tool'])
####        pl.draw_frame(wobj, size=0.2)
####        break
##    pl.draw_joint_paths(path)
##    pl.draw_joint_velocities(joint_vel)
###    pl.draw_jacobian_determinants(J)
##    pl.show()
