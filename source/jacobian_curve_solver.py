import random

import numpy, numpy as n
from numpy import pi
from numpy.linalg import norm, det, inv
from helperfunctions_plot import *
from helperfunctions_math import rand_range,\
                                 rotation_matrix_rot_tilt_skew as ori,\
                                 homogenous_matrices, nzip, nmap
from pylab import axhline
from plane_relative import generate_symmetric_curve,\
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
rot  = numpy.linspace(0,180,num_points)
tilt = numpy.linspace(-40,40,num_points)
skew = numpy.linspace(0,0,num_points)
angles = nzip(rot, tilt, skew)

normalize = lambda x: x / norm(x)

def calc_jacobian(*joints):
    J = jacobian_from_joints(*joints)
    return J

def calc_joint_velocities(joints, vw):
    J = calc_jacobian(*joints)
    joint_vel = inv(J).dot(vw)
    return joint_vel * 180.0 / pi

if __name__ == '__main__':
    dh_table['tool'] = hom(0,0,0,[0.0,0.0,0.1])
    wobj = hom(-90,180,0,[0.3,0,0.5])

    # generate a curve in the last global robot-frame
    curve = create_circle_curve(diameter=0.3, num_p=num_points)

    oris = orientation_frames(angles)
    frames = attach_frames(oris, curve)

    # tansform frames - paper -> robot
    trajectory = attach_to_base_frame(wobj, *frames)
    center = wobj[:3,3]
    r = trajectory[:,:3,3] - center
    normal = normalize( n.cross(r[0], r[12]) )
    # traverse curve in 10 seconds
    # speed = circumf / 10s = 0.3*pi / 10s
    speed = 0.0942477796076937
    w = normal * speed
    velocity = n.cross(w,r) * speed
    velocity = nmap(lambda x: list(x)+[0,0,0], velocity)
    result = inverse_kinematics_curve(trajectory)
    path = find_single_path(result)
    joint_vel = nmap(lambda x: calc_joint_velocities(*x), zip(path, velocity))
    J = nmap(lambda x: calc_jacobian(*x), path)
    print path
    print joint_vel

    for count in xrange(0,50,11):
        print count
        pl = StPlot()
        joints = path[count]
                
        robot_info = forward_kinematics(*joints, **dh_table)
        pl.draw_robot(robot_info['robot_geometry_global'])
        pl.draw_trajectory(trajectory)
        pl.draw_tool(robot_info['flange'],
                           dh_table['tool'])
        break
    pl.draw_joint_paths(path)
    pl.draw_joint_velocities(joint_vel)
    pl.draw_jacobian_determinants(J)
    pl.show()
