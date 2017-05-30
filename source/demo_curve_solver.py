import random

import numpy
from numpy.linalg import norm
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

#from pyqtplot import QtPlot
from standardplot import StPlot
import pylab as plt

import itertools as ited

import sys
import time

#test

import utils

##sys.path.append('../int/djikstra/')
##from graph import shortestPath as shortest_path


numpy.set_printoptions(precision=2)
numpy.set_printoptions(suppress=True)

# generate angles
rot  = numpy.linspace(0,180)
tilt = numpy.linspace(-40,40)
skew = numpy.linspace(0,0)
angles = nzip(rot, tilt, skew)

if __name__ == '__main__':
    dh_table['tool'] = hom(0,0,0,[0.0,0.0,0.1])
    wobj = hom(-90,180,0,[0.3,0,0.5])

    # generate a curve in the last global robot-frame
    curve = create_circle_curve(diameter=0.3)

    oris = orientation_frames(angles)
    frames = attach_frames(oris, curve)

    # tansform frames - paper -> robot
    trajectory = attach_to_base_frame(wobj, *frames)
    result = inverse_kinematics_curve(trajectory)
    path = find_single_path(result)
    print path

    for count in xrange(0,50,11):
        print count
        pl = StPlot()
        robot_info = forward_kinematics(*path[count], **dh_table)
        pl.draw_robot(robot_info['robot_geometry_global'])
        pl.draw_trajectory(trajectory)
        pl.draw_tool(robot_info['flange'],
                           dh_table['tool'])
        pl.show()
    pl.draw_joint_paths(path)
