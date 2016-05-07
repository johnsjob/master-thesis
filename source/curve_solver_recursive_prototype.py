import random

import numpy
from numpy.linalg import norm
from helperfunctions_plot import *
from helperfunctions_math import rand_range, rotation_matrices,\
                                 rotation_matrix_rot_tilt_skew as ori,\
                                 homogenous_matrices
from pylab import axhline
from plane_relative import generate_symmetric_curve,\
                           get_transformed_points, apply_transform_on_frames

from denavit_hartenberg140 import forward_kinematics, DH_TABLE as dh_table,\
     inverse_kinematics_curve, find_single_path, \
     calc_valid_invkin_irb140 as invkin

from denavit_hartenberg import homogenous_matrix as hom

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


numpy.set_printoptions(precision=2)
numpy.set_printoptions(suppress=True)

if __name__ == '__main__':
    for count in xrange(0,50,11):
        dh_table['tool'] = hom(0,0,0,[0.0,0,0.1])
        T44 = hom(-90,180,0,[0.3,0,0.5])

        # generate a curve in the last global robot-frame
        num_p = 50
        point_matrix = generate_symmetric_curve(num_points=num_p,
                                                ampl_factor=0.3)
        point_matrix_tf = get_transformed_points(T44, point_matrix)

        # generate angles
        rot  = numpy.linspace(0,180)
        tilt = numpy.linspace(-40,40)
        skew = numpy.linspace(0,0)

        # generate frames
        angles = zip(rot, tilt, skew)
        R = rotation_matrices(angles)
        
        pre_frames = zip(R, point_matrix)
        frames = homogenous_matrices(pre_frames)

     # tansform frames - paper -> robot
        transf_frames = apply_transform_on_frames(T44, frames)

        total = []
        total_time = 0

        result = inverse_kinematics_curve(transf_frames)
        path = find_single_path(result)
        print path
        
        robot_info = forward_kinematics(*path[count], **dh_table)

        pl = StPlot()
        pl.draw_robot(robot_info['robot_geometry_global'])
        pl.draw_trajectory(transf_frames)
        pl.draw_tool(robot_info['flange'],
                           dh_table['tool'])
        pl.show()
