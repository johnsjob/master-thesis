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
     inverse_kinematics_curve,\
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

    
def __find_solution_path(res, result, curr_point=0, curr_ind=0, tol=20.0):
    p_curr = res[curr_point] #get the solutions for current point
    solution = p_curr[curr_ind]
    if curr_point+1 >= len(res):
        result.append(solution)
        return 
    else:
        p_dest = res[curr_point+1]
    sol_diff = map(norm, solution - p_dest)
    z = zip(p_dest, sol_diff, range(len(p_dest)))
    z_sorted = sorted(z, key=lambda x: x[1])
    z_sorted, _, index = zip(*z_sorted)
    z, _, _ = zip(*z)
    sel_z = z_sorted[0]
    sel_ind = index[0]
    if sol_diff[sel_ind] > tol:
        del result[:]
        return
    else:
        result.append(solution)
    print curr_point
    __find_solution_path(res, result, curr_point = curr_point+1,curr_ind=sel_ind)

def __find_path(ik_curve, index):
    path = []
    __find_solution_path(ik_curve, path,
                         curr_point=0, curr_ind=index)
    return path

def find_paths(ik_curve):
    result = ik_curve
    total = []
    with utils.timing.Timer() as t:
        for index in xrange(len(result[0])):
            path = __find_path(result, index)
            if path:
                print 'FOUND ONE' 
                total.append(list(path))
        return mat(total)

def find_single_path(ik_curve):
    result = ik_curve
    total = []
    with utils.timing.Timer() as t:
        for index in xrange(len(result[0])):
            path = __find_path(result, index)
            if path:
                return mat(path)

if __name__ == '__main__':
    for count in xrange(1):
        dh_table['tool'] = hom(0,-40,0,[0,0,0.05])
        T44 = hom(ori(-90,90,0),[0.5,0,0.5])

        # generate a curve in the last global robot-frame
        num_p = 50
        point_matrix = generate_symmetric_curve(num_points=num_p,
                                                ampl_factor=0.1)
        point_matrix_tf = get_transformed_points(T44, point_matrix)

        # generate angles
        rot  = numpy.linspace(0,0)
        tilt = numpy.linspace(0,0)
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
        
        robot_info = forward_kinematics(*path[37], **dh_table)

        pl = StPlot()
        pl.draw_robot(robot_info['robot_geometry_global'])
        pl.draw_trajectory(transf_frames)
        pl.draw_tool(robot_info['flange'],
                           dh_table['tool'])
        pl.show()
