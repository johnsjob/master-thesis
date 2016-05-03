import random

from helperfunctions_plot import *
from pylab import axhline
from plane_relative import *
from denavit_hartenberg140 import *

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


def apply_along_axis(M, func, axis=1):
    return n.apply_along_axis(func, axis, arr=M)

def plot_robot_geometry(ax, global_robot_frames, color='k'):
        for robot_frame in global_robot_frames:
            plot_plane(ax,robot_frame, '--', scale_factor=0.1)
        
        ax.plot(global_robot_frames[:,0,3],
                global_robot_frames[:,1,3],
                global_robot_frames[:,2,3], color, linewidth=2)

def plot_curve(ax, point_matrix):
        ax.scatter(point_matrix[:,0],
                   point_matrix[:,1],
                   point_matrix[:,2])

def plot_robot(ax, color='k', *joint_values):
    T44, debug = forward_kinematics(*joint_values, **DH_TABLE)
    robot_frames = construct_robot_geometry(debug)
    plot_robot_geometry(ax, robot_frames, color)

def construct_robot_geometry(fk_debug_info):
        plane0 = define_plane_from_angles([0,0,0],0, 0, 0)
        global_robot_frames = matmul_series(*fk_debug_info)
        global_robot_frames.insert(0, plane0)
        global_robot_frames = mat( global_robot_frames )
        return global_robot_frames


def plot_robot_from_angles(plot, *args):
    s = forward_kinematics(*args, **DH_TABLE)
    plot.draw_robot(s['robot_geometry_global'])
    return
    
def do_it(res, result, curr_point=0, curr_ind=0, tol=20.0):
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
    do_it(res, result, curr_point = curr_point+1,curr_ind=sel_ind)

if __name__ == '__main__':
    for count in xrange(1):
        j1 =  rand_range(-120,120)
        j2 =  rand_range(-90, 110)
        j3 =  rand_range(-230, 50)
        j4 =  rand_range(-200, 200)
        j5 =  rand_range(-115, 115)
        j6 =  rand_range(-400, 400)

        j1 =  0
        j2 =  0
        j3 =  0
        j4 =  10
        j5 =  20
        j6 =  30

        joint_values = j1,j2,j3,j4,j5,j6

        robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6, **DH_TABLE)
        T44 = robot_info['T44']
        robot_frames = robot_info['robot_geometry_global']

        # generate a curve in the last global robot-frame
        num_p = 50
        point_matrix = generate_symmetric_curve(num_points=num_p,
                                                ampl_factor=0.50)
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
        for index in xrange(31):
            # inverse knematics over curve
            with utils.timing.Timer() as t:
                result = inverse_kinematics_curve(transf_frames)
            print result[0].shape
            res = result
            res = list(res)

            result = []
            with utils.timing.Timer() as t:
                do_it(res, result,curr_point = 0, curr_ind=index)

            if result:
                print 'FOUND ONE!!'
                total.append(list(result))
        print 'total_paths: {}'.format(len(total))
