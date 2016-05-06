import random

from helperfunctions_plot import *
from pylab import axhline
from plane_relative import *
from denavit_hartenberg140 import forward_kinematics, DH_TABLE,\
     inverse_kinematics_curve

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
        j4 =  0
        j5 =  0
        j6 =  0

        joint_values = j1,j2,j3,j4,j5,j6

        robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6, **DH_TABLE)
        T44 = robot_info['tcp']
        robot_frames = robot_info['robot_geometry_global']

        # generate a curve in the last global robot-frame
        num_p = 50
        point_matrix = generate_symmetric_curve(num_points=num_p,
                                                ampl_factor=0.05)
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
        # inverse knematics over curve
        #with utils.timing.Timer() as t:
        result = inverse_kinematics_curve(transf_frames)

        with utils.timing.Timer() as t:
            for index in xrange(len(result[0])):
                _res = []
                __find_solution_path(result, _res,
                                     curr_point=0, curr_ind=index)
                if _res:
                    print 'FOUND ONE!!'
                    total.append(list(_res))
                    #break
        total = mat(total)
        print 'total_paths: {}'.format(len(total))
        pl = StPlot()
        pl.draw_robot(robot_info['robot_geometry_global'])
        pl.draw_trajectory(transf_frames)
        pl.show()
