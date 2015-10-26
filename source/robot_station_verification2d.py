import random

from helperfunctions_plot import *
from pylab import axhline
from plane_relative import *
from denavit_hartenberg140 import *

import itertools as it

def apply_along_axis(M, func, axis=1):
    return n.apply_along_axis(func, axis, arr=M)

def plot_robot_geometry(ax, global_robot_frames, color='k'):
        
        ax.plot(global_robot_frames[:,0,3],
                global_robot_frames[:,2,3], color ,linewidth=2)

def plot_curve(ax, point_matrix):
        ax.scatter(point_matrix[:,0],
                   point_matrix[:,1],
                   point_matrix[:,2])

def construct_robot_geometry(fk_debug_info):
        plane0 = define_plane_from_angles([0,0,0],0, 0, 0)
        global_robot_frames = matmul_series(*fk_debug_info)
        global_robot_frames.insert(0, plane0)
        global_robot_frames = mat( global_robot_frames )
        return global_robot_frames

def get_closest_solutions_pair(s0, s1, norm_func,**kwargs):
    data = []
    for i, s0i in enumerate(s0):
        for j, s1j in enumerate(s1):
####            print norm_func(s0i - s1j, **kwargs)
            data.append([norm_func(s0i - s1j, **kwargs), i, j])
####    print ''
    data = mat(data)

    ret = []
    solution_col_row_pairs = n.argwhere(data == data.min(axis = 0)[0])
    solution_indices = solution_col_row_pairs[:,0]
####    print len(data[solution_indices])
    for solution_data in data[solution_indices]:
        norm_value, i, j = solution_data
        pair = mat([s0[i], s1[j]])
####        print 'small: ' + str(n.linalg.norm(n.diff(pair, axis=0)))
        return pair

def merge_solutions(*args):
    result = []
    for m in args:
        result += zip(*m)
    return mat(zip(*result))

def __modulo_solutions(solution_matrix, index, modulo=360.0):
    for s in solution_matrix.T:
        result = s.copy()
        value = result[index]
        result[index] = value + modulo
        yield result

def generate_modulo_solutions(solution_matrix, index, modulo=360.0):
    return mat(zip(*__modulo_solutions(solution_matrix, index, modulo)))

def extract_closest_solutions(all_solutions, norm_func, **kwargs):
        # obtain closest solutions
        chosen_solutions = []
        num_solutions = len(all_solutions)
        for k in xrange(num_solutions-1):
####            print 'INDEX: ' + str(k)
            if k == 0:
                o = all_solutions[k]
            else:
                o = chosen_solutions[-1]

            pair = get_closest_solutions_pair(o, all_solutions[k+1], norm_func, **kwargs)

            if k == 0:
                chosen_solutions.append(mat([pair[0]]))
                chosen_solutions.append(mat([pair[1]]))
            else:
                chosen_solutions.append(mat([pair[1]]))

        chosen_solutions = mat(chosen_solutions).reshape(num_solutions, 6)
        return chosen_solutions
def my_norm(x, **kwargs):
    max_num = len(x)
    factors = mat([max_num-k for k in xrange(max_num)])
    return norm(factors*x)
        
if __name__ == '__main__':
    for count in xrange(10):
        ax, fig = init_plot()
        fig.clear()

        j1 =  0
        j2 =  180
        j3 =  -90+5
        j4 =  0
        j5 =  0
        j6 =  0

        joint_values = mat([j1,j2,j3,j4,j5,j6])

        # get forward kinematics i.e. last global robot-frame
        T44, debug = forward_kinematics(*joint_values, **DH_TABLE)
        IK_angles = inverse_kinematics_irb140(DH_TABLE, T44)

        elbow_up = mat(zip(IK_angles[:,0],IK_angles[:,4],IK_angles[:,8]))
        elbow_down = mat(zip(IK_angles[:,0+1],IK_angles[:,4+1],IK_angles[:,8+1]))
        elbow_up_fl = mat(zip(IK_angles[:,0+2],IK_angles[:,4+2],IK_angles[:,8+2]))
        elbow_down_fl = mat(zip(IK_angles[:,0+3],IK_angles[:,4+3],IK_angles[:,8+3]))
        IK_angles = mat(zip(elbow_up, elbow_down, elbow_up_fl, elbow_down_fl)).reshape(6,12)
        mode = ['up','down','up_fl','down_fl']
        # sanity check of forward kinematics
        for i,angles in enumerate(IK_angles.T):
            t44, _ = forward_kinematics(*angles, **DH_TABLE)
            print mode[i/3]
            print '\nERROR:'
            print round(n.abs(angles-joint_values))
            print '\nFK-ERROR:'
            print norm(T44 - t44)
            print '\nANGLE-ERROR:'
            print norm(angles - mat(joint_values))
            assert(norm(T44 - t44) < 1e-10)
            print '---'
            assert(norm(T44 - t44) < 1e-7)
        print "forward kinematics ok!"

        # list of global-robot-frames
        global_robot_frames = construct_robot_geometry(debug)
            
        # plot robot frames
        ax = fig.add_subplot(2,2,1)
        plot_robot_geometry(ax, global_robot_frames)

        T44_IK, debug_IK = forward_kinematics(*elbow_up[:,0], **DH_TABLE)
        global_IK_frames = construct_robot_geometry(debug_IK)

        plot_robot_geometry(ax, global_IK_frames,'g--')
##        plot_equal_perspective(ax,
##                                   [-0.5,0.5],
##                                   [-0.5,0.5],
##                      +             [0,0.5])
        # plot robot frames
        ax = fig.add_subplot(2,2,2)
        plot_robot_geometry(ax, global_robot_frames)

        T44_IK, debug_IK = forward_kinematics(*elbow_down[:,0], **DH_TABLE)
        global_IK_frames = construct_robot_geometry(debug_IK)

        plot_robot_geometry(ax, global_IK_frames,'r--')
##        plot_equal_perspective(ax,
##                                   [-0.5,0.5],
##                                   [-0.5,0.5],
##                                   [0,0.5])

        # plot robot frames
        ax = fig.add_subplot(2,2,3)
        plot_robot_geometry(ax, global_robot_frames)

        T44_IK, debug_IK = forward_kinematics(*elbow_up_fl[:,0], **DH_TABLE)
        global_IK_frames = construct_robot_geometry(debug_IK)

        plot_robot_geometry(ax, global_IK_frames,'g--')
##        plot_equal_perspective(ax,
##                                   [-0.5,0.5],
##                                   [-0.5,0.5],
##                                   [0,0.5])
        # plot robot frames
        ax = fig.add_subplot(2,2,4)
        plot_robot_geometry(ax, global_robot_frames)

        T44_IK, debug_IK = forward_kinematics(*elbow_down_fl[:,0], **DH_TABLE)
        global_IK_frames = construct_robot_geometry(debug_IK)

        plot_robot_geometry(ax, global_IK_frames,'r--')
##        plot_equal_perspective(ax,
##                                   [-0.5,0.5],
##                                   [-0.5,0.5],
##                                   [0,0.5])
        show()
        break
