import random

from helperfunctions_plot import *
from plane_relative import *
from denavit_hartenberg140 import *

import itertools as it

def apply_along_axis(M, func=n.diff, axis=1):
    return n.apply_along_axis(func, axis, arr=M)

def plot_robot_geometry(ax, global_robot_frames):
        for robot_frame in global_robot_frames:
            plot_plane(ax,robot_frame, '--',scale_factor=0.1)
        
        ax.plot(global_robot_frames[:,0,3],
                global_robot_frames[:,1,3],
                global_robot_frames[:,2,3], 'k',linewidth=2)

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

def get_closest_solutions_pair(s0, s1):
    data = []
    #import pdb; pdb.set_trace()
    s0 = filter_solutions(s0)
    s1 = filter_solutions(s1)
    for i, s0i in enumerate(s0.T):
        for j, s1j in enumerate(s1.T):
            print norm(s0i-s1j, ord=inf)
            data.append([norm(s0i - s1j, ord = inf), i, j])
    data = mat(data)
    print '+++'

    ret = []
    solution_col_row_pairs = n.argwhere(data == data.min(axis = 0)[0])
    solution_indices = solution_col_row_pairs[:,0]
    print len(data[solution_indices])
    for solution_data in data[solution_indices]:
        print solution_data
        norm_value, i, j = solution_data
        pair = mat([s0.T[i], s1.T[j]])
        print '###'
        print norm(pair[0]-pair[1], ord=inf)
        print
        return pair


if __name__ == '__main__':
    for count in n.linspace(-180,180,10):
        ax, fig = init_plot()
        fig.clear()
        j1 =  count
        j2 =  0#rand_range(-90, 110)
        j3 =  0#rand_range(-230, 50)
        j4 =  0#rand_range(-200, 200)
        j5 =  0#rand_range(-115, 115)
        j6 =  0#rand_range(-400, 400)

        joint_values = j1,j2,j3,j4,j5,j6

        # get forward kinematics i.e. last global robot-frame
        T44, debug = forward_kinematics(*joint_values, **DH_TABLE)
        IK_angles = inverse_kinematics_irb140(DH_TABLE, T44)

        # sanity check of forward kinematics
        for angles in IK_angles.T:
            t44, _ = forward_kinematics(*joint_values, **DH_TABLE)
            assert(norm(T44 - t44) < 1e-7)

        # list of global-robot-frames
        global_robot_frames = construct_robot_geometry(debug)
            
        # generate a curve in the last global robot-frame
        num_p = 50
        point_matrix = generate_symmetric_curve(num_points=num_p)
        point_matrix_tf = get_transformed_points(T44, point_matrix)

        # plot robot frames
        ax = fig.add_subplot(1,2,1, projection='3d')
        plot_robot_geometry(ax, global_robot_frames)
        plot_curve(ax, point_matrix_tf)
        plot_equal_perspective(ax,
                               [-0.5,0.5],
                               [-0.5,0.5],
                               [0,1])
        #show()

        # rename some variables for convenience
        plane = global_robot_frames[-1]
        global_plane_curve = point_matrix_tf

        # perform inverse kinematics over a curve and collect all solutions
        all_solutions = []
        for point in global_plane_curve:
            fk_p = homogenous_matrix(plane[:3,:3],
                                     point[:3])
            angle_solutions = inverse_kinematics_irb140(DH_TABLE, fk_p)
            #angle_solutions = filter_solutions(angle_solutions)
            #angle_solutions = angle_solutions.T
            for k in angle_solutions.T:
                T44, debug = forward_kinematics(*k, **DH_TABLE)
                # sanity check
                err = norm(fk_p - T44)
                #assert( err < 1e-10)
            print angle_solutions.shape
            all_solutions.append(angle_solutions)
        all_solutions = mat(all_solutions)
        print all_solutions.shape

        # obtain closest solutions
        chosen_solutions = []
        num_solutions = len(all_solutions)
        for k in xrange(num_solutions-1):
            if k == 0:
                o = all_solutions[k]
            else:
                o = chosen_solutions[-1]

            pair = get_closest_solutions_pair(o, all_solutions[k+1])

            if k == 0:
                chosen_solutions.append(mat([pair[0]]).T)
                chosen_solutions.append(mat([pair[1]]).T)
            else:
                chosen_solutions.append(mat([pair[1]]).T)

        chosen_solutions = mat(chosen_solutions).reshape(num_solutions,6)
        diff_solutions = apply_along_axis(chosen_solutions, func=n.diff, axis=0)
        max_err_solutions = n.max(n.abs(diff_solutions), axis=1)
        max_err_solutions = apply_along_axis(apply_along_axis(chosen_solutions, func=diff, axis=0),func=norm, axis=1)

        ax = fig.add_subplot(1,2,2)
        plot(max_err_solutions)
        show()
        break
