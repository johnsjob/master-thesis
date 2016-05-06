import random

from helperfunctions_plot import *
from pylab import axhline
from plane_relative import *
from denavit_hartenberg140 import forward_kinematics,\
     DH_TABLE, inverse_kinematics_irb140, generate_modulo_solutions,\
     merge_solutions, filter_solutions

import itertools as it

def apply_along_axis(M, func, axis=1):
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
        
if __name__ == '__main__':
    for count in xrange(1000):
        ax, fig = init_plot()
        fig.clear()
        j1 =  rand_range(-120,120)
        j2 =  rand_range(-90, 110)
        j3 =  rand_range(-230, 50)
        j4 =  rand_range(-200, 200)
        j5 =  rand_range(-115, 115)
        j6 =  rand_range(-400, 400)

        j1 =  10
        j2 =  20
        j3 =  30
        j4 =  40
        j5 =  50
        j6 =  60

        joint_values = j1,j2,j3,j4,j5,j6

        # get forward kinematics i.e. last global robot-frame
        robot_info = forward_kinematics(*joint_values, **DH_TABLE)
        T44 = robot_info['tcp']
        IK_angles = inverse_kinematics_irb140(DH_TABLE, T44)

        # sanity check of forward kinematics
        for angles in IK_angles.T:
            info = forward_kinematics(*joint_values, **DH_TABLE)
            t44 = info['tcp']
            assert(norm(T44 - t44) < 1e-7)

        # list of global-robot-frames
        global_robot_frames = construct_robot_geometry(info['robot_geometry_local'])
            
        # generate a curve in the last global robot-frame
        num_p = 50
        point_matrix = generate_symmetric_curve(num_points=num_p, ampl_factor=0.20)
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

            extra = [angle_solutions]
            for index in xrange(3,6):
                extra.append( generate_modulo_solutions(angle_solutions, index, 360.0))
                extra.append( generate_modulo_solutions(angle_solutions, index, -360.0))
            angle_solutions = merge_solutions(*extra)
            angle_solutions = filter_solutions(angle_solutions)

            all_solutions.append(angle_solutions.T)
        all_solutions = mat(all_solutions)
        print all_solutions.shape

        #check so that all chosen solutions are within angle-ranges
        try:
            chosen_solutions = extract_closest_solutions(all_solutions, norm)
        except:
            continue
        all_solutions_valid = mat(filter_solutions(chosen_solutions.T).shape) - mat(chosen_solutions.T.shape)
        all_solutions_valid = n.sum(all_solutions_valid) == 0.0
        assert(all_solutions_valid == True)
        print 'All solutions within valid ranges!'

        max_norm = lambda x,**kwargs: norm(x,ord=inf,**kwargs)
        diff_solutions = apply_along_axis(chosen_solutions, func=n.diff, axis=0)
        solution_distance = apply_along_axis(apply_along_axis(chosen_solutions, func=diff, axis=0),func=norm, axis=1)
        solution_distance_max = apply_along_axis(apply_along_axis(chosen_solutions, func=diff, axis=0),func=max_norm, axis=1)

##        stop_running = False
##        if not (n.max(solution_distance) < 20.0):
##            print 'too large deviation: '+str(n.max(solution_distance))
##            continue
        print apply_along_axis(n.abs(apply_along_axis(chosen_solutions, func=diff, axis=0)),func=n.max, axis=0)

        ax = fig.add_subplot(1,2,2)
        plot(solution_distance)
        plot(solution_distance_max)
        show()
        fig.clear()
        
        ax, fig = init_plot()
        fig.clear()

        joint_ranges = [[-180, 180],
                         [-90, 110],
                         [-230, 50],
                         [-200, 200],
                         [-115, 115],
                         [-400, 400]]

        for k in xrange(6):
            ax = fig.add_subplot(3,2,k+1)
            plot(chosen_solutions[:,k])
            axhline(joint_ranges[k][0])
            axhline(joint_ranges[k][1])
            
            legend(['j'+str(k+1)])
        show()
        break
