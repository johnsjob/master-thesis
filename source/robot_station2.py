import random

from helperfunctions_plot import *
from plane_relative import *
from denavit_hartenberg140 import *

import itertools as it

def apply_along_axis(M, func=n.diff, axis=1):
    return n.apply_along_axis(func, axis, arr=M)

def get_closest_solutions_pair(s0, s1):
    compares = []
    for i, s0i in enumerate( s0 ):
        for j, s1j in enumerate( s1 ):
            #print "%d, %d, %d" % (norm(s0i-s1j), i , j)
            compares.append( [norm(s0i-s1j), i , j])
    comp = mat( compares )

    ret = []
    try:
        wh = n.argwhere(comp == comp.min(0)[0])
        for k in comp[wh[:,0],:]:
            i,j = k[1:3]

            #norm_value = k[0]
            pair = [s0[i], s1[j]]
            ret.append( pair )
#            import pdb; pdb.set_trace()
    except:
        pass
##    if len(ret) > 1:
##        import pdb; pdb.set_trace()
    return ret

def add_solutions(solutions, solution_value, index=5):
    for s in solutions.T:
        tmp1 = s.copy()
        tmp2 = s.copy()
        old_val = s[index]
        tmp1[index] = old_val + solution_value
        yield tmp1
        tmp2[index] = old_val - solution_value
        yield tmp2

def traverse_solutions(*args):
    for solutions in args:
        for s in solutions.T:
            yield s

def make_array(list_of):
    return mat(list_of).T

		
if __name__ == '__main__':
    for count in n.linspace(-180,180,10):
        ax, fig = init_plot()
        fig.clear()
        j1 =  count
        j2 =  rand_range(-90, 110)
        j3 =  rand_range(-230, 50)
        j4 =  rand_range(-200, 200)
        j5 =  rand_range(-115, 115)
        j6 =  rand_range(-400, 400)

        joint_values = j1,j2,j3,j4,j5,j6

        # get forward kinematics i.e. last global robot-frame
        T44, debug = forward_kinematics(*joint_values, **DH_TABLE)
        IK_angles = inverse_kinematics_irb140(DH_TABLE, T44)

        # sanity check of forward kinematics
        for angles in IK_angles.T:
            t44, _ = forward_kinematics(*joint_values, **DH_TABLE)
            assert(norm(T44 - t44) < 1e-7)

        # the base frame which defines the world-coordinate system
        plane0 = define_plane_from_angles([0,0,0],0, 0, 0)

        # list of global-robot-frames
        global_robot_frames = matmul_series(*debug)
        global_robot_frames.insert(0, plane0)
        global_robot_frames = mat( global_robot_frames )
        global_robot_points = global_robot_frames[:, :3, 3]

        # generate a curve in the last global robot-frame
        num_p = 50
        point_matrix = generate_symmetric_curve(num_points=num_p)
        point_matrix_tf = get_transformed_points(T44, point_matrix)

    ######
        ax = fig.add_subplot(1,2,1, projection='3d')
        for robot_frame in global_robot_frames:
            plot_plane(ax,robot_frame, '--',scale_factor=0.1)

        ax.scatter(point_matrix_tf[:,0],
                   point_matrix_tf[:,1],
                   point_matrix_tf[:,2])
        
        ax.plot(global_robot_points[:,0],
                global_robot_points[:,1],
                global_robot_points[:,2], 'k',linewidth=2)

        plot_equal_perspective(ax,
                               [-0.5,0.5],
                               [-0.5,0.5],
                               [0,1])
        #show()
    ######
        plane = global_robot_frames[0]
        global_plane_curve = point_matrix_tf

        lost_p = 0
        all_solutions = []
        for point in global_plane_curve:
            FK_p = homogenous_matrix(plane[:3,:3],
                                    point[:3])
            angle_solutions = inverse_kinematics_irb140(DH_TABLE, FK_p)
            angle_solutions = filter_solutions( angle_solutions )
            angle_solutions = angle_solutions.T
            if n.sum( angle_solutions.shape) == 0.0:
                lost_p += 1
                continue
            print angle_solutions.shape
            all_solutions.append(angle_solutions)

        chosen_solutions = []
        for k in xrange(1, len(all_solutions)):
            if k == 1:
                o = all_solutions[k-1]
            else:
                o = chosen_solutions[-1]
            pairs = mat(get_closest_solutions_pair(o, all_solutions[k]))

            if k==1:
                chosen_solutions.append(pairs[0,0,:].reshape(1,6))
                chosen_solutions.append(pairs[0,1,:].reshape(1,6))
            else:
                chosen_solutions.append(pairs[0,1,:].reshape(1,6))

        chosen_solutions = mat(chosen_solutions).reshape(num_p - lost_p,6)
        diff_solutions = apply_along_axis(chosen_solutions, func=n.diff, axis=0)
        max_err_solutions = n.max(n.abs(diff_solutions), axis=1)

        ax = fig.add_subplot(1,2,2)
        plot(max_err_solutions)
    show()
