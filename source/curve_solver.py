import random

from helperfunctions_plot import *
from pylab import axhline
from plane_relative import *
from denavit_hartenberg140 import *

import itertools as it

import sys
sys.path.append('../int/djikstra/')


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

def construct_robot_geometry(fk_debug_info):
        plane0 = define_plane_from_angles([0,0,0],0, 0, 0)
        global_robot_frames = matmul_series(*fk_debug_info)
        global_robot_frames.insert(0, plane0)
        global_robot_frames = mat( global_robot_frames )
        return global_robot_frames

def plot_robot(ax, color='k', *joint_values):
    T44, debug = forward_kinematics(*joint_values, **DH_TABLE)
    robot_frames = construct_robot_geometry(debug)
    plot_robot_geometry(ax, robot_frames, color)

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

def my_norm(x, **kwargs):
    max_num = len(x)
    factors = mat([max_num-k for k in xrange(max_num)])
    return norm(factors*x)

def calc_pair_norms(p0, p1):
    res = []
    for s0 in p0:
        tmp = []
        for s1 in p1:
            N = norm(s0 - s1)
###            print N
            tmp.append(N)
        res.append(tmp)
###    print ''
    return res
    
def map_norms(solutions):
    res = []
    num_solutions = len(solutions)
    for i in xrange(num_solutions-1):
        p_i = solutions[i]
        p_j = solutions[i+1]
        pair_norms = calc_pair_norms(p_i, p_j)
        res.append(pair_norms)
##        break
    return res

def map_edge_connections(i,res_i, dict_res):
    res_i = mat(res_i)
    num_edges_j, num_edges_k = res_i.shape
    for j in xrange(num_edges_j):
        for k in xrange(num_edges_k):
            glob = 'p('+str(i)+','+str(j)+')'
            loc = 'p('+str(i+1)+','+str(k)+')'
            if not dict_res.has_key(glob):
                dict_res[glob] = {}
            dict_res[glob][loc] = res_i[j,k]
        dict_res[glob] = dict(dict_res[glob])
    #import pdb; pdb.set_trace()

def extract_solution_from_node_id(str_node_id, all_solutions):
    glob, loc = str_node_id[2:-1].split(',')
    glob = int(glob)
    loc = int(loc)
    return all_solutions[glob][loc]

def get_solutions_from_node_ids(all_solutions, *node_ids):
    res = []
    for node_id in node_ids:
        res.append( extract_solution_from_node_id(node_id, all_solutions) )
    return res


#### code saved for test-case
####        all_solutions = [
####            [[0],[0],[0]],
####            [[2],[1]],
####            [[3],[4],[5]]
####            ]
####        all_solutions = mat(all_solutions)
####        for k in xrange(len(all_solutions)):
####            all_solutions[k] = mat(all_solutions[k])


def inverse_kinematics_joints(*joint_values, **DH_TABLE):
    # get forward kinematics i.e. last global robot-frame
    T44, geometry_info = forward_kinematics(*joint_values, **DH_TABLE)
    ik_angles = inverse_kinematics_irb140(DH_TABLE, T44)

    # sanity check of forward kinematics
    for angles in ik_angles.T:
        t44, _ = forward_kinematics(*angles, **DH_TABLE)
        assert(norm(T44 - t44) < 1e-7)
    return T44, geometry_info

def inverse_kinematics_point(*args, **DH_TABLE):
    if len(args) == 1:
        input_type = type(args[0])
        if not  input_type in [list, numpy.array]:
            except ArithmeticError('one argument, assuming matrix - only list or numpy.array accepted for the T44 matrix.')
        elif input_type == list:
            T44_point = mat(list)
        
    import pdb; pdb.set_trace()
    # get forward kinematics i.e. last global robot-frame
    ik_angles = inverse_kinematics_irb140(DH_TABLE, T44_point)

    # sanity check of forward kinematics
    for angles in ik_angles.T:
        t44, _ = forward_kinematics(*angles, **DH_TABLE)
        assert(norm(T44_point - t44) < 1e-7)
    return ik_angles


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

        j1 =  180
        j2 =  0
        j3 =  0
        j4 =  0
        j5 =  0
        j6 =  0

        joint_values = j1,j2,j3,j4,j5,j6
        break
 

            
####        # generate a curve in the last global robot-frame
####        num_p = 50
####        point_matrix = generate_symmetric_curve(num_points=num_p, ampl_factor=0.30)
####        point_matrix_tf = get_transformed_points(T44, point_matrix)
####
####        # plot robot frames
####        ax = fig.add_subplot(1,2,1, projection='3d')
####        plot_robot_geometry(ax, global_robot_frames,'b--')
####        plot_curve(ax, point_matrix_tf)
####        plot_equal_perspective(ax,
####                               [-0.5,0.5],
####                               [-0.5,0.5],
####                               [0,1])
####        #show()
####
####        # rename some variables for convenience
####        plane = global_robot_frames[-1]
####        global_plane_curve = point_matrix_tf
####
####        # perform inverse kinematics over a curve and collect all solutions
####        all_solutions = []
####        for point in global_plane_curve:
####            fk_p = homogenous_matrix(plane[:3,:3],
####                                     point[:3])
####            angle_solutions = inverse_kinematics_irb140(DH_TABLE, fk_p)
####            extra = [angle_solutions]
####            for index in xrange(3,6):
####                extra.append( generate_modulo_solutions(angle_solutions, index, 360.0))
####                extra.append( generate_modulo_solutions(angle_solutions, index, -360.0))
####            angle_solutions = merge_solutions(*extra)
####            angle_solutions = filter_solutions(angle_solutions)
####            all_solutions.append(angle_solutions.T)
####        all_solutions = mat(all_solutions)
####
####
####        all_solutions = mat(all_solutions)
####        print all_solutions.shape
####        print mat(all_solutions[0]).shape
####        res = map_norms(all_solutions)
####        print '#1'
####        d = {}
####        for i in xrange(len(res)):
####            res_i = res[i]
####            map_edge_connections(i, res_i, d)
#######        for k in sorted(d.keys()):
#######            print k+':'
#######            for l in sorted(d[k].keys()):
#######                print '\t'+l+' = '+str(d[k][l])
####        print '#2'
####        #fix the ends that are not connected to anything
####        glob_ends = ['p('+str(len(res))+','+str(i)+')' for i in xrange(len(all_solutions[-1]))]
####        for k in glob_ends:
####            d[k] = {}
####        graph = d
####        #graph['p(2,2)'] = {}
####        from graph import shortestPath as sp
####        num_starts = len(all_solutions[0])
####        num_ends = len(all_solutions[-1])
####        chosen_solutions = []
####        print '#3'
####        import time
####        _start = time.time()
####        for s in xrange(num_starts):
####            print s
####            for e in xrange(num_ends):
####                #R = sp(graph, 'p(0,'+str(s)+')','p(49,'+str(e)+')')
####                R = sp(graph, 'p(0,'+str(s)+')','p('+str(len(all_solutions) - 1)+','+str(e)+')')
####                S = get_solutions_from_node_ids(all_solutions, *R)
####                S = mat(S)
####                chosen_solutions.append(S)
####        _stop = time.time()
####        print 'time: ' + str(_stop-_start)
####        print '#4'
####        all_solution_distances = apply_along_axis(apply_along_axis(chosen_solutions, func=diff, axis=1),func=norm, axis=2)
####        ax = fig.add_subplot(1,2,2)
####        for solution_distance in all_solution_distances:
####            plot(solution_distance)
####        print 'total paths available:' + str(len(chosen_solutions))
####        count = 0
####        for i,k in enumerate(all_solution_distances):
####            if n.max(abs(k)) < 20:
####                count = count + 1
####                print 'max-err: ' + str(n.max(abs(k)))
####                print 'index: ' + str(i)
####        print 'valid paths: ' + str(count)
##########        ax = fig.add_subplot(1,3,3, projection='3d')
##########        for k in chosen_solutions:
##########            for l in k:
##########                #print len(l)
##########                t44, deb = forward_kinematics(*l, **DH_TABLE)
##########                frames = construct_robot_geometry(deb)
##########                plot_robot_geometry(ax, frames)
##########                assert(norm(T44[:3,:3] - t44[:3,:3]) < 1e-7)
##########            show()
####        show()
####        print mat(chosen_solutions)[175,10,:]
####        plot(mat(chosen_solutions)[175,:,:3]); show()
####        show()
####        break
