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
            #N = abs(s0[0]-s1[0])
            print N
            tmp.append(N)
        res.append(tmp)
    print ''
##        break
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

        j1 =  0
        j2 =  90
        j3 =  0
        j4 =  0
        j5 =  0
        j6 =  0

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
####        all_solutions = [
####            [[0],[0],[0]],
####            [[2],[1]],
####            [[3],[4],[5]]
####            ]
####        all_solutions = mat(all_solutions)
####        for k in xrange(len(all_solutions)):
####            all_solutions[k] = mat(all_solutions[k])
########        #1/0
        all_solutions = mat(all_solutions)
        print all_solutions.shape
        print mat(all_solutions[0]).shape
        res = map_norms(all_solutions)

        d = {}
        for i in xrange(len(res)):
            res_i = res[i]
            map_edge_connections(i, res_i, d)
        for k in sorted(d.keys()):
            print k+':'
            for l in sorted(d[k].keys()):
                print '\t'+l+' = '+str(d[k][l])
        #fix the ends that are not connected to anything
        glob_ends = ['p('+str(len(res))+','+str(i)+')' for i in xrange(len(all_solutions[-1]))]
        for k in glob_ends:
            d[k] = {}
        graph = d
        #graph['p(2,2)'] = {}
        from graph import shortestPath as sp
        R = sp(graph, 'p(0,0)','p(49,0)')
        print R
        S = get_solutions_from_node_ids(all_solutions, *R)
        S = mat(S)
        print n.round(n.diff(S,axis=0))
        #import pdb; pdb.set_trace()
        break
####        #check so that all chosen solutions are within angle-ranges
####        try:
####            chosen_solutions = extract_closest_solutions(all_solutions, norm)
####        except:
####            continue
####        all_solutions_valid = mat(filter_solutions(chosen_solutions.T).shape) - mat(chosen_solutions.T.shape)
####        all_solutions_valid = n.sum(all_solutions_valid) == 0.0
####        assert(all_solutions_valid == True)
####        print 'All solutions within valid ranges!'
####
####        max_norm = lambda x,**kwargs: norm(x,ord=inf,**kwargs)
####        diff_solutions = apply_along_axis(chosen_solutions, func=n.diff, axis=0)
####        solution_distance = apply_along_axis(apply_along_axis(chosen_solutions, func=diff, axis=0),func=norm, axis=1)
####        solution_distance_max = apply_along_axis(apply_along_axis(chosen_solutions, func=diff, axis=0),func=max_norm, axis=1)
####
######        stop_running = False
######        if not (n.max(solution_distance) < 20.0):
######            print 'too large deviation: '+str(n.max(solution_distance))
######            continue
####        print apply_along_axis(n.abs(apply_along_axis(chosen_solutions, func=diff, axis=0)),func=n.max, axis=0)
####
####        ax = fig.add_subplot(1,2,2)
####        plot(solution_distance)
####        plot(solution_distance_max)
####        show()
####        fig.clear()
####        
####        ax, fig = init_plot()
####        fig.clear()
####
####        joint_ranges = [[-180, 180],
####                         [-90, 110],
####                         [-230, 50],
####                         [-200, 200],
####                         [-115, 115],
####                         [-400, 400]]
####
####        for k in xrange(6):
####            ax = fig.add_subplot(3,2,k+1)
####            plot(chosen_solutions[:,k])
####            axhline(joint_ranges[k][0])
####            axhline(joint_ranges[k][1])
####            
####            legend(['j'+str(k+1)])
####        show()
####        break
