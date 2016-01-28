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

sys.path.append('../int/djikstra/')
from graph import shortestPath as shortest_path


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

def calc_pair_norms(pi, pj):
    """
    Calculate all the norms for all the solutions
    between points pi and pj.
    """
    norms_i = []
    for si in pi:
        norm_ij = []
        for sj in pj:
            norm_ij.append( norm(si - sj) )
        norms_i.append( norm_ij )
    return norms_i
    
def map_point_norms(point_solutions):
    """
    Calculate all the norms for all the solutions
    between points pi and pj, for all points.
    """
    res = []
    num_solutions = len(point_solutions)
    for i in xrange(num_solutions-1):
        p_i = point_solutions[i]
        p_j = point_solutions[i+1]
        pair_norms_ij = calc_pair_norms(p_i, p_j)
        res.append(pair_norms_ij)
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



def __inverse_kinematics_sanity_check(T44_point, ik_angles):
    """
    sanity check of forward kinematics, just making sure
    that all solutions are valid.
    """
    for angles in ik_angles.T:
        t44, _ = forward_kinematics(*angles, **DH_TABLE)
        try:
            norm_value = norm(T44_point - t44)
            if not numpy.isnan(norm_value):
                assert(norm_value < 1e-10)
        except AssertionError:
            ##TODO: move this into utils for pretty printing matrices
            st = ''
            for row in T44_point:
                st += '\t'+str(row)+'\n'
            raise Exception('\n\tInverse-kinematics failed for point-frame: \n{0}\n\tNorm: {1}'.format(st, norm_value))
    return

def inverse_kinematics_joints(*joint_values, **DH_TABLE):
    # get forward kinematics i.e. last global robot-frame
    T44_point, geometry_info = forward_kinematics(*joint_values, **DH_TABLE)

    # get forward kinematics i.e. last global robot-frame
    ik_angles = inverse_kinematics_irb140(DH_TABLE, T44_point)

    # perform solution check of end-effector
    __inverse_kinematics_sanity_check(T44_point,ik_angles)

    return T44_point, geometry_info

def inverse_kinematics_point(*args, **DH_TABLE):
    """
    Input is a point-frame ( 4x4 matrix of type [[R,t],[0,1]] ),
    the function allows input on the forms:
    1:    (rot, tilt, skew, x,y,z)                (3x, 3x)
    2:    (rot, tilt, skew, t)                    (3x, 1x)
    3:    (angles, x,y,z)                         (1x, 3x)
    3.1:  (R, x,y,z)                              (1x, 3x)
    4:    (R, t), where R is list or numpy.array  (1x, 1x)
    """
    T44_point = homogenous_matrix(*args)

    # get inverse kinematics i.e. valid joint-value configurations
    ik_angles = inverse_kinematics_irb140(DH_TABLE, T44_point)

    # perform solution check of end-effector
    __inverse_kinematics_sanity_check(T44_point,ik_angles)

    return ik_angles

def inverse_kinematics_curve(trans_frames):
    # perform inverse kinematics over a curve and collect all solutions
    all_solutions = []
    for point_frame in trans_frames:
        angle_solutions = inverse_kinematics_irb140(DH_TABLE, point_frame)

        extra = [angle_solutions]
        for index in xrange(6):
            extra.append( generate_modulo_solutions(angle_solutions, index, 360.0))
            extra.append( generate_modulo_solutions(angle_solutions, index, -360.0))
            pass
        angle_solutions = merge_solutions(*extra)
        angle_solutions = filter_solutions(angle_solutions)
        all_solutions.append(angle_solutions.T)
    return mat(all_solutions)

def generate_solutions_graph(all_solutions):
    point_norms = map_point_norms(all_solutions)
    node_edges = {}
    for i in xrange( len(point_norms) ):
        norms_i = point_norms[i]
        map_edge_connections(i, norms_i, node_edges)

    #fix the ends that are not connected to anything
    glob_ends = ['p('+str(len(point_norms))+','+str(i)+')' for i in xrange(len(all_solutions[-1]))]
    for k in glob_ends:
        node_edges[k] = {}

    graph = node_edges
    return graph

def map_solution_paths(solution_graph, all_solutions):
    num_starts = len(all_solutions[0])
    num_ends = len(all_solutions[-1])
    solution_paths = []

    for s in xrange(num_starts):
        for e in xrange(num_ends):
            path = shortest_path(solution_graph,
                   'p(0,{0})'.format(str(s)),
                   'p({0},{1})'.format( str(len(all_solutions) - 1), str(e)))
            S = get_solutions_from_node_ids(all_solutions, *path)
            solution_paths.append( mat(S) )
    return solution_paths

def find_inverse_kinematics_paths_from_curve(trans_frames):
    start = time.time()
    all_solutions = inverse_kinematics_curve(trans_frames)
    print 'inverse-kinematics, curve: {0}s'.format(time.time() - start)

    start = time.time()
    solution_graph = generate_solutions_graph(all_solutions)
    print 'solution graph: {0}s'.format(time.time() - start)
    
    
    start = time.time()
    solution_paths = map_solution_paths(solution_graph, all_solutions)
    print 'Djikstra: {0}s'.format(time.time() - start)
    

    start = time.time()
    all_solution_distances = apply_along_axis(apply_along_axis(solution_paths, func=diff, axis=1),func=norm, axis=2)
    result = {
        'solutions_per_point' : all_solutions,
        'solution_graph' : solution_graph,
        'solution_paths': solution_paths,
        'solution_path_nodes_differences' : all_solution_distances
        }
    print 'Collect result: {0}s'.format(time.time() - start)
    return result

def plot_robot_from_angles(plot, *args):
    s = forward_kinematics(*args, **DH_TABLE)
    plot.draw_robot(s['robot_geometry_global'])
    return

def plot_path(plot, paths, index):
    for i in xrange(6):
        #plot.add_subplot(6,1,i+1, title='j'+str(i+1))
        plot.plot(paths[index,:,i])
        #plt.legend(['j{0}'.format(i+1)])
    #fig.show()
    

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

        # inverse knematics over curve
        start = time.time()
        result = find_inverse_kinematics_paths_from_curve(transf_frames)
        stop = time.time()
        
        # results
        print 'Time: {0}'.format(stop - start)
        print '\n'

        print result.keys()
        print '\n'
        
        print 'total paths available:' + str( len(result['solution_paths']) )
        count = 0

        valid = []
        for i,k in enumerate( result['solution_path_nodes_differences'] ):
            if n.max(abs(k)) < 20.0:
                count = count + 1
##                print 'max-err: ' + str(n.max(abs(k)))
##                print 'index: ' + str(i)
                valid.append( i )
        print 'valid paths: ' + str(count)

##        for solution_distance in result['solution_path_nodes_differences']:
##            plot(solution_distance)
##        show()
        paths = mat(result['solution_paths'])

        k = valid[-1]
        path = result['solution_paths'][k]
        path_diff = result['solution_path_nodes_differences'][k]
        plot(path_diff)
        show()

        # plotting
        plot = StPlot()
        plot_path(plot, paths, k)
        plot.draw_robot(robot_frames)
##        for i in xrange(0,len(path),len(path)/5):
##            plot_robot_from_angles(plot, *path[i])
        plot.draw_trajectory(transf_frames)
        plot.show()
