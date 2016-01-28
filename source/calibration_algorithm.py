# -*- coding: cp1252 -*-
from __future__ import division
#----------------------------------------#
import numpy
import time

from numpy.linalg import solve, det, inv, cond
from numpy import array as mat, log10
#----------------------------------------
# custom imports
from helperfunctions_math import *
from helperfunctions_plot import *
import plane_relative as plane_tools
#----------------------------------------#
# num_points = 12 absolute minimum, actually 12+1
num_points = 160
#========================================#
# placing paper origin
o = mat([1000*rand(), 1000*rand(), 1000*rand()])

### defining paper orientation
r,t,s = -rand_range(-180,180), rand_range(-180,180), rand_range(-180, 180)

# define the paper-orientation in global (robot-system) directions
# however, the paper uses "euclidean" orientation ('local' = euclidean)
plane = plane_tools.define_plane_from_angles(o, r, t, s, 'local')
#################################################
# define delta vector which we want to find (in tool-space / local space)
L = 100
local_delta_vector = mat([1,2,3])
local_delta_vector = (local_delta_vector / norm(local_delta_vector))*L #length L
# Orientation of the tool in tool (local) coordinate system 
local_tool_orientation = rotation_matrix_rot_tilt_skew(-10, 20, 30)
#----------------------------------------
def merge_dicts(*list_of_dicts):
    # init
    ret = {}
    keys = []

    # get all unique keys
    for d in list_of_dicts:
        keys += d.keys()
    keys = set().union(keys)

    # for all keys ...
    for k in keys:
        # prepare a k:th-list if none exists
        if not ret.has_key(k):
            ret[k] = []
        # for all dicts ...
        for d in list_of_dicts:
            # if dict has key ...
            if d.has_key(k):
                # check so that the key is not an empty list ...
                empty = False
                try:
                    empty = len(d[k]) == 0
                except:
                    # not a list/array-type, equivalent to non-empty list
                    pass
                # append item or non-empty list
                if not empty:
                    ret[k].append( d[k] )
    # for all keys ...
    for k in keys:
        # if we only got one item for this key from al the dicts ...
        if len(ret[k]) == 1:
            # un-list it
            ret[k] = ret[k][0]
        # remove empy lists if any manage to get here
        elif len(ret[k]) == 0:
            del ret[k]
        else:
            # turn remaining lists into numpy-arrays
            ret[k] = mat(ret[k])
    return ret
#----------------------------------------
def rad_to_ang(v):
    return v*180/pi
#----------------------------------------
def problem_formulation(dx, dy, dR):
    r11,r12,r13,r21,r22,r23,r31,r32,r33 = (-dR).reshape(9)
    S1 = [dx, dy, r11, 0,   0, r12, 0, 0,   r13]
    S2 = [0,   0, r21, dx, dy, r22, 0, 0,   r23]
    S3 = [0,   0, r31, 0,   0, r32, dx, dy, r33]
    row_value = 3
    col_value = 9
    return mat([S1, S2, S3]),  row_value, col_value
#----------------------------------------
def vec_diff(v1, v2):
    err = norm(v1 - v2)
    norm_err = abs(norm(v1) - norm(v2))
    angle_err = rad_to_ang(acos( (v1/norm(v1)).dot((v2/norm(v2))) ))
    return err, norm_err, angle_err
#----------------------------------------
def generate_random_Anoto_Point(L):
    px = L*rand()-L/2.0
    py = L*rand()-L/2.0
    return px, py
#----------------------------------------
def solve_tool0_tip_alt(array_forward_kinematics_T44,
                        array_anoto2D,
                        array_lhs_sys_eq = None):
    try:
        num_points,m,n = array_forward_kinematics_T44.shape
    except Exception as e:
        print 'solve_tool0_tip:\n\tWrong shape or type for input parameter: array_forward_kinematics_T44'
    try:
        m,n = array_anoto2D.shape
    except Exception as e:
        print 'solve_tool0_tip:\n\tWrong shape or type for input parameter: array_anoto2D'

    l_xtcp = array_forward_kinematics_T44[:, 0:3, 3]
    l_R = array_forward_kinematics_T44[:, 0:3, 0:3]

    dxtcp = diff(l_xtcp, axis=0)
    dR = diff(l_R, axis=0)
    danoto2D = diff(array_anoto2D, axis=0)

    lhs = []
    rhs = []
    l_cond = []
    l_err = []

    for i in xrange(0, num_points-1): #one less after forward-differences....
        A, row_value, col_value = array_lhs_sys_eq(danoto2D[i,0], danoto2D[i,1], dR[i])
        b = dxtcp[i]
        lhs.append(A)
        rhs.append(b)
    lhs = mat(lhs).reshape(((num_points-1) * row_value, col_value))

    #shape the rhs depending on shape-info from lhs
    if row_value != 1:
        rhs = mat(rhs).reshape((num_points-1) *  row_value)
    else:
        rhs = mat(rhs)
        
    L = lhs.T.dot(lhs)
    R = lhs.T.dot(rhs)

    X, Y, D = solve(L, R).reshape(3,3).T
    Z = cross(X, Y)

    X = X / numpy.linalg.norm(X)
    Y = Y / numpy.linalg.norm(Y)
    Z = Z / numpy.linalg.norm(Z)
    result = mat([X,Y,Z,D]).T
    condition = cond(L)

    return result, condition
#----------------------------------------
def generate_Xflange_orientation(plane,rot, tilt, skew):
    """
        Generate Xtcp-orientation in world coordinates, using Anoto-paper
        orientation formulation
        Planes are homoenous matrices, if we want the orientation
        we need the 0:3,0:3 submatrix.
    """
    return plane_tools.define_plane_relative_from_angles(plane, (0,0,0),
                                                         rot, tilt, skew,'global')[:3,:3]
#----------------------------------------
def setup_geometry(current_plane, point_spread, num_points):
    global local_delta_vector, local_tool_orientation

    geometry_info = {'plane':current_plane}
    geometry_info['local_tool_orientation'] = local_tool_orientation
    geometry_info['local_delta_vector'] = local_delta_vector

    geometry_info['correct_solution_geometry'] = mat(list(geometry_info['plane'][:3,:3].T.flatten()) +
               list(geometry_info['local_delta_vector'])).reshape(4,3).T


    #generating points and "forward-kinematics"
    collected_data = []
    for k in xrange(0,num_points):
        info = {}
        info['angles'] = \
        {
            'rot':  rand_range(-180,180),
            'tilt': rand_range(-60, 60),
            'skew': rand_range(-180,180) 
        }
        # Xtcp (flange) orientation in global space, generated relative to the paper plane
        info['Xflange_orientation_relative_to_paper_plane'] = \
                                                generate_Xflange_orientation(geometry_info['plane'],**info['angles'])
        #generate pen-tip position in Anoto2d in mm
        px,py = generate_random_Anoto_Point(point_spread)
        info['pentip_2d'] = [px,py]

        # generate global Xtcp position in mm
        info['Xtcp0'] = (plane_tools.get_plane_point(geometry_info['plane'], px, py)[:3] - \
                 matmul( info['Xflange_orientation_relative_to_paper_plane'], local_delta_vector[:3]) )
        # ^OK

        # generate relative-tool-orientation in world coordinates
        info['global_tool_orientation'] = matmul( info['Xflange_orientation_relative_to_paper_plane'],
                                                  local_tool_orientation )
        # ^OK

        info['forward_kinematics'] = homogenous_matrix( info['Xflange_orientation_relative_to_paper_plane'],
                                                        info['Xtcp0'] )
        # ^OK
        collected_data.append(info)
    geometry_info['data'] = merge_dicts(*collected_data)
    geometry_info['data']['angles'] = merge_dicts(*geometry_info['data']['angles'])

    return geometry_info
#----------------------------------------
def find_solution_pen_tip(geometry_info, included_solutions_from_start = -1):
    start_time = time.clock()
    
    result, cond_num = solve_tool0_tip_alt(geometry_info['data']['forward_kinematics'][:included_solutions_from_start],
                                      geometry_info['data']['pentip_2d'][:included_solutions_from_start],
                                      problem_formulation)
    return result, cond_num

def find_solution_pen_ori(geometry_info, included_solutions_from_start = -1):
    # solve for orientation s which should be same as local_tool_orientation
    l,m,n = geometry_info['data']['Xflange_orientation_relative_to_paper_plane'].shape
    flange_orientation_reshaped = geometry_info['data']['Xflange_orientation_relative_to_paper_plane'].reshape(l*m,n)

    lhs = matmul( flange_orientation_reshaped.T,
                  flange_orientation_reshaped)

    l,m,n = geometry_info['data']['global_tool_orientation'].shape
    rhs = matmul( flange_orientation_reshaped.T,
                  geometry_info['data']['global_tool_orientation'].reshape(l*m,n))

    solved_tool_orientation = linalg.solve(lhs, rhs)

    #normalize result
    solved_tool_orientation[:,0] = solved_tool_orientation[:,0] / norm(solved_tool_orientation[:,0])
    solved_tool_orientation[:,1] = solved_tool_orientation[:,1] / norm(solved_tool_orientation[:,1])
    solved_tool_orientation[:,2] = solved_tool_orientation[:,2] / norm(solved_tool_orientation[:,2])
    return solved_tool_orientation, cond(lhs)
#----------------------------------------
def perform_solution_run(geometry_info):
    interval = range(3,num_points)
    list_of_solving = []
    for k in interval:
        solve_info = {}
        solve_info['result'], solve_info['solved_tool_orientation'], solve_info['cond_num'], solve_info['time_spent'] = \
                                find_solution(geometry_info, k)                    
        solve_info['err'] = norm(geometry_info['correct_solution_geometry'] - solve_info['result']) + norm(geometry_info['local_tool_orientation'] - solve_info['solved_tool_orientation'])
        list_of_solving.append(solve_info)
    solving_data = merge_dicts(*list_of_solving)
    solving_data['interval'] = interval
    print 'solution error = ' + str( solving_data['err'][-1] )
    return solving_data
#----------------------------------------
def make_plots(solving_data):
    logcond = log10( solving_data['cond_num'] )
    logerr  = log10( solving_data['err'] )
    plot(solving_data['interval'], logcond, label='Condition number',       linewidth=2);
    plot(solving_data['interval'], logerr,  label='Error (frobenious norm)',linewidth=2);
    plot(solving_data['interval'], log10( solving_data['time_spent'] ),'r', label='time spent',linewidth=2)
    hlines(-1, solving_data['interval'][0], solving_data['interval'][-1], label='Tolerance = 10^-1')
    xlim(solving_data['interval'][0], solving_data['interval'][-1])
    xlabel('Number of points collected', fontsize=14)
    ylabel('log10', fontsize=14)

    index = 4-3
    plt.annotate("number of points = 4",
                xy=(solving_data['interval'][index]+0.01, logerr[index]+0.2), xycoords='data',
                xytext=(solving_data['interval'][index]+0.8, logerr[index]+4.5), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"),
                )
    grid()
    title('Calibration algorithm verification using simulated geometry')
    legend()
    show()
#----------------------------------------
if __name__ == '__main__':
    print "Sampling points..."    
    geometry_info = setup_geometry(plane, 300, num_points)

    print "Solving for dirx, diry, dirz and local_delta_vector..."
    result, s, cond_num, time_spent = find_solution(geometry_info)

    print
    print 'Time spent solving '+str(num_points)+' points: ' + str(time_spent) +' seconds.'

    print
    print 'Collecting solving information...'
    solving_data = perform_solution_run(geometry_info)

    print
    print 'Preparing plots...'
    make_plots(solving_data)
