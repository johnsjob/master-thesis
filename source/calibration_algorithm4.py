# -*- coding: cp1252 -*-
from __future__ import division
#----------------------------------------#
import numpy
import time
mat = numpy.array
from numpy.linalg import solve, det, inv, cond
#----------------------------------------
# custom imports
from helperfunctions_math import *
from helperfunctions_plot import *
import plane_relative as plane_tools
#----------------------------------------#
# num_points = 12 absolute minimum, actually 12+1
num_points = 120
#========================================#
# placing paper origin
o = mat([1000*rand(), 1000*rand(), 1000*rand()])

### defining paper orientation
r,t,s = -rand_range(-180,180), rand_range(-180,180), rand_range(-180, 180)

# define the paper-orientation in global (robot-system) directions
plane = plane_tools.define_plane_from_angles(o, r, t, s, 'local')
dirx = plane[:3,0]
diry = plane[:3,1]

# Orientation of the tool in tool (local) coordinate system 
local_tool_orientation = rotation_matrix_rot_tilt_skew(-10, 20, 30)
#################################################
# define delta vector which
# we want to find (in tool-space / local space)
L = 100
local_delta_vector = mat([1,2,3])
local_delta_vector = (local_delta_vector / norm(local_delta_vector))*L #length L
#----------------------------------------
def merge_dicts(*list_of_dicts):
    ret = {}
    keys = []
    for d in list_of_dicts:
        keys += d.keys()
    keys = set().union(keys)
    for k in keys:
        if not ret.has_key(k):
            ret[k] = []
        for d in list_of_dicts:
            if d.has_key(k):
                empty = False
                try:
                    empty = len(d[k]) == 0
                except:
                    pass
                if not empty:
                    ret[k].append( d[k] )
    for k in keys:
        if len(ret[k]) == 1:
            ret[k] = ret[k][0]
        elif len(ret[k]) == 0:
            del ret[k]
        else:
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
def solve_tool0_tip_alt(array_forward_kinematics_T44, array_anoto2D, array_lhs_sys_eq = None):
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
    result = mat([X,Y,Z,D]).T
    condition = cond(L)
    return result, condition
#----------------------------------------
def generate_Xflange_orientation(plane,rot, tilt, skew):
    """
        Generate Xtcp-orientation in world coordinates (Rrel)
        Planes are homoenous matrices, if we want the orientation
        we need the 0:3,0:3 submatrix.
    """
    return plane_tools.define_plane_relative_from_angles(plane, (0,0,0),
                                                         rot, tilt, skew,'global')[:3,:3]
#----------------------------------------
def setup_geometry(current_plane=None, point_spread=300, num_points=120):
    global local_delta_vector

    collected_info = {'plane':current_plane}
    collected_data = []
    #generating points and "forward-kinematics"
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
                                                generate_Xflange_orientation(collected_info['plane'],**info['angles'])
        #generate pen-tip position in Anoto2d in mm
        px,py = generate_random_Anoto_Point(point_spread)
        info['pentip_2d'] = [px,py]

        # generate global Xtcp position in mm
        info['Xtcp0'] = (plane_tools.get_plane_point(plane, px, py)[:3] - \
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
    collected_info['data'] = merge_dicts(*collected_data)
    return collected_info
#----------------------------------------
def find_solution(collected_info):
    start_time = time.clock()

    result, cond_num = solve_tool0_tip_alt(collected_info['data']['forward_kinematics'],
                                      collected_info['data']['pentip_2d'],
                                      problem_formulation)
    stop_time = time.clock()
    time_spent = stop_time - start_time


    # solve for orientation s which should be same as local_tool_orientation
    flange_orientation_reshaped = collected_info['data']['Xflange_orientation_relative_to_paper_plane'].reshape(360,3)

    lhs = matmul( flange_orientation_reshaped.T,
                  flange_orientation_reshaped)

    rhs = matmul( flange_orientation_reshaped.T,
                  collected_info['data']['global_tool_orientation'].reshape(360,3))

    solved_tool_orientation = linalg.solve(lhs, rhs)

    return result, solved_tool_orientation, cond_num, time_spent
#----------------------------------------
if __name__ == '__main__':
    print "Sampling points..."    
    collected_info = setup_geometry(plane)

    print "Solving for dirx, diry, local_delta_vector..."
    result, s, cond_num, time_spent = find_solution(collected_info)

    print
    print 'Time spent solving '+str(num_points)+' points: ' + str(time_spent) +' seconds.'

    print
    print 'Preparing plots...'
    comp = mat(list(plane[:3,:3].T.flatten()) +
               list(local_delta_vector)).reshape(4,3).T
    l_cond = []
    l_err = []
    for k in xrange(3,num_points):
        res, cond_num = solve_tool0_tip_alt(collected_info['data']['forward_kinematics'][0:k,:,:],
                                      collected_info['data']['pentip_2d'][0:k],
                                      problem_formulation)
        l_cond.append(cond_num)

        err = abs(comp-res)
        l_err.append(norm(err))    
    print 'solution error = ' + str(norm(err))

    t = range(3,num_points)
    logcond = numpy.log10(l_cond)
    logerr = numpy.log10(l_err)
    plot(t, logcond, label='Condition number',linewidth=2);
    plot(t, logerr, label='Error (frobenious norm)',linewidth=2);
    hlines(-1, t[0], t[-1], label='Tolerance = 10^-1')
    xlim(t[0], t[-1])
    xlabel('Number of points collected', fontsize=14)
    ylabel('log10', fontsize=14)

    index = 4-3
    plt.annotate("number of points = 4",
                xy=(t[index]+0.01, logerr[index]+0.2), xycoords='data',
                xytext=(t[index]+0.8, logerr[index]+4.5), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"),
                )
    grid()
    title('Calibration algorithm verification using simulated geometry')
    legend()
    show()
