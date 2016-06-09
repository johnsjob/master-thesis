# -*- coding: cp1252 -*-
from __future__ import division
#----------------------------------------#
import numpy
import numpy as n
import utils

from pylab import xlim, ylim
from numpy.linalg import solve, lstsq, det, inv, cond, svd
from numpy import array as mat, log10, diag
#----------------------------------------
# custom imports
from helperfunctions_math import *
from helperfunctions_plot import *
import plane_relative as plane_tools
#----------------------------------------#
##      --Refactoring information--
##        geometry_info = {'plane':current_plane}
##        geometry_info['local_tool_orientation']    : 3x3-dim numpy array
##        geometry_info['local_delta_vector']        : 3-dim numpy array
##        geometry_info['correct_solution_geometry'] : solution_tensor: 3x4-dim numpy array
##
##
##        geometry_info['angles'] = list_of:
##                                    {
##                                        rot_val : float
##                                        tilt_Val : float
##                                        skew_val : float
##                                    }
##        # Xtcp (flange) orientation in global space, generated relative to the paper plane
##        geometry_info['Xflange_orientation_relative_to_paper_plane']:
##                list_of: hom
##
##        #generate pen-tip position in Anoto2d in mm
##        geometry_info['pentip_2d']:
##                list_of: 2d tuple
##
##        # generate global Xtcp position in mm
##        geometry_info['Xtcp0']:
##            list_of: 3-dim numpy array
##
##        # generate relative-tool-orientation in world coordinates
##        geometry_info['global_tool_orientation']:
##            list_of: 3x3-dim of numpy array
##
##        geometry_info['forward_kinematics']:
##            list_of: hom
##            from: global_tool_orientation, Xtcp0

# num_points = 12 absolute minimum, actually 12+1
num_points = 500
#========================================#
# placing paper origin
o = mat([1000*rand(), 1000*rand(), 1000*rand()])

### defining paper orientation
r,t,s = -rand_range(-180,180), rand_range(-180,180), rand_range(-180, 180)

# define the paper-orientation in global (robot-system) directions
# however, the paper uses "euclidean" orientation ('local' = euclidean)
plane = plane_tools.define_plane_from_angles(o, r, t, s, 'local')

#################################################
# define the units so less mistakes are made
chosen_unit = 'mm'

unit_descriptor = {
    'm' : 1.0,
    'mm': 1000.0
    }

unit = unit_descriptor[chosen_unit]
metre = unit
millimetre = unit / 1000.0
#----------------------------------------
# define delta vector which we want to find (in tool-space / local space)
# in unit lengths
L = 100 * millimetre
local_delta_vector = mat([1,2,3])
local_delta_vector = (local_delta_vector / norm(local_delta_vector))*L #length L

# Orientation of the tool in tool (local) coordinate system 
local_tool_orientation = rotation_matrix_rot_tilt_skew(-10, 20, 30)
#----------------------------------------
# define the anoto point spread in unit lengths
#plane_point_spread = 47 * millimetre
plane_point_spread = 200 * millimetre

# define max-tilt of pen
pen_max_tilt = 40
#pen_max_tilt = 10

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
def vec_ang(v,w):
    res = matmul(v,w) / (norm(v) * norm(w))
    return numpy.arccos(res) * 180 / numpy.pi
#----------------------------------------
def vec_diff(v1, v2):
    err = norm(v1 - v2)
    norm_err = abs(norm(v1) - norm(v2))
    angle_err = rad_to_ang(acos( (v1/norm(v1)).dot((v2/norm(v2))) ))
    return err, norm_err, angle_err
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
def generate_random_Anoto_Point(L):
    px = L*rand()
    py = L*rand()
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
    #X, Y, D = lstsq(lhs,rhs)[0].reshape(3,3).T
    X = X / numpy.linalg.norm(X)
    Y = Y / numpy.linalg.norm(Y)
    Z = cross(X, Y)
#    Z = Z / numpy.linalg.norm(Z)
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
def setup_geometry(current_plane, point_spread, num_points, perturbations=None):
    global local_delta_vector, local_tool_orientation,\
           millimetre, pen_max_tilt

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
            'tilt': rand_range(0, pen_max_tilt),
            'skew': rand_range(-90, 180) 
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
        # perturbations
        if type(perturbations) in [list, tuple]:
            if 'tip' in perturbations:
                tilt = info['angles']['tilt']
                r = lambda: (0.6*rand()-0.3) * millimetre
                if tilt >= 0:
                    r2 = lambda: (0.02*rand()-0.01)*(1-abs(info['angles']['tilt'])/numpy.max(abs(info['angles']['tilt']))) * millimetre
                    info['pentip_2d'] = [px + r() + r2(), py + r() + r2()]
                else:
                    r2 = lambda: (0.02*rand()-0.01) * millimetre
                    info['pentip_2d'] = [px + r() + r2(), py + r() + r2()]
        # ^OK
        collected_data.append(info)
    geometry_info['data'] = merge_dicts(*collected_data)
    geometry_info['data']['angles'] = merge_dicts(*geometry_info['data']['angles'])

    return geometry_info
#----------------------------------------
def find_solution_pen_tip(geometry_info, included_solutions_from_start = -1):    
    result, cond_num = solve_tool0_tip_alt(geometry_info['data']['forward_kinematics'][:included_solutions_from_start],
                                      geometry_info['data']['pentip_2d'][:included_solutions_from_start],
                                      problem_formulation)
    return result, cond_num

##def find_solution_pen_ori(geometry_info, included_solutions_from_start = -1):
##    # solve for orientation s which should be same as local_tool_orientation
##    l,m,n = geometry_info['data']['Xflange_orientation_relative_to_paper_plane'].shape
##    flange_orientation_reshaped = geometry_info['data']['Xflange_orientation_relative_to_paper_plane'].reshape(l*m,n)
##
##    lhs = matmul( flange_orientation_reshaped.T,
##                  flange_orientation_reshaped)
##
##    l,m,n = geometry_info['data']['global_tool_orientation'].shape
##    rhs = matmul( flange_orientation_reshaped.T,
##                  geometry_info['data']['global_tool_orientation'].reshape(l*m,n))
##
##    solved_tool_orientation = linalg.solve(lhs, rhs)
##
##    #normalize result
##    solved_tool_orientation[:,0] = solved_tool_orientation[:,0] / norm(solved_tool_orientation[:,0])
##    solved_tool_orientation[:,1] = solved_tool_orientation[:,1] / norm(solved_tool_orientation[:,1])
##    solved_tool_orientation[:,2] = solved_tool_orientation[:,2] / norm(solved_tool_orientation[:,2])
##    return solved_tool_orientation, cond(lhs)

def _solve_orientation(As, Bs):
    Ac    = As - n.mean(As, axis=0)
    Bec    = Bs - n.mean(Bs, axis=0)
    h     = lop(n.outer, Ac, Bec)
    H     = n.sum(h, axis=0)
    U,S,V = svd(H)
    # solve for solution tensor of order 9x9
    # in this tensor the best solution resides in index 0,4,8 in the tensor form
    # of 9x3x3 which corresponds to diagonal elements of the 3x3x3x3 solution tensor
    solution_tensor = U.dot(V)
    
    # best column-solutions resides in the diagonal of the 3x3x3x3 solution tensor
    C1, C2, C3 = solution_tensor[0::3, 0::3],\
                 solution_tensor[1::3, 1::3],\
                 solution_tensor[2::3, 2::3]

    D1, D2, D3 = solution_tensor[0:3,0:3],\
                 solution_tensor[3:6,3:6],\
                 solution_tensor[6:9,6:9]
    solution1 = mat([C1.T, C2.T, C3.T])
    solution2 = mat([D1, D2, D3])
    return solution1, solution2, solution_tensor, [U,S,V]

def find_solution_pen_ori(geometry_info, included_solutions_from_start = -1):
    # solve for orientation s which should be same as local_tool_orientation
    #l,m,n = geometry_info['data']['Xflange_orientation_relative_to_paper_plane'].shape

    flange_orientation = geometry_info['data']['forward_kinematics'][:,:3,:3]
    pen_orientation = geometry_info['data']['global_tool_orientation'][:,:3,:3]

    solved_tool_orientation = _solve_orientation(flange_orientation, pen_orientation)

##    #normalize result
##    solved_tool_orientation[:,0] = solved_tool_orientation[:,0] / norm(solved_tool_orientation[:,0])
##    solved_tool_orientation[:,1] = solved_tool_orientation[:,1] / norm(solved_tool_orientation[:,1])
##    solved_tool_orientation[:,2] = solved_tool_orientation[:,2] / norm(solved_tool_orientation[:,2])
    return solved_tool_orientation, 0.0
#----------------------------------------
def perform_solution_run(geometry_info):
    interval = range(3,num_points)
    list_of_solving = []
    for k in interval:
        solve_info = {}
        tip_wobj_res = find_solution_pen_tip(geometry_info, k)

##        solve_info['point_spread_x'] = numpy.std(geometry_info['data']['pentip_2d'][:k], axis=0)[0]
##        solve_info['point_spread_y'] = numpy.std(geometry_info['data']['pentip_2d'][:k], axis=0)[1]
        
        solve_info['tipwobj-result']   = tip_wobj_res[0]
        solve_info['tip-result']       = tip_wobj_res[0][:,3]
        solve_info['wobj-result']      = tip_wobj_res[0][:,:3]
        solve_info['tip-cond_num']     = tip_wobj_res[1]
        
        solve_info['orientation-result'], solve_info['orientation-cond_num'] = find_solution_pen_ori(geometry_info, k)
        sol1, sol2, tens, (_u,_s,_v) = solve_info['orientation-result']
        solve_info['orientation-result'] = sol2

        solve_info['err-tipwobj'] = abs(geometry_info['correct_solution_geometry'] - solve_info['tipwobj-result'])
        solve_info['err-tip']     = numpy.linalg.norm(solve_info['err-tipwobj'][:,3])
        solve_info['err-wobj']    = numpy.linalg.norm(solve_info['err-tipwobj'][:,:3])
        
        solve_info['err-ori']     = numpy.linalg.norm(geometry_info['local_tool_orientation'] - sol2)
        list_of_solving.append(solve_info)

    solving_data = merge_dicts(*list_of_solving)
    solving_data['interval'] = interval
    print 'solution max  error tip(1) = {}\n'.format( numpy.max( abs( solving_data['err-tip'][1:]) ))

    print 'solution max  error tip(20) = {}'.format( numpy.max( abs( solving_data['err-tip'][20:]) ))
    print 'solution mean error tip(20) = {}\n'.format( numpy.mean( abs( solving_data['err-tip'][1:21]) ))

    print 'solution max  error tip(40) = {}'.format( numpy.max( abs( solving_data['err-tip'][40:]) ))
    print 'solution mean error tip(40) = {}\n'.format( numpy.mean( abs( solving_data['err-tip'][1:41]) ))
    print 'solution error ori = {}'.format( numpy.max( abs( solving_data['err-ori'][1:]) ))
    return solving_data
#----------------------------------------
def make_plots(solving_data):
    global chosen_unit
    
    logcond = log10( solving_data['tip-cond_num'] )
    plot(solving_data['interval'], logcond,
                       'b--',label='Condition number tip/wobj',
                       linewidth=2)
    
    logerr  = log10( solving_data['err-tip'] )
    plot(solving_data['interval'], logerr,
                      'b',label='Error tip (frobenious norm)',
                      linewidth=2);

    logerr  = log10( solving_data['err-wobj'] )
    plot(solving_data['interval'], logerr,
                      'g',label='Error wobj (frobenious norm)',
                      linewidth=2);

##    logerr  = log10( solving_data['err-ori'] )
##    plot(solving_data['interval'], logerr,
##                      'r',label='Error ori (frobenious norm)',
##                      linewidth=2);
    
    if chosen_unit == 'mm':
        tol = -1
        hlines(tol, solving_data['interval'][0],
                   solving_data['interval'][-1],
                   label='Tolerance = 10^{}'.format(tol))
    else:
        tol = -4
    hlines(-4, solving_data['interval'][0],
               solving_data['interval'][-1])
    
    xlim(solving_data['interval'][0], solving_data['interval'][-1])
    xlabel('Number of points collected', fontsize=14)
    ylabel('log10'.format(chosen_unit), fontsize=14)

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
    res = []
    ks = range(1)
    for k in ks:
        print 'Run {} % complete!'.format(100* k / len(ks))
        with utils.timing.Timer() as timer:
            try:
                print "Sampling points..."    
                geometry_info = setup_geometry(plane, plane_point_spread,
                                               num_points, perturbations=['tip'])
                print 'Collecting solving information...'
                solving_data = perform_solution_run(geometry_info)
                res.append(solving_data)
                print
                print 'Preparing plots...'
            except Exception as e:
                print str(e)
                raise
  #              continue

    make_plots(solving_data)
    if chosen_unit == 'mm':
        length_tol = -1 # 1/10th millimetre
    else:
        length_tol = -4 # 1/10th millimetre (in meters)

    for key,tol, un in zip(['err-tip', 'err-wobj'],
                           [length_tol, -4],
                           [chosen_unit, '']):
        val = mat( [x[key] for x in res] )
        maxval = numpy.max(val, axis=0)
        minval = numpy.min(val, axis=0)
        meanval = numpy.mean(val, axis=0)
        if un:
            unit_str = '[{}]'.format(un)
        else:
            unit_str = ''
        plot(log10(maxval),'b', label = 'max {} {}'.format(key, unit_str))
        plot(log10(meanval),'g', label = 'mean {} {}'.format(key, unit_str))
        plot(log10(minval),'r', label = 'min {} {}'.format(key, unit_str))
##        sx = numpy.mean(mat( [x['point_spread_x'] for x in res] ), axis=0)
##        sy = numpy.mean(mat( [x['point_spread_y'] for x in res] ), axis=0)
##        plot(log10(sx),'k', label = 'sx')
##        plot(log10(sy),'k', label = 'sy')
        hlines(tol, solving_data['interval'][0],
                    solving_data['interval'][-1],
                    label='Tolerance = 10^{} {}'.format(tol, unit_str))
        legend()
        xlim(solving_data['interval'][0], solving_data['interval'][-1])
        xlabel('Number of measured points', fontsize=14)
        ylabel('log10', fontsize=14)
        if 'tip' in key:
            ylim(-4, 4)
        elif 'wobj' in key:
            ylim(-6, 1)
        title('Calibration algorithm verification with repetition using simulated geometry')
        grid()
        show()
