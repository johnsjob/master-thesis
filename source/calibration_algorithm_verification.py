import random

import numpy, numpy as n
from numpy import pi, linspace, meshgrid as mesh, zeros,\
     arccos as acos, log10
from numpy.linalg import norm, det, inv
from helperfunctions_plot import *
import calibration_algorithm as cal
from helperfunctions_math import rand_range,\
                                 rotation_matrix_rot_tilt_skew as ori,\
                                 rotation_matrix_x_y_z as xyz,\
                                 homogenous_matrices, nzip, nmap,\
                                 quat_to_rot, rot_to_quat,\
                                 rot_tilt_skew as rts
from pylab import axhline
from plane_relative import generate_symmetric_curve,\
                           get_transformed_points, attach_to_base_frame,\
                           create_circle_curve, place_curve, attach_frames,\
                           orientation_frames

from denavit_hartenberg140 import forward_kinematics, DH_TABLE as dh_table,\
     calc_valid_invkin_irb140 as invkin

from denavit_hartenberg import homogenous_matrix as hom
from jacobian import jacobian_from_joints

from standardplot import StPlot
import pylab as plt

import itertools as it

import sys
import time
import json

import utils

from collections import OrderedDict


numpy.set_printoptions(precision=4)
numpy.set_printoptions(suppress=True)

# generate angles
num_points = 50
rot  = numpy.linspace(0,180,num_points)
tilt = numpy.linspace(-40,40,num_points)
skew = numpy.linspace(0,0,num_points)
angles = nzip(rot, tilt, skew)

normalize = lambda x: x / norm(x)

def print_dict(d):
    print json.dumps(d, indent=4)
    
def targets(xmin,xmax,ymin,ymax, num):
    '''
    Generates a ('num'X'num') rectangle grid of targets in the range of
    [xmin, xmax], [ymin, ymax]. Depending on application these can either be
    relative to a work boject or absolute.

    return:
        ndarray : (num**2, 4, 4)
    '''
    x = linspace(xmin,xmax,num)
    y = linspace(ymin,ymax,num)
    xy = mat(mesh(x,y)).T
    res = nmap(lambda row: nmap(lambda pair: hom(0,0,0,[pair[0],pair[1],0]) ,row), xy)
    a,b,c,d = res.shape
    res = res.reshape((a*b,c,d))
    return res

def ik(pose):
    return invkin(pose)[0]

def flange(pose):
    '''
    Given a pose, calulate the forward-kinematics frame obtained
    from reaching the point by inverse-kinematics.
    '''
    return forward_kinematics(*ik(pose), **dh_table)['flange']

def pen_pos(pose):
    return pose[:2,3]

def e_hom_flange(mag):
    '''
        noise function used to generate a homogenous matrix with(-mag, mag)
        noise on the translation part.
    '''
    res = zeros((4,4))
    mag = mag*1e-3
    pos = mat([rand_range(-mag, mag) for x in range(3)])
    res[:3,3] = pos
    return res

def solution_run(geometry_info):
    num_meas = len(geometry_info['data']['forward_kinematics'])
    res = [ cal.find_solution_pen_tip(geometry_info, k)
            for k in range(3, num_meas) ]
    return mat(res)

def solution_errors(solution, comparison):
    c = solution[1]
    solution = solution[0]
    wobj = comparison[:,:3]
    tool = comparison[:,3]

    tool_err = abs(tool - solution[:,3]) * 1000
    wobj_err = abs(wobj - solution[:,:3])
    errors = {
        'wobj': {
            'norm': norm(wobj_err),
            'x_ang': acos(wobj[:,0].dot(normalize(solution[:,0])))*180.0/pi,
            'y_ang': acos(wobj[:,1].dot(normalize(solution[:,1])))*180.0/pi,
            'n_ang': acos(wobj[:,2].dot(normalize(solution[:,2])))*180.0/pi,
            'unit': 'deg'
            },
        'tool': {
            'norm': norm(tool_err),
            'x': tool_err[0],
            'y': tool_err[1],
            'z': tool_err[2],
            'unit': 'mm'
            },
        'cond': c
        }
    return errors

if __name__ == '__main__':
    # init - tool
    joints = (0,0,0,0,0,0)
    dh_table['tool'] = hom(0,0,0,
                           [0.05,0.05,0.1])
    dh_table['tool'][:3,:3] = ori(50,-45,0).dot(dh_table['tool'][:3,:3])
    dh_table['tool'][:3,:3] = dh_table['tool'][:3,:3].dot(ori(0,180,90)).dot(xyz(5,-10,0))
    tool = dh_table['tool']

    # init - wobj
    wobj = hom(180-45,180,0,[0.2,-0.2,-0.2])
    wobj[:3,:3] = wobj[:3,:3].dot(ori(0,180,90))

    # init reference
    ref = zeros((3,4))
    ref[:3,:3] = wobj[:3,:3]
    ref[:,3] = tool[:3,3]

    # init - relative targets
    tars = targets(0, 10e-3,
                   0, 10e-3, num=10)
    tars = nmap(lambda x: x.dot(hom(ori(0,rand_range(1,20),rand_range(-10,10)))), tars)

    # init - global targets
    res = nmap(wobj.dot, tars)

    # calculate flanges from global targets
    flanges = nmap(flange, res)
    flanges = nmap(lambda x: x + e_hom_flange(0.03), flanges)

    # obtain relative pen positions from local targets
    penxy = nmap(pen_pos, tars)

    # convert to 'anoto coords'
    penxy = nmap(lambda x: x + rand_range(0,1e-4)-0.5e-4, penxy)

    #penxy = nmap(lambda x: (rand_range(-100e-3,100e-3),rand_range(-100e-3,100e-3)), penxy)
    penxy = penxy*1000/0.3 + 1e9
    
    # supply information to calibration algorithm
    geometry_info = {
            'data':
                {
                    'forward_kinematics': flanges,
                    'pentip_2d': penxy
                }
        }
    # perform calibration, present results
    solution = cal.find_solution_pen_tip(geometry_info)
    d = solution_errors(solution, ref)
    print_dict(d)
    run = solution_run(geometry_info)
    errors = map(lambda x: solution_errors(x, ref),run)
    plot(log10([e['cond'] for e in errors]))
    plot(log10([e['tool']['norm'] for e in errors]))
    axhline(-1,color='k')
    grid()
    show()
    

##    # plot setup
##    pl = StPlot()
##    robot_info = forward_kinematics(*joints, **dh_table)
##    pl.draw_robot(robot_info['robot_geometry_global'])
##    pl.draw_tool(robot_info['flange'],
##                 dh_table['tool'])
##    pl.draw_frame(wobj, size=0.1)
##    map(lambda x: pl.draw_frame(x, size=0.01), res)
##    pl.show()

