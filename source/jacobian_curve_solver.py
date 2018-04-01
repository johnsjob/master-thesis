import random

from plotsettings import PlotSettings

import numpy, numpy as n
from numpy import pi, linspace
from numpy.linalg import norm, det, inv
from helperfunctions_plot import *
from helperfunctions_math import rand_range,\
                                 rotation_matrix_rot_tilt_skew as ori,\
                                 homogenous_matrices, nzip, nmap,\
                                 quat_slerp

from pylab import axhline, xticks, yticks, savefig, clf
from pylab import minorticks_on, subplot, subplots_adjust

from plane_relative import generate_symmetric_curve,\
                           generate_curve,\
                           get_transformed_points, attach_to_base_frame,\
                           create_circle_curve, place_curve, attach_frames,\
                           orientation_frames

from denavit_hartenberg140 import forward_kinematics, DH_TABLE as dh_table,\
     inverse_kinematics_curve, find_single_path, \
     calc_valid_invkin_irb140 as invkin

from denavit_hartenberg import homogenous_matrix as hom
from jacobian import jacobian_from_joints

#from pyqtplot import QtPlot
from standardplot import StPlot
import pylab as plt

import itertools as it

import sys
import time

#test

import utils

##sys.path.append('../int/djikstra/')
##from graph import shortestPath as shortest_path

#===============================================================================

numpy.set_printoptions(precision=4)
numpy.set_printoptions(suppress=True)

# generate angles
num_points = 50
rot  = numpy.linspace(0,0,num_points)
tilt = numpy.linspace(0,0,num_points)
skew = numpy.linspace(0,0,num_points)
angles = nzip(rot, tilt, skew)

normalize = lambda x: x / norm(x)

#===============================================================================

def calc_robot_tcp(*joints):
    return forward_kinematics(*joints, **dh_table)['tcp']


def custom_plot(*args, **kwargs):
    plot(*args, **kwargs)
    grid(b=True, which='major', color='k', linestyle='-', linewidth=0.3)
    grid(b=True, which='minor', color='k', linestyle='-', linewidth=0.1)
    minorticks_on()

def calc_robot_curve_j1():
    N = 100.0
    T = 10.0
    L = pi/2.0
    djs = linspace(0,L*180/pi,N)
    dts = linspace(0,T,N)

    djdts = djs/dts
    print djdts
    print djdts*dts

    dt = T / N # seconds / dx
    rads_per_second = L / T

    tcps = mat([calc_robot_tcp(v,10,20,0,90,0) for v in linspace(0, L*180/pi, N)])
    #w = -tcps[:,:3,2]*rads_per_second
    w = mat([[0,0,1]*int(N)]).reshape(N,3)*rads_per_second
    print w.shape
    r = tcps[:,:3,3]
    v = nmap(lambda x: reduce(n.cross, x), zip(w, r))
    return tcps, nzip(v,w).reshape(N,6),dts

#===============================================================================

XLABEL = 'Time [s]'
YLABEL_W = 'Joint angular velocity [deg / s]'
YLABEL_Q = 'Joint angle [deg]'

#===============================================================================

def main():
    dh_table['tool'] = hom(0,0,0,[0.0,0.0,0.1])
    wobj = hom(-90,180,0,[0.6,0,0.5])
    # robot movement
    ##robot_movement_j1, vw, dts = calc_robot_curve_j1()
    ##trajectory = robot_movement_j1

    curve = generate_symmetric_curve()
    N = len(curve)
    T = 12.56
    L = 2*pi
    dt = T / N # seconds / dx
    rads_per_second = L / T

    frame = calc_robot_tcp(0,0,0,0,90,0)
    center = frame[:3,3]
    curve = nmap(frame.dot, curve)
    r = curve[:,:3] - center
    frames = n.tile(frame,(50,1)).reshape(50,4,4)
    frames[:,:,3] = curve
    trajectory = frames
    t = linspace(0, T, N)
    w = mat([[0,0,-1]]*len(trajectory))*rads_per_second
    v = nmap(lambda x: reduce(n.cross, x), zip(w,r))
    w = mat([[0,0,0]]*len(trajectory))
    vw = n.hstack((v,w))

    #inverse kinematics over curve
    result = inverse_kinematics_curve(trajectory)
    path = find_single_path(result)
    J = nmap(lambda x: jacobian_from_joints(*x), path)
    Jinv = nmap(inv, J)
    joint_angular_vel = nmap(lambda x: reduce(n.dot, x), zip(Jinv, vw))*180/pi

    print 'Path found: \n {}'.format(path)
    print 'Joint angular velocities: \n {}'.format(joint_angular_vel)

    index=0
    subplots_adjust(hspace=0.55)
    subplot(211)
    custom_plot(n.round(path[:,index],decimals=3), 'b', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_Q, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$q_'+str(index+1)+'$'], fontsize=PlotSettings.legend_size*1.25)
    subplot(212)
    custom_plot(n.round(joint_angular_vel[:,index],decimals=3), 'r', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_W, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$\dot{q}_'+str(index+1)+'$'], loc='lower right', fontsize=PlotSettings.legend_size*1.25)
    maximize_plot()
    savefig('C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\inverse_kinematics_over_curve\\q1.png',
            bbox_inches='tight', pad_inches=0)
    clf()

    index=1
    subplots_adjust(hspace=0.55)
    subplot(211)
    custom_plot(n.round(path[:,index],decimals=3), 'b', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_Q, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$q_'+str(index+1)+'$'], fontsize=PlotSettings.legend_size*1.25)
    subplot(212)
    custom_plot(n.round(joint_angular_vel[:,index],decimals=3), 'r', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_W, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$\dot{q}_'+str(index+1)+'$'], loc='upper right', fontsize=PlotSettings.legend_size*1.25)
    maximize_plot()
    savefig('C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\inverse_kinematics_over_curve\\q2.png',
            bbox_inches='tight', pad_inches=0)
    clf()

    index=2
    subplots_adjust(hspace=0.55)
    subplot(211)
    custom_plot(n.round(path[:,index],decimals=3), 'b', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_Q, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$q_'+str(index+1)+'$'], loc='lower right', fontsize=PlotSettings.legend_size*1.25)
    subplot(212)
    custom_plot(n.round(joint_angular_vel[:,index],decimals=3), 'r', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_W, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$\dot{q}_'+str(index+1)+'$'], loc='lower right', fontsize=PlotSettings.legend_size*1.25)
    maximize_plot()
    savefig('C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\inverse_kinematics_over_curve\\q3.png',
            bbox_inches='tight', pad_inches=0)
    clf()

    index=3
    subplots_adjust(hspace=0.55)
    subplot(211)
    custom_plot(n.round(path[:,index],decimals=3), 'b', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_Q, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$q_'+str(index+1)+'$'], fontsize=PlotSettings.legend_size*1.25)
    subplot(212)
    custom_plot(n.round(joint_angular_vel[:,index],decimals=3), 'r', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_W, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$\dot{q}_'+str(index+1)+'$'], loc='upper right', fontsize=PlotSettings.legend_size*1.25)
    maximize_plot()
    savefig('C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\inverse_kinematics_over_curve\\q4.png',
            bbox_inches='tight', pad_inches=0)
    clf()

    index=4
    subplots_adjust(hspace=0.55)
    subplot(211)
    custom_plot(n.round(path[:,index],decimals=3), 'b', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_Q, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$q_'+str(index+1)+'$'], fontsize=PlotSettings.legend_size*1.25)
    subplot(212)
    custom_plot(n.round(joint_angular_vel[:,index],decimals=3), 'r', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_W, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$\dot{q}_'+str(index+1)+'$'], loc='lower right', fontsize=PlotSettings.legend_size*1.25)
    maximize_plot()
    savefig('C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\inverse_kinematics_over_curve\\q5.png',
            bbox_inches='tight', pad_inches=0)
    clf()

    index=5
    subplots_adjust(hspace=0.55)
    subplot(211)
    custom_plot(n.round(path[:,index],decimals=3), 'b', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_Q, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$q_'+str(index+1)+'$'], fontsize=PlotSettings.legend_size*1.25)
    subplot(212)
    custom_plot(n.round(joint_angular_vel[:,index],decimals=3), 'r', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel(YLABEL_W, fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    legend(['$\dot{q}_'+str(index+1)+'$'], loc='lower right', fontsize=PlotSettings.legend_size*1.25)
    maximize_plot()
    savefig('C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\inverse_kinematics_over_curve\\q6.png',
            bbox_inches='tight', pad_inches=0)
    clf()

    custom_plot(nmap(det, J), 'b', linewidth=1.5)
    xlabel(XLABEL, fontsize=PlotSettings.label_size)
    ylabel('Jacobian determinant value', fontsize=PlotSettings.label_size)
    xticks(range(0,50,10)+[49],linspace(0,T,6), fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    xlim((0,49))
    maximize_plot()
    savefig('C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\inverse_kinematics_over_curve\\J.png',
            bbox_inches='tight', pad_inches=0)
    clf()

    pl = StPlot()
    joints = path[13]
    robot_info = forward_kinematics(*joints, **dh_table)
    pl.draw_robot(robot_info['robot_geometry_global'])
    pl.draw_trajectory(trajectory)
    pl.draw_tool(robot_info['flange'],
                 dh_table['tool'])
#    pl.show()

if __name__ == '__main__':
    main()
