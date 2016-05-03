import random

from pylab import axhline, plot, show, axes
from helperfunctions_math import mat
from denavit_hartenberg140 import forward_kinematics,\
     calc_valid_raw_invkin_irb140,\
     calc_valid_invkin_irb140,\
     calc_invkin_irb140,\
     DH_TABLE as dh_table

def plot_robot_geometry(robot_info, color='k'):
        global_robot_frames = mat(robot_info['robot_geometry_global'])
        plot(global_robot_frames[:,0,3],
                global_robot_frames[:,2,3], color ,linewidth=2)

if __name__ == '__main__':

        j1 =  0
        j2 =  10
        j3 =  10
        j4 =  0
        j5 =  10
        j6 =  0

        joint_values = [j1, j2, j3, j4, j5, j6]

        info = forward_kinematics(*joint_values, **dh_table)
        pose = info['T44']

        s = calc_invkin_irb140(pose, raw_solutions=True)
        ik_up = forward_kinematics(*s[0], **dh_table)
        ik_down = forward_kinematics(*s[5], **dh_table)
        ik_up_back = forward_kinematics(*s[10], **dh_table)
        ik_down_back = forward_kinematics(*s[15], **dh_table)
        
        plot_robot_geometry(info)
        plot_robot_geometry(ik_up,'b')
        plot_robot_geometry(ik_up_back,'b--')
        plot_robot_geometry(ik_down,'r')
        plot_robot_geometry(ik_down_back,'r--')
        axes().set_aspect('equal', 'datalim')
        show()
