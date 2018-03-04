import random

import os.path as path
from plotsettings import PlotSettings
from pylab import axhline, plot, show, \
                  axes, grid, xlabel, ylabel, title, xticks, yticks, savefig

from helperfunctions_math import mat, homogenous_matrix as hom, nzip,\
     rotation_matrix_skew_tilt_rot as ori

from denavit_hartenberg140 import forward_kinematics,\
     calc_valid_raw_invkin_irb140,\
     calc_valid_invkin_irb140,\
     calc_invkin_irb140,\
     DH_TABLE as dh_table

def plot_robot_geometry(robot_info, color='k'):
        global_robot_frames = mat(robot_info['robot_geometry_global'])
        plot(global_robot_frames[:,0,3],
                global_robot_frames[:,2,3], color ,linewidth=2)

        tool_pos = nzip(robot_info['tcp'][:3, 3],
                        robot_info['flange'][:3, 3]).T
        plot(tool_pos[:,0],
             tool_pos[:,2], color='g' ,linewidth=3)

def main():
        j1 =  0
        j2 =  0
        j3 =  0
        j4 =  0
        j5 =  0
        j6 =  0

        joint_values = [j1, j2, j3, j4, j5, j6]
        tool = hom(0,0,0,[0.1,0,0])

        dh_table['tool'] = tool
        info = forward_kinematics(*joint_values, **dh_table)
        # multiplication from left: rotation in base to base
        # multiplication from right: rotation in tool to base
        tcp = info['tcp'].dot(hom(ori(-10,30,40)))
        pose = hom(tcp[:3,:3],[0.6, 0, 0.3])

        s = calc_invkin_irb140(pose, raw_solutions=True)
        ik_up = forward_kinematics(*s[0], **dh_table)
        ik_down = forward_kinematics(*s[5], **dh_table)
        ik_up_back = forward_kinematics(*s[10], **dh_table)
        ik_down_back = forward_kinematics(*s[15], **dh_table)

        #plot_robot_geometry(info)
        plot_robot_geometry(ik_up,'b')
        plot_robot_geometry(ik_up_back,'b--')
        plot_robot_geometry(ik_down,'r')
        plot_robot_geometry(ik_down_back,'r--')
        xlabel('x [m]', fontsize=PlotSettings.label_size)
        ylabel('z [m]', fontsize=PlotSettings.label_size)
        xticks(fontsize=PlotSettings.tick_size)
        yticks(fontsize=PlotSettings.tick_size)
        grid()
        axes().set_aspect('equal', 'datalim')
        # save figure plot
        figpath = r"C:\Users\***REMOVED***\Dropbox\exjobb\rapport\images"
        savefig(path.join(
            figpath, "invkin-elbowupdown2.png"), bbox_inches='tight')

if __name__ == '__main__':
    main()
