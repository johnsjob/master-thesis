# -*- coding: utf-8 -*-
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import axes

from helperfunctions_plot import plot_equal_perspective

import numpy as np
from numpy import array as mat
from numpy import linalg as li
from numpy.linalg import det

import sys


#debug

class StPlot:
    
    def __init__(self, length_unit='m'):
        self.length_factors = {
            'm': 1.0,
            'mm': 1000.0
            }
        self.length_unit = length_unit
        self.length_factor = self.length_factors[self.length_unit]

        self.app = None
        self.plots = []        
        self.fig = self.__initWindow()
        self.__setupEnvironment(self.fig)
        return

    def __initWindow(self):
        # Twice as wide as it is tall.
        fig = plt.figure(1)#,figsize=plt.figaspect(0.5))
        return fig

    def __setupEnvironment(self, fig):
        #---- First subplot
        ax_3d = fig.add_subplot(1, 1, 1, projection='3d')

        cube_limits = mat([[0, 0.55],
                          [-0.25, 0.25],
                          [0.0, 0.55]])*self.length_factor

        plot_equal_perspective(ax_3d, cube_limits[0],
                                      cube_limits[1],
                                      cube_limits[2])
        self.plots.append(ax_3d)
        return

    def draw_joint_paths(self, J):
        fig = plt.figure(2)
        for i in range(6):
            ax = fig.add_subplot(6, 1, i+1)
            ax.plot(J[:,i])
            plt.legend(['j{}'.format(i+1)])

    def draw_joint_velocities(self, Jvel):
        fig = plt.figure(3)
        for i in range(6):
            ax = fig.add_subplot(6, 1, i+1)
            ax.plot(Jvel[:,i])
            plt.legend(['w{}'.format(i+1)])

    def draw_jacobian_determinants(self, jacobs):
        fig = plt.figure(4)
        ax = fig.add_subplot(1, 1, 1)
        D = map(det, jacobs)
        ax.plot(D)
        plt.legend(['det(J)'])
        
        
    def draw_frames(self, frames, size=1.0, **kwargs):
        for frame in frames:
            self.draw_frame(frame, size, **kwargs)
        return

    def draw_trajectory(self, poses, **kwargs):
        self.draw_frames(poses, size = 0.005, **kwargs)
        self.draw_curve(poses[:,:3,3], linewidth = 2, **kwargs)

    def draw_robot(self, robot_geometry):
        for i,robot_frame in enumerate(robot_geometry):
            self.draw_frame(robot_frame, size=0.02, linewidth=2)
            if i>0:
                a = robot_geometry[i-1][:3,3]
                b = robot_geometry[i][:3,3]
                k = 1.0/7*(i+1)
                c = (k, 0, k)
                self.draw_line2(a, b, color = c, linewidth = 2)
        return

    def draw_solution(self, points):
        ax = self.plots[1]
        ax.plot(*points.T)

    def draw_curve(self, points, **kwargs):
        ax = self.plots[0]
        ax.plot(points[:,0],
                points[:,1],
                points[:,2], color='0.75', **kwargs)
        return

    def draw_frame(self, frame_matrix, size=1.0, **kwargs):
        """
        Renders 4x4 homogenous matrices of the form [[R, t]
                                                     [0, 1]]
        """
        size = size * self.length_factor

        dirx = frame_matrix[:3,0] * size
        diry = frame_matrix[:3,1] * size
        dirz = frame_matrix[:3,2] * size
        px, py, pz = frame_matrix[:3,3]

        self.draw_direction(px, py, pz, *dirx, color='r', **kwargs)
        self.draw_direction(px, py, pz, *diry, color='g', **kwargs)
        self.draw_direction(px, py, pz, *dirz, color='b', **kwargs)
        return

        
    def draw_tool(self, robot_flange, tool):
        ax = self.plots[0]
        tcp = robot_flange.dot(tool)

        self.draw_line2(robot_flange[:3,3],
                        tcp[:3,3], linewidth=2)
#        ax.plot([x],[y],[z],'bo')
#        ax.plot([x2],[y2],[z2],'ro')
        self.draw_frame(tcp, size=0.01, linewidth=2)
        return
        
    
    def draw_line(self, x0, y0, z0,
                  x01, y01, z01, **kwargs):

        line = np.array(zip([x0,x01],
                            [y0,y01],
                            [z0,z01]))
        ax = self.plots[0]
        ax.plot(*line.T, **kwargs)
        return

    def draw_line2(self, v0, v1,
                  color='b', **kwargs):
        x0, y0, z0 = v0
        x01, y01, z01 = v1
        line = np.array(zip([x0,x01],
                            [y0,y01],
                            [z0,z01]))
        ax = self.plots[0]
        ax.plot(*line.T, **kwargs)
        return

    def draw_direction(self, x,y,z, dx,dy,dz, **kwargs):
        self.draw_line(x, y, z,
                       x + dx,
                       y + dy,
                       z + dz, **kwargs)
        return
        

    def plot(self, y_values):
##            def draw_curve(self, points, **kwargs):
        ax = self.plots[-1]
        ax.plot(*points.T, color='0.75', **kwargs)


    def show(self):
        plt.show()
        return



if __name__ == '__main__':
    plot = StPlot()
    plot.show()
