# -*- coding: utf-8 -*-
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from helperfunctions_plot import plot_equal_perspective


import numpy as np
from numpy import array as mat
from numpy import linalg as li

import sys


class StPlot:

    def __init__(self):
        self.app = None

        self.plots = []        
        self.fig = self.__initWindow()
        self.__setupEnvironment(self.fig)

        return

    def draw_frames(self, frames, size=1.0, **kwargs):
        for frame in frames:
            self.draw_frame(frame, size, **kwargs)
        return

    def draw_trajectory(self, poses, **kwargs):
        self.draw_frames(poses, size = 0.01, **kwargs)
        self.draw_curve(poses[:,:3,3], linewidth = 2, **kwargs)

    def draw_robot(self, robot_geometry):
        for i,robot_frame in enumerate(robot_geometry):
            self.draw_frame(robot_frame, size=0.1)
            if i>0:
                a = robot_geometry[i-1][:3,3]
                b = robot_geometry[i][:3,3]
                k = 1.0/7*(i+1)
                c = (k, 0, k)
                self.draw_line2(a, b, color = c, linewidth = 4)
        return

    def draw_curve(self, points, **kwargs):
        ax = self.plots[0]
        ax.plot(*points.T, color='0.75', **kwargs)
        return

    def draw_frame(self, frame_matrix, size=1.0, **kwargs):
        """
        Renders 4x4 homogenous matrices of the form [[R, t]
                                                     [0, 1]]
        """
        dirx = frame_matrix[:3,0] * size
        diry = frame_matrix[:3,1] * size
        dirz = frame_matrix[:3,2] * size
        px, py, pz = frame_matrix[:3,3]

        self.draw_direction(px, py, pz, *dirx, color='b', **kwargs)
        self.draw_direction(px, py, pz, *diry, color='g', **kwargs)
        self.draw_direction(px, py, pz, *dirz, color='r', **kwargs)
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
        
    def __initWindow(self):
        # Twice as wide as it is tall.
        fig = plt.figure(1)#,figsize=plt.figaspect(0.5))
        return fig

    def __setupEnvironment(self, fig):
        #---- First subplot
        ax_3d = fig.add_subplot(1, 2, 1, projection='3d')
        plot_equal_perspective(ax_3d,[-0.7, 0.7],
                                     [-0.7, 0.7],
                                     [-0.5, 1])
        self.plots.append(ax_3d)

        #---- Second subplot
        ax_2d = fig.add_subplot(1, 2, 2)
        self.plots.append(ax_2d)
        return


    def show(self):
        plt.figure(2,figsize=plt.figaspect(0.5))
        plt.show()
        return



if __name__ == '__main__':
    plot = StPlot()
    plot.show()
