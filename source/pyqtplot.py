# -*- coding: utf-8 -*-
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import pyqtgraph as pg

import numpy as np
from numpy import array as mat
from numpy import linalg as li

import sys


class QtPlot:

    def __init__(self):
        self._plot_color = {
                        'b': pg.glColor(0,0,255),
                        'g': pg.glColor(0,255,0),
                        'r': pg.glColor(255,0,0),
                        'y': pg.glColor(255,255,0),
                        'gr': pg.glColor(255*3/4,255*3/4,255*3/4)
                      } 

        self.app = QtGui.QApplication([])
        self.window = self.__initWindow()
        self.__setupEnvironment(self.window)

        win = pg.GraphicsWindow(title="Basic plotting examples")
        win.setWindowTitle('pyqtgraph example: Plotting')
        self.win = win
        # Enable antialiasing for prettier plots
        pg.setConfigOptions(antialias=True)

        self.plots = []

##        self.add_subplot(2,2,1)
##        self.plot([4,5,6])
##        self.add_subplot(2,2,4)
##        self.plot([6,5,4])
        #self.plot(2,2,0)
        return

    def draw_frames(self, frames, size=1.0, **kwargs):
        for frame in frames:
            self.draw_frame(frame, size, **kwargs)
        return

    def draw_trajectory(self, poses, **kwargs):
        self.draw_frames(poses, size = 0.01, **kwargs)
        self.draw_curve(poses[:,:3,3], width = 2, **kwargs)

    def draw_robot(self, robot_geometry):
        for i,robot_frame in enumerate(robot_geometry):
            self.draw_frame(robot_frame, size=0.1)
            if i>0:
                a = robot_geometry[i-1][:3,3]
                b = robot_geometry[i][:3,3]
                k = 255/7*(i+1)
                c = (k, 0, k)
                self.draw_line2(a, b, col = c, width=6)
        return

    def draw_curve(self, points, **kwargs):
        if not kwargs.has_key('col'):
            kwargs['col'] = 'gr'

        col = kwargs['col']
        kwargs.pop('col')

        plt = gl.GLLinePlotItem(pos=points, color=self._plot_color[col], antialias=False, **kwargs)
        self.window.addItem(plt)        
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

        self.draw_direction(px, py, pz, *dirx, col='b', **kwargs)
        self.draw_direction(px, py, pz, *diry, col='g', **kwargs)
        self.draw_direction(px, py, pz, *dirz, col='r', **kwargs)
        return

    def draw_line(self, x0, y0, z0,
                  x01, y01, z01,
                  col='b', **kwargs):

        line = np.array(zip([x0,x01],
                            [y0,y01],
                            [z0,z01]))

        plt = gl.GLLinePlotItem(pos=line, color=self._plot_color[col], antialias=False, **kwargs)
        self.window.addItem(plt)
        return

    def draw_line2(self, v0, v1,
                  col='b', **kwargs):
        x0,y0,z0 = v0
        x01,y01,z01 = v1
        line = np.array(zip([x0,x01],
                            [y0,y01],
                            [z0,z01]))
        if type(col) is str:
            plt = gl.GLLinePlotItem(pos=line, color=self._plot_color[col], antialias=False, **kwargs)
        else:
            plt = gl.GLLinePlotItem(pos=line, color=pg.glColor(*col), antialias=False, **kwargs)
        self.window.addItem(plt)
        return

    def draw_direction(self, x,y,z, dx,dy,dz, **kwargs):
        self.draw_line(x, y, z,
                       x + dx,
                       y + dy,
                       z + dz, **kwargs)
        return

    def add_subplot(self, row, col, index, title=''):
    ###       sudplots seperate window
        index = index - 1
        num = row*col
        tmp = title

        index = index % num
            
        n = index / row
        m = index / col
        ##import pdb; pdb.set_trace()

        p = self.win.addPlot(title=tmp, row=m, col=n)
        self.plots.append(p)


    def plot(self, y_values):
        curve = self.plots[-1].plot(pen='y')
        curve.setData(y_values)

        
    def __initWindow(self):
        w = gl.GLViewWidget()
        w.opts['distance'] = 4
        w.setWindowTitle('pyqtgraph example: GLLinePlotItem')
        w.show()
        return w

    def __setupEnvironment(self,w):
        gz = gl.GLGridItem()
        s = 0.5
        gz.setSize(4,4,4)
        gz.setSpacing(s, s, s)
        gz.translate(0, 0, 0)
        w.addItem(gz)

        self.draw_line(0,0,0,
                     1,0,0, 'b')
        self.draw_line(0,0,0,
                     0,1,0, 'g')
        self.draw_line(0,0,0,
                     0,0,1, 'r')

    def show(self):
        self.__startQt()

    def __startQt(self):
        ## Start Qt event loop unless running in interactive mode.
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()


if __name__ == '__main__':
##    from pylab import *
##    f = figure()
##    f.add_subplot(2,2,0)
##    x = [1,2,3]
##    y = [4,5,6]
##    plot(x,y)
##    show()
    plot = QtPlot()
    plot.show()
