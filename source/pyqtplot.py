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
                        'gr': pg.glColor(128,128,128)
                      } 

        self.app = QtGui.QApplication([])
        self.window = self.__initWindow()

        self.__setupEnvironment(self.window)



###       sudplots seperate window
##        win = pg.GraphicsWindow(title="Basic plotting examples")
##        win.resize(1000,600)
##        win.setWindowTitle('pyqtgraph example: Plotting')
##
##        # Enable antialiasing for prettier plots
##        pg.setConfigOptions(antialias=True)
##
##        p1 = win.addPlot(title="Basic array plotting", y=np.random.normal(size=100), pen=(0,255,0))
##
##        p2 = win.addPlot(title="Multiple curves")
##        p2.plot(np.random.normal(size=100), pen=(255,0,0), name="Red curve")
##        win.nextRow()
##
##        p3 = win.addPlot(title="Drawing with points")
##        p3.plot(np.random.normal(size=100), pen=(200,200,200), symbolBrush=(255,0,0), symbolPen='w')
##        p3 = win.addPlot(title="Drawing with points")
##        p3.plot(np.random.normal(size=100), pen=(200,200,200), symbolBrush=(255,0,0), symbolPen='w')


        return
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

    def draw_curve(self, points, col='g', **kwargs):
        plt = gl.GLLinePlotItem(pos=points, color=self._plot_color[col], antialias=False, **kwargs)
        self.window.addItem(plt)        
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

    def draw_frame(self, frame_matrix, size=1, **kwargs):
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
    plot = QtPlot()
    plot.show()
