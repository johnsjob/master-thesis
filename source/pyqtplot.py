# -*- coding: utf-8 -*-
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import pyqtgraph as pg

import numpy as np
from numpy import array as mat
from numpy import linalg as li


class QtPlot:
    def __init__(self):
        self._plot_color = {
                        'b': pg.glColor(0,0,255),
                        'g': pg.glColor(0,255,0),
                        'r': pg.glColor(255,0,0)
                      } 
        self.__doQt()
        return
        
    def draw_line(self, x0, y0, z0,
                  x01, y01, z01,
                  col='b'):
        line = np.array(zip([x0,x01],
                            [y0,y01],
                            [z0,z01]))
        plt = gl.GLLinePlotItem(pos=line, color=self._plot_color[col], antialias=False)
        self.window.addItem(plt)
        return
        
    def __initWindow(self):
        w = gl.GLViewWidget()
        w.opts['distance'] = 40
        w.show()
        w.setWindowTitle('pyqtgraph example: GLLinePlotItem')
        return w

    def __setupEnvironment(self,w):
        gz = gl.GLGridItem()
        gz.translate(0, 0, 0)
        w.addItem(gz)
        self.draw_line(0,0,0,
                     10,0,0, 'b')
        self.draw_line(0,0,0,
                     0,10,0, 'g')
        self.draw_line(0,0,0,
                     0,0,10, 'r')

    def __startQt(self):
        ## Start Qt event loop unless running in interactive mode.
        import sys
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()

    def __doQt(self):
        self.app = QtGui.QApplication([])
        self.window = self.__initWindow()
        self.__setupEnvironment(self.window)
        self.__startQt()


if __name__ == '__main__':
    plot = QtPlot()
