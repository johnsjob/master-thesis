# -*- coding: utf-8 -*-
import serial, os

is_windows = (os.name == 'nt') or (os.name == 'win32')
is_linux = not is_windows

if is_linux:
        STANDARD_PORT = "/dev/ttyACM1"
else:
        STANDARD_PORT = "COM5"

def connect_port(port = STANDARD_PORT, SERIAL_TIMEOUT = 1.0):
                dbg = None
                major,minor = None, None
                print '\n[  connect_port  ]'
                try:
                        dbg = ComDebugChannel(port, SERIAL_TIMEOUT)
                        dbg.connect()
                        if dbg.channel:
                                print "Port["+port+"]: active."
                except Exception as e:
                        if not (str(e) ==  "ComDebugChannel: Could not configure port: (5, 'Input/output error')"):
                                print "Port["+port+"]: unknown exeption catched: \n\t"+str(e)
                if dbg.channel is not None:
                        try:
                                (major, minor) = dbg.get_version()
                        except Exception as e:
                                print "Port["+port+"]: error communicating with port."
                                return None
                        else:
                                print "Port["+port+"]: connected."
                                connected = True
                                print('DebugChannel protocol version: %d.%d' % (major, minor))
                                print(time.asctime())
                        if major is not None:
                                try:
                                        dotpos_ver     = dbg.get_pps_version()
                                        print('DotPos revision: %d' % dotpos_ver)
                                except Exception as e:
                                        print('Error finding dotpos revision: ' + str(e))
                                        #pass
                if dbg is not None:
                        print "Turning on coord-flow..."
                        dbg.dbg_coord_flow_on();                
                        print "Listening..."
                        return dbg
                return None
#--------------------------------------------------------------------------#
def read_apr(dbg):
        data = None
        try:
                data = dbg.get_pps_apr_analysis() #check for packet
                if data is not None: #packet obtained
                        T,frameinfo,apr = data
        except IOError as e: #lost connection temporary or invalid data
                        pass
        except Exception as e: #check for permission denied errors and print them
                        if "Permission" in str(e):
                                print e
                        elif "None" in str(e):
                                raise(Exception("Pen not connected or unresponsive."))
                        else:
                                print "UNSUPPORTED PACKET:"
                                print "\t"+str(e)
        if data is not None:
                return T,frameinfo,apr
        return None
#--------------------------------------------------------------------------#
def get_coords(apr):
        Xmm = apr.intX + apr.fracX/256.0
        Ymm = apr.intY + apr.fracY/256.0
        return mat([Xmm / 39.391, Ymm / 35.906])
#--------------------------------------------------------------------------#
def get_fsr(frameinfo):
        return frameinfo.fsr_val
#====================================#s
from numpy import array as mat
from numpy import linalg
from numpy.linalg import norm
#====================================#s
__file__ = './tmatrix.py'
# Should use the files from submodules
import os, sys
from random import random as rand
filePath, fileName=os.path.split(__file__) #you are here
sys.path.append(os.path.normpath(os.path.join(filePath, '../../../../penorientation/source/')))

L = []
from apranalysis import *
from FrameInfo import *
from debugchannel import *

from RfromT import get_R_from_T as get_orientation
from RfromT import rot_tilt_skew

sys.path.append(os.path.normpath(os.path.join(filePath, 'C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\python\\source\\')))
import helperfunctions as help_funcs
import plane_relative as pl_rel_funcs
from calibration_algorithm3 import solve_tool0_tip as solve_tip
from calibration_algorithm3 import extract_solutions
#--------------------------------------------------------------------------#
def get_T44(plane, anoto2D, rot, tilt, skew, L):
        T44 = numpy.zeros((4,4))
        T44[0:3,0:3] = pl_rel_funcs.get_plane_relative_R(plane, rot, tilt, skew)
        T44[0:3, 3] = pl_rel_funcs.get_plane_relative_point(plane, anoto2D[0], anoto2D[1],
                                                            rot, tilt, skew, L)
        T44[3, 0:4] = [0, 0, 0, 1]
        return T44
#--------------------------------------------------------------------------#
#defining plane and unknwon delta-vector
L = 100 #mm
o = [rand()*3600,rand()*3600,rand()*3600]
R = help_funcs.rotation_matrix_rot_tilt_skew( (rand()*360-180), rand()*90-60, rand()*360-180 )
dirx = R[:,0]
diry = R[:,1]
do = mat([1,2,3])
do = (do / norm(do))*L
plane = help_funcs.define_plane(o, dirx, diry)
#defining T44 orientation
rot, tilt, skew = (45, 40, 30)
#comparing values
comp = mat([dirx, diry, do, do, do])
#--------------------------------------------------------------------------#
pen_values = []
forward_kinematics = []
if __name__ == '__main__':
        print "Arguments(python): " + str(sys.argv)
        if len(sys.argv) < 2: 
                dbg = connect_port()
        else:
                dbg = connect_port(sys.argv[1])
        first_it = True
        origin = None
        while True:
                res = read_apr(dbg)
                if res is not None:
                        T, frameinfo, apr = res
                        T = T / 2**16.0
                        ret = get_coords(apr)
                        if first_it is True:
                                origin = ret
                                first_it = False

                        pen_values.append(ret - origin)
                        forward_kinematics.append( get_T44(plane, pen_values[-1], rot, tilt, skew, L) )
                        print "PENVALUES = " + str(pen_values[-1])
                        #The pen orientation is emulating the orientation of the flange
                        R, debug = get_orientation(T)
                        rot, tilt, skew = rot_tilt_skew(R)
#--------------------------------------------------------------------------#
        log('Exiting...')
