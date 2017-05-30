# imports - core
import os, sys
import logging
import ast
import cPickle as pickle
from random import random as rand
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
import time

from numpy import array as mat, pi, cos, sin, \
                  arctan2 as atan2, linspace, meshgrid as mesh
from numpy.linalg import norm
from numpy.random import uniform, seed
import itertools as it

# imports - anoto
filePath, fileName=os.path.split(__file__) #you are here
sys.path.append(os.path.normpath(os.path.join(filePath, '../int')))

import robot_com
from robot_com.ftp.robotFTP import robotFileTransfer
from robot_com.serial.arap import robotSerial

# imports - master thesis
sys.path.append(os.path.normpath(os.path.join(filePath, '../int/master_thesis_code/source')))
from denavit_hartenberg140 import forward_kinematics, DH_TABLE as irb140_params
from helperfunctions_math import rotation_matrix_rot_tilt_skew as R_rts, matmul
from standardplot import StPlot

from threading import Thread, Lock, RLock

# imports - newer thesis files
from ed_probe import PenED, rx_ry_rz as rot_xyz

class Robot:

    # inits
    def _upload_files(self):
        ftp = robotFileTransfer(debug_state = debug)
        ftp.connect()
        ftp.upload('ARAPMOD.prg', robot_com.paths['serial'])
        ftp.upload('hptr_control.mod', robot_com.paths['server'])
        ftp.upload('common.mod', robot_com.paths['server'])
        ftp.upload('calibrationdata.mod', robot_com.paths['server'])
        ftp.disconnect()
        self._ftp = ftp

    def _init_serial(self):
        pass
        arap = robotSerial(debug_state = debug)
        self._arap = arap
        arap.load('hptr_control.mod')
        arap.load('hptr_control.mod') # perform 2 times if pen is inited before robot
        arap.load('common.mod')
        arap.load('calibrationdata.mod')
        arap.command('initHaptor')
        self.move_to_door()

    # constructor
    def __init__(self,lock = None, angle=45, pen_interface=None, num_data=16):
        self._upload_files()
        self._init_serial()

        self.lock = lock
        self.all_data = []

        self.start_time = time.strftime('%H%M')

        self.alive = False
        self._finished = False
        self.num_data_points = num_data
        if pen_interface:
            self.pen = pen_interface
            self.pen_hit_thread = Thread(target=pen.check_hit)
            self.pen_hit_thread.start()
            time.sleep(1)
            with self.lock:
                self.alive = self.pen.alive
            if pen.alive:
                self.move_to_ready(angle)
                self.save_tool_pos()

    # destructor
    def __del__(self):
        print('Robot shutdown...')
        if self.all_data:
            print('Saving data...')
            if self._finished:
                filename = "measurement_%d%m%Y_{}_%H%M.pickle".format(self.start_time)
            else:
                filename = "measurement_%d%m%Y_{}_%H%M_abort.pickle".format(self.start_time)

            # remove the first coordinate,
            # since it will become dupliate coordinate
            # during multiple runs, pen angles are also ill-conditioned for rot/skew
            self.all_data.pop(0)
            with open(time.strftime(filename),'wb+') as fp:
                pickle.dump(self.all_data, fp)
        else:
            print('No data to save!')

        # prepare shutdown
        self.move_to_door()
        self.pen_hit_thread.join()

        # call dels
        del self._arap
        self._arap = None

        del self._ftp
        self._tp = None
        print('Shutdown complete!')
        return

    # properties
    @property
    def arap(self):
        return self._arap

    # functions
    def __abort(self):
        if not self.pen:
            return
        pen._abort = True

    def move_to_door(self):
        self.arap.command('changePen')
        self.save_tool_pos()

    def set_vel_molusk(self):
        if not self.alive:
            return
        self.set_vel(0.01)
        self.set_ang_vel(30)

    def set_vel_parrot(self):
        if not self.alive:
            return
        self.set_vel(60)
        self.set_ang_vel(30)

    def set_vel(self, vel):
        if not self.alive:
            return
        self.arap.command('setVel', vel)

    def set_ang_vel(self, vel):
        if not self.alive:
            return
        self.arap.command('setAngVel', vel)

    def set_rel_tcp_z(self, z):
        if not self.alive:
            return
        self.arap.command('SetRelTCPz', '{:0.3f}'.format(z))

    def save_tool_pos(self):
        self.arap.command('SavePos')

    def move_to_saved_pos(self):
        if not self.alive:
            return
        self.arap.command('MoveSavedPos')

    def move_to_saved_pos_ori(self):
        if not self.alive:
            return
        self.arap.command('MoveSavedPosOri')

    def move_to_saved_ori(self):
        if not self.alive:
            return
        self.arap.command('MoveSavedOri')

    def move_tcp(self, x,y,z, speed):
        if not self.alive:
            return
        if speed == 'normal':
            self.arap.command('MoveTCPNormal', '[{0:0.3f},{1:0.3f},{2:0.3f}]'.format(x,y,z))
        elif speed == 'slow':
            self.arap.command('MoveTCPSlow', '[{0:0.3f},{1:0.3f},{2:0.3f}]'.format(x,y,z))

    def move_wobj_tcp(self, x,y,z):
        if not self.alive:
            return
        self.arap.command('MoveWObjTCP', '[{0:0.3f},{1:0.3f},{2:0.3f}]'.format(x,y,z))

    def move_flange(self, x,y,z):
        if not self.alive:
            return
        self.arap.command('MoveFlange', '[{0:0.3f},{1:0.3f},{2:0.3f}]'.format(x,y,z))

    def rel_tool_ori(self, rx, ry, rz):
        if not self.alive:
            return
        self.arap.command('RelToolOri', '[{0:0.2f},{1:0.2f},{2:0.2f}]'.format(rx, ry, rz))

    def rel_tool_dir(self, x,y,z):
        if not self.alive:
            return
        self.arap.command('RelToolPos', '[{0:0.3f},{1:0.3f},{2:0.3f}]'.format(x,y,z))

    def rel_tool_z(self, z):
        if not self.alive:
            return
        self.arap.command('RelZPos', '{:0.3f}'.format(z))

    def abs_tool_dir(self, x,y,z):
        if not self.alive:
            return
        self.arap.command('RelBasePos', '[{0:0.3f},{1:0.3f},{2:0.3f}]'.format(x,y,z))

    def move_to_ready(self, angle):
        if not self.alive:
            return
        self.arap.command('SetJ123', [0,0,0])
        self.arap.command('SetJzero')
        self.arap.command('SetJ5', angle)
        print self.get_tcp()

    def get_flange(self):
        while True:
            try:
                res = self.arap.command('getFlange')
                i0 = res.find('[')
                i1 = res.find(']')+1
                res = res[i0 : i1]
                res = ast.literal_eval(res)
                if not(len(res) == 3):
                    raise Exception('Not valid FLANGE format (length: {}'.format(len(res)))
                return mat(res)
            except Exception as e:
                print 'Failed to obtain Flange-values!'
                print str(e)

    def get_tcp(self):
        while True:
            try:
                res = self.arap.command('getTCP')
                i0 = res.find('[')
                i1 = res.find(']')+1
                res = res[i0 : i1]
                res = ast.literal_eval(res)
                if not(len(res) == 3):
                    raise Exception('Not valid TCP format (length: {}'.format(len(res)))
                return mat(res)
            except Exception as e:
                print 'Failed to obtain TCP-values!'
                print str(e)


    def get_joints(self):
        while True:
            try:
                res = self.arap.command('getJ')
                i0 = res.find('[')
                i1 = res.find(']')+1
                res = res[i0 : i1]
                res = ast.literal_eval(res)
                if not(len(res) == 6):
                    raise Exception('Not valid joints format (length: {}'.format(len(res)))
                return mat(res)
            except Exception as e:
                print 'Failed to obtain joint-values!'
                print str(e)

    def __grab_data(self):
        with self.lock:
            #these might need a few tries to get a valid result
            J = self.get_joints()
            tcp = self.get_tcp()
            flange = self.get_flange()

            pos = self.pen.position
            if pos is None:
                return None
            pos = mat(pos)

            ori = self.pen.orientation
            if ori is None:
                return None
            ori = mat(ori)

            fsr = self.pen.fsr
            if fsr is None:
                return None

            fsr_adc = self.pen.fsr_adc
            if fsr_adc is None:
                return None

            perspective = self.pen.perspective
            if perspective is None:
                return None

            data = {
            'pen_pos': pos,
            'pen_ori': ori,
            'pen_fsr': fsr,
            'pen_fsr_adc': fsr_adc,
            'pen_perspective': perspective,
            'robot_joints': J,
            'robot_tcp': tcp,
            'robot_flange': flange
            }
            return data

    def __find_start_pos(self, amount, height):
        if (not self.pen) or (not self.arap):
            return
        elif self.pen:
            if not self.pen.alive:
                return

        log.info('preparing SEARCH MODE.')
        self.save_tool_pos()
        while True:
            if not self.pen.alive:
                log.error('LOST CONNCTION TO PEN, ABORTING.')
                break
            if self.pen.hit:
                self.tool_height = height
                self.rel_tool_dir(0, 0, height) # move back 10mm
                self.save_tool_pos()
                log.info('Found starting point!')
                return
            else:
                self.arap.command('DispDistSaved')
                self.set_vel(10)
                self.rel_tool_dir(0,0, -amount) # move forward by 'amount'
        return

    def __search_tool_z(self, amount, _grid, _angles, _start_index):
        if (not self.pen) or (not self.arap):
            return
        elif self.pen:
            if not self.pen.alive:
                return

        if (_start_index is None):
            index = 0
        else:
            index = _start_index

        if (_grid is None):
            disp_pos = mat([[10,0],[0,10],[-10,0],[0,-10]])
        else:
            disp_pos = _grid

        if (_angles is None):
            angles = [[uniform(-20,20), uniform(-20,20), uniform(-180, 180)]for _ in xrange(len(disp_pos))]
        else:
            angles = _angles
        log.info('SEARCH MODE.')
        print 'Number of angles: {}'.format(len(angles))

        self.save_tool_pos()
        self._finished = True # This remains if everything works

        assert(len(angles) == len(disp_pos))

        while not (len(self.all_data) == len(angles)-_start_index):
            if not self.pen.alive:
                print 'CURRENT INDEX: {}'.format(index)
                log.error('LOST CONNCTION TO PEN, ABORTING.')
                self._finished = False
                break

            if self.pen.hit:
                data = self.__grab_data()
                if data is None:
                        print 'FAILED MEASUREMENT - TRYING AGAIN!'
                        continue
                self.all_data.append(data)
                print 'HIT!\a'

                print 'CURRENT INDEX: {}'.format(index)
                for key in data:
                    print '\n{}: {}'.format(key, data[key])
                self.set_vel(60)

                self.move_to_saved_pos()
                self.move_to_saved_ori()

                pos = disp_pos[index % len(disp_pos)]
                ang = angles[index % len(angles)]
                if (_grid is None):
                    if (index+1) % 4 == 0:
                            perturb = mat((rand()*10-5, rand()*10-5))
                            pos = pos + perturb

                if (_grid is None):
                    self.abs_tool_dir(pos[0],pos[1], 0)
                else:
                    self.move_wobj_tcp(pos[0], pos[1], self.tool_height)

                self.save_tool_pos()
                self.rel_tool_ori(*ang)
                index = index+1
            else:
                self.arap.command('DispDistSaved')
                self.set_vel(10)
                self.rel_tool_dir(0,0, -amount)
        self.move_to_saved_pos()
        self.__abort()
        return

    def __rel_toolz_ori(self, rel_z, rx,ry,rz):
        self.set_rel_tcp_z(rel_z)
        self.rel_tool_ori(rx, ry, rz)
        self.set_rel_tcp_z(-rel_z)

    def paper_verify(self):
        self.__find_start_pos(0.1, height=30)
        self.__rel_toolz_ori(30, 0, -20, 0)
        self.rel_tool_z(10)
        self.rel_tool_z(10)
        self.move_to_saved_pos_ori()
        self.__abort()


    def paper_search(self, grid=None, angles=None,
                           start_index=None, start_height=None):
        self.__find_start_pos(0.1, height=start_height)

        self.__search_tool_z(0.1, _grid = grid, _angles = angles,
                                  _start_index = start_index)
        return

    def parrot(self, prints=False):
        if (not self.pen) or not(self.arap):
            return
        elif self.pen:
            if not self.pen.alive:
                return

        log.info('PARROT MODE.')
        while True:
            if self.pen.hit:
                with self.lock:
                    mxy = mat((353455349.28125, 363.38671875))
                    pos = (self.pen.pos-mxy)*0.3
                L = norm(pos)
                x,y = pos
                if y < -74:
                    self.abs_tool_dir(0, 0, -x)
                else:
                    self.abs_tool_dir(x, -y, 0)

                if prints:
                    J = self.get_joints()
                    tcp = self.get_tcp()
                    flange = self.get_flange()
                    print 'Joints: {}'.format(J)
                    print 'TCP: {}'.format(tcp)
                    print 'Flange: {}'.format(flange)
            else:
                pass
        return

def position_grid(xmin, xmax, ymin, ymax, num):
    x = linspace(xmin,xmax,num)
    y = linspace(ymin,ymax,num)
    xy = mat(mesh(x,y)).T
    #shape: x_grid, y_grid, grid_index/layer
    return xy

def convert_grid(grid):
    a,b,c = grid.shape
    xy = grid.reshape(a*b,c)
    return xy

def define_angles(tilt, skew):
    return [
        [tilt,0,90+skew], [-tilt,0,90+skew], [0,tilt,90+skew],  [0,-tilt,90+skew],
        [tilt,0,-90-skew],[-tilt,0,-90-skew],[0,tilt,-90-skew], [0,-tilt,-90-skew],
        [tilt,0,45+skew], [-tilt,0,45+skew], [0,tilt,45+skew],  [0,-tilt,45+skew],
        [tilt,0,45-180+skew], [-tilt,0,45-180+skew], [0,tilt,45-180+skew],  [0,-tilt,45-180+skew],
        [tilt,0,-45+180-skew],[-tilt,0,-45+180-skew],[0,tilt,-45+180-skew], [0,-tilt,-45+180-skew]
        ]

if __name__ == '__main__':
    seed(123456789*3)
    debug = False
    shared_lock = RLock()
    pen = PenED(hover_mode=True,lock = shared_lock) #Put after robot inits, or perform LOAD 2 times the first time
    robot = Robot(angle=90, pen_interface=pen, lock = shared_lock)

    ## generate grids
    f = 7.0 *(25.4/600.0)
    gridmm = position_grid(0, f*(4572.2 - 3619.8), #282.22786666666656
                           0, f*(3917.4 - 3262.7), #194.00943333333342
                           11)
    gridad = position_grid(3619.8, 4572.2,
                           3262.7, 3917.4,
                           11)

    ## extract subgrids and convert from (a,b,c) shape to (a*b,c) shape
#    gridmm = gridmm[1:-1, 2:-3, :]
#    gridad = gridad[1:-1, 2:-3, :]
    gridmm = convert_grid(gridmm)
    gridad = convert_grid(gridad)
#    gridmm = gridmm[1::2]
#    gridad = gridad[1::2]

    ## comment/uncomment these lines
    ## if we want distorted grids
#    gridmm_pert = mat([uniform(-5,5) for _ in xrange(len(gridmm)*2)]).reshape(len(gridmm),2)
#    gridad_pert = gridmm_pert/f
#    gridmm = gridmm + gridmm_pert
#    gridad = gridad + gridad_pert

    angles = None
    ## comment/uncomment these lines
    ## if we want symmtrical angles
#    num_chunks = len(define_angles(0,0))
#    num_angle_chunks = int(len(gridmm) / num_chunks)
#    angles = mat([define_angles(uniform(-20,20), uniform(-90,90)) for _ in xrange(num_angle_chunks+1)])
#    a,b,c = angles.shape
#    angles = angles.reshape(a*b,c)
#    angles = angles[:len(gridmm)]
#    print '\nAngles:\n{}'.format(angles)

    ## uncommnt this for
    ## checking calibration
    ## These were used for measuring accuracy of wobj calibration
#    robot.rel_tool_ori(20,0,-45) #xtilt
#    robot.rel_tool_ori(0,20,0) #ytilt
#    robot.rel_tool_ori(20,20,-45) #xytilt
    robot.move_wobj_tcp(*(list(gridmm[-1])+[-0.6]))
    while True:
        pass

    ## uncommnt this for
    ## checking grid placement
#    for p in gridmm:
#        pos = list(p) + [5];
#        robot.move_wobj_tcp(*pos)
#    1/0

    # choose starting index
    index=0
    print 'pos mm: {}'.format(gridmm[index])
    print 'pos ad: {}'.format(gridad[index])
    # start measuring

    robot.move_wobj_tcp(*(list(gridmm[index])+[3]))
    robot.paper_search(start_height = 3, start_index = index,
                       grid = gridmm,    angles = angles)

