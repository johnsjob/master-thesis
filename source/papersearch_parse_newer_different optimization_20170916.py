from __future__ import division
# imports - core
import os, sys
import logging
import cPickle as pickle
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

import numpy
import numpy as n
numpy.set_printoptions(precision=2)
numpy.set_printoptions(suppress=True)

from numpy import array as mat, pi, cos, sin, arctan2 as atan2,\
                  arccos as acos, log10, degrees as deg, radians as rad
from numpy.linalg import norm, svd, inv

from numpy.random import randint, seed
from copy import deepcopy
import itertools as it
import pylab

# imports - anoto
filePath, fileName=os.path.split(__file__) #you are here
sys.path.append(os.path.normpath(os.path.join(filePath, '../int')))

# imports - master thesis
sys.path.append(os.path.normpath(os.path.join(filePath, '../int/master_thesis_code/source')))
from denavit_hartenberg140 import forward_kinematics, DH_TABLE as dh_table
from helperfunctions_math import homogenous_matrix as hom, matmul,\
                                quat_to_rot, rot_to_quat, nmap, \
                                rot_tilt_skew,\
                                nzip
import calibration_algorithm_newer as calib
from orientation_verification import solve_ori as SVD
from standardplot import StPlot


# imports - newer thesis files
#from ed_probe import PenED, rx_ry_rz as rot_xyz

# imports - 3rd party
import yaml

def normalize(x):
    return x / norm(x)

def tool_svd(A,B):
    points = mat(range(30)).reshape(3,10)
    a = A.dot(points)
    b = B.dot(points)
    ac    = a.T - n.mean(a.T, axis=0)
    bc    = b.T - n.mean(b.T, axis=0)
    h     = mat(map(lambda x: n.outer(x[0],x[1]),zip(ac,bc)))
    H     = n.sum(h, axis=0)
    U,S,V = svd(H)
    return U.dot(V)

class Calibration(object):
    def __init__(self, data_parser, num_points=-1, random=False):
        self.parser = data_parser
        self.__calibrate(num_points, random)

    def __calibrate(self, num_points=-1, random=False):
        self.robot_infos = map(lambda x: forward_kinematics(*x, **dh_table), self.parser.robot_joints)
        fk = self.flange

        # anoto coords
        pentip_2d = self.parser.pen_pos

        # convert to mm, first point is origin
        #pentip_2d = (pentip_2d - pentip_2d[0])*0.3
        f_ad = 7.0 *(25.4/600.0)
        pentip_2d = (pentip_2d) * f_ad
        if not num_points == -1:
            if random:
                unique_rand = list(set(randint(0, len(fk), num_points)))
                fk = fk[unique_rand]
                pentip_2d = pentip_2d[unique_rand]
        self.geometry_info = {
            'data':
                {
                'forward_kinematics':fk,
                'pentip_2d':pentip_2d
                }
            }
        self._res, self.cond = calib.find_solution_pen_tip(self.geometry_info, num_points)

    @property
    def wobj_ori(self):
        return self._res[:,:3]
    @property
    def wobj_ori_x(self):
        return self.wobj_ori[:,0]
    @property
    def wobj_ori_y(self):
        return self.wobj_ori[:,1]
    @property
    def wobj_ori_z(self):
        return self.wobj_ori[:,2]

    @property
    def wobj_error(self):
        return abs(90 - acos(self.wobj_ori[:,0].dot(self.wobj_ori[:,1]))*180.0/pi)

    @property
    def wobj_error_x(self):
        if (self.parser.yaml_wobj is None):
            return numpy.NaN
        acos_val = acos(normalize(self.wobj_ori_x).dot(self.parser.yaml_wobj_ori_x))
        return deg(acos_val)

    @property
    def wobj_error_y(self):
        if (self.parser.yaml_wobj is None):
            return numpy.NaN
        return acos(normalize(self.wobj_ori_y).dot(self.parser.yaml_wobj_ori_y))*180.0/pi

    @property
    def wobj_error_z(self):
        if (self.parser.yaml_wobj is None):
            return numpy.NaN
        return acos(normalize(self.wobj_ori_z).dot(self.parser.yaml_wobj_ori_z))*180.0/pi

    @property
    def tip(self):
        return self._res[:,3]

    @property
    def tip_error(self):
        """
        Used for optimizing the calibration result.
        Calculates the tip error using a metric function.
        """
        return self.global_tip_error
#        if (self.parser.yaml_tip is None):
#            return  numpy.NaN
#        return norm(self._res[:,3] - self.parser.yaml_tip)

    def tip_error_max(self):
        return self.global_tip_error

    @property
    def global_tip_error(self):
        print self.parser.num_points
        print len(self.parser.pen_pos_pattern_wobj_global)
        print len(self.tcp)
        k = min(len(self.tcp), self.parser.num_points)
        print "k = {}".format(k)
        print "..."
        global_error = abs(self.tcp[:k] - self.parser.pen_pos_pattern_wobj_global[:k])[:,:,3]
        max_xyz = norm(n.max(global_error,0)[:3])
        mean_xyz = norm(n.mean(global_error,0)[:3])
        min_xyz = norm(n.min(global_error,0)[:3])
        err_metric = (max_xyz, mean_xyz, min_xyz)        
        return sum(err_metric)

    @property
    def tool(self):
        return hom(0,0,0,self.tip)

    @property
    def tcp(self):
        tcp  = matmul(
                self.flange,
                self.tool
                )
        return tcp

    @property
    def flange(self):
        flange  = mat([x['flange'] for x in self.robot_infos])
        flange[:,:3,3] = flange[:,:3,3]*1000 #convert to mm
        return flange

    @property
    def flange_pos(self):
        flange_pos  = self.flange[:,:3,3]
        return flange_pos

    @property
    def flange_pos_error(self):
        return nmap(norm, self.flange_pos - self.parser.robot_flange)

    @property
    def wobj(self):
        return hom(self._res[:,:3])

    @property
    def robot_frames_global(self):
        frames = mat([x['robot_geometry_global'] for x in self.robot_infos])
        frames[:,:,:3,3] = frames[:,:,:3,3] * 1000.0
        return frames
    
class DataParser(object):

    # constructor
    def __init__(self, new_data, clean_data=False):
        print "Initializing data parser..."
        if type(new_data) == str:
            with open(new_data) as fp:
                self.data = pickle.load(fp)
        elif type(new_data) == list or type(new_data) == tuple:
            self.data = new_data
        else:
            self.data = None

        if clean_data:
            self.__clean()

        self._yaml_data = None
        self._yaml_file = None
        self._subparse_len = 20

    def __clean(self):
        tmp = []
        all_pos = self.pen_pos
        for index, pos in enumerate(all_pos):
            if not pos.dtype == numpy.object:
                tmp.append(index)
        self.data = [x[1] for x in enumerate(self.data) if (x[0] in tmp)]

    def to_pickle(self, filepath):
        with open(filepath,'wb+') as fp:
            pickle.dump(self.data, fp)

    def merge(self, other_parser):
        self.data += other_parser.data
        return self

    def emptycopy(self):
        self.subparser = 0
        return self.subparser

    def jumblecopy(self):
        self.subparser = self.numpoints
        return self.subparser

    @property
    def subparser(self):
        tmp = deepcopy(self.data)
        points = [self.data.pop(randint(0, self.num_points-1)) for _ in xrange(self._subparse_len)]
        self.data = tmp
        new_parser = DataParser(points)
        new_parser.yaml_file = self.yaml_file
        return new_parser

    @subparser.setter
    def subparser(self,value):
        self._subparse_len = value

    @property
    def num_points(self):
        return len(self.data)

    @property
    def robot_joints(self):
        return mat([x['robot_joints'] for x in self.data])

    @property
    def robot_flange(self):
        return mat([x['robot_flange'] for x in self.data])

    @property
    def robot_tcp(self):
        return mat([x['robot_tcp'] for x in self.data])

    @property
    def pen_pos(self):
        return mat([x['pen_pos'] for x in self.data])

    def gen_robot_flanges(self):
        robot_infos = map(lambda x: forward_kinematics(*x, **dh_table), self.robot_joints)
        flange  = mat([x['flange'] for x in robot_infos])
        flange[:,:3,3] = flange[:,:3,3]*1000 #convert to mm
        return flange

    @property
    def pen_pos_pattern_mm(self):
        """
        Pen position relative to paper origin, which
        is the UL position of the pattern with crosses
        used in calibration. The coordinates are converted from ad to mm.
        """
        ad_to_mm = 7.0 *(25.4 / 600.0)
        paper_origin = mat([385911998.7,327159331.8])
        return mat([x['pen_pos'] - paper_origin \
                    for x in self.data]) * ad_to_mm
    @property
    def pen_pos_pattern_mm_3d(self):
        """
        Pen position relative to paper origin, which
        is the UL position of the pattern with crosses
        used in calibration. The coordinates are converted from ad to mm.
        """
        return nmap(lambda x: x + [0], self.pen_pos_pattern_mm.tolist())

    @property
    def pen_pos_pattern_mm_hom(self):
        return nmap(lambda x: hom(0,0,0,x), self.pen_pos_pattern_mm_3d)

    @property
    def pen_pos_pattern_wobj_global(self):
        return nmap(lambda x: self.yaml_wobj.dot(x), self.pen_pos_pattern_mm_hom)

    @property
    def pen_ori(self):
        return mat([x['pen_ori'] for x in self.data])

    @property
    def pen_angles(self):
        return mat([rot_tilt_skew(x['pen_ori']) for x in self.data])

    @property
    def pen_fsr(self):
        return mat([x['pen_fsr'] for x in self.data])

    @property
    def pen_fsr_adc(self):
        return mat([x['pen_fsr_adc'] for x in self.data])

    @property
    def yaml_tip(self):
        if not self.__yaml_data:
            return None
        return mat(self.__yaml_data['tool']['pos'])

    @property
    def yaml_tool_q(self):
        if not self.__yaml_data:
            return None
        return self.__yaml_data['tool']['q']

    @property
    def yaml_wobj(self):
        """
        The coordinate system of the paper pattern.
        """
        if not self.__yaml_data:
            return None
        return hom(self.yaml_wobj_ori, self.yaml_wobj_pos)

    @property
    def yaml_wobj_ori(self):
        """
        The coordinate system of the paper pattern.
        This function converts the orientation from global target
        coordiante system to global pattern oordinate sytem
        by the transformation: [ 0  1  0
                                 1  0  0
                                 0  0 -1 ]
        """
        if not self.__yaml_data:
            return None
        return matmul(quat_to_rot(self.__yaml_data['wobj']['q']),
                      mat([[0,1,0],
                           [1,0,0],
                           [0,0,-1]]))
    @property
    def yaml_wobj_q(self):
        return self.__yaml_data['wobj']['q']
    
    @property
    def yaml_wobj_ori_x(self):
        if not self.__yaml_data:
            return None
        return self.yaml_wobj_ori[:,0]

    @property
    def yaml_wobj_ori_y(self):
        if not self.__yaml_data:
            return None
        return self.yaml_wobj_ori[:,1]

    @property
    def yaml_wobj_ori_z(self):
        if not self.__yaml_data:
            return None
        return self.yaml_wobj_ori[:,2]

    @property
    def yaml_wobj_pos(self):
        if not self.__yaml_data:
            return None
        return mat(self.__yaml_data['wobj']['pos'])

    @property
    def yaml_file(self):
        if self.__yaml_data:
            return self._yaml_file
        else:
            return None

    @property
    def yaml_dir(self):
        if self.__yaml_data:
            return os.path.join(os.path.dirname(self.yaml_file))
        else:
            return None

    @yaml_file.setter
    def yaml_file(self,filepath):
        self.__yaml_data = filepath
        if self.__yaml_data:
            self._yaml_file = filepath
        else:
            self._yaml_file = None

    @property
    def __yaml_data(self):
        return self._yaml_data

    @__yaml_data.setter
    def __yaml_data(self, filepath):
        try:
            with open(filepath) as fp:
                self._yaml_data = yaml.load(fp)
        except Exception as e:
            self._yaml_data = None
            log.warning(str(e))

    # destructor
    def __del__(self):
        pass

def parser_merge(x,y):
    x.merge(y)
    return x

def create_parser(x):
    return DataParser(x)

def calibration_optimization_rand(p, error):
    curr_error = error
    curr_calib = Calibration(p)
    for k in range(p.num_points*50):
        try:
            c = Calibration(p.subparser)
        except:
            continue
        if c.tip_error < curr_error:
            curr_calib = deepcopy(c)
            curr_error = deepcopy(c.tip_error)
            print curr_error
    return curr_calib.parser, error

def try_remove_one(parser):
    c = Calibration(parser)
    for k in range(parser.num_points):
        tmp = parser.data.pop(k)
        new_c = Calibration(parser)
        if new_c.tip_error < c.tip_error:
            return tmp
        else:
            parser.data.insert(k, tmp)

def try_remove_largest(parser):
    c = Calibration(parser)
    curr_k = None
    curr_error = c.tip_error
    for k in range(parser.num_points):
        tmp = parser.data.pop(k)
        new_c = Calibration(parser)
        if new_c.tip_error <= curr_error:
            curr_k = k
            curr_error = new_c.tip_error
        parser.data.insert(k, tmp)
    if not (curr_k is None):
        return parser.data.pop(curr_k)

def calibration_optimization(parser):
    curr_parser = parser.emptycopy()
    curr_parser.data = [parser.data.pop(0) for _ in xrange(6)]
    curr_calib = Calibration(curr_parser)
    counter = 0
    limit = parser.num_points
    print '{} | {}'.format(curr_calib.tip_error,parser.num_points)
    while parser.num_points > 0:
        curr_parser.data.append(parser.data.pop(0))
        new_calib = Calibration(curr_parser)
        if new_calib.tip_error < curr_calib.tip_error:
            curr_calib = new_calib
            counter = 0
            print '{} | {}'.format(curr_calib.tip_error,parser.num_points)
        else:
            ret = try_remove_one(curr_parser)
            #ret = try_remove_largest(curr_parser)
            if not (ret is None) and (counter < limit):
                parser.data.append(ret)
            else:
                print '{} | {}'.format(curr_calib.tip_error,parser.num_points)
                counter=0
        counter += 1

    print curr_calib.tip_error
    while try_remove_largest(curr_parser):
        pass
    return curr_parser

if __name__ == '__main__':
    _debug = False
    _optimize = False
    if not len(sys.argv) > 1:
        log.error('\n\t\tUsage: please supply file!')
        #sys.argv.append( raw_input('Filepath: ') )
        ##sys.argv.append( 'C:\\Users\\***REMOVED***\\Dropbox\\exjobb\\results\\measurements\\measurements_automatic\\26_27052016\\4_measurement_27052016_0227.pickle' )
        wobj_origin = mat([385911998.7,327159331.8])
        o = wobj_origin
        tmp_path = r"C:\Users\***REMOVED***\Desktop\final_meas_20160724_2016025\all\measurement_24072016_2251_2321_abort.pickle"
        sys.argv.append(tmp_path)
        _debug = True

    print sys.argv
    
    if len(sys.argv) == 2:
        parser = DataParser(sys.argv[1])
    else:
        parsers = map(create_parser, sys.argv[1:])
        parser = reduce(parser_merge, parsers)

    yaml_filepath = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])),
                                 'tool_wobj.yaml')
    parser.yaml_file = yaml_filepath

    if raw_input('Optimize results?: ').lower() == 'y':
        _optimize = True
        print '(be patient, this will take some time ...)'
        org_parser = deepcopy(parser)
        print parser.num_points
        cal = Calibration(parser)
        print cal.tip_error
        curr_parser = calibration_optimization(parser)
        print "===[ Calibration Result ]==="
        print curr_parser.num_points
        print "tip error metric: {}".format(Calibration(curr_parser).tip_error)
        #pos = matmul(org_parser.gen_robot_flanges(), cal.tool)[:,:3,3]
        #meas_pos = org_parser.pen_pos_pattern_wobj_global[:,:3,3]
        #res = meas_pos - pos
        parser = curr_parser
        yaml_output = {
            "tool": {
                "pos":cal.tip.tolist(),
                "q": org_parser.yaml_tool_q,
                "name": "Optimized calibration tip"
                },
            "wobj": {
                "pos": org_parser.yaml_wobj_pos.tolist(),
                "q": org_parser.yaml_wobj_q
                }
            }
        with open(os.path.join(parser.yaml_dir, "tool_wobj_opt.yaml"), 'w+') as fp:
            yaml.dump(yaml_output, fp)
        
        if not _debug:
            parser.to_pickle(os.path.join(parser.yaml_dir,'optimized_measurements.pickle'))

    # plotting - misc
    pos = parser.pen_pos
    pos = pos - pos[0]
    pylab.grid()
    pylab.plot(pos[:,0],
               pos[:,1],'k.')
    pylab.xlabel('x [mm]')
    pylab.ylabel('y [mm]')
    if not _debug:
        filename = 'pen_coords'
        if _optimize:
            filename += '_opt'
        pylab.savefig(os.path.join(parser.yaml_dir,'{}.png'.format(filename)))
        pylab.clf()
    else:
        pylab.show()
    calibration_run = []
    x_val = range(2, parser.num_points)
    for p in x_val:
        cal = Calibration(parser, num_points=p)
        calibration_run.append(cal)
    calibration = Calibration(parser)

    from pandas import DataFrame

    df = DataFrame(data=zip(
                            [x.tip_error for x in calibration_run],
                            [x.wobj_error_x for x in calibration_run],
                            [x.wobj_error_y for x in calibration_run],
                            [x.wobj_error_z for x in calibration_run],
                            calibration.flange_pos_error,
                            ),
                   columns=['tip_error','wobj_error_x',
                            'wobj_error_y',
                            'wobj_error_z',
                            'flange_error'],
                   index=x_val)
    print df
    if not _debug:
        filename = 'calibration_errors'
        if _optimize:
            filename += '_opt'
        df.to_csv(os.path.join(parser.yaml_dir,'{}.csv'.format(filename)))
        df.to_latex(os.path.join(parser.yaml_dir,'{}.tex'.format(filename)))

    tip_error = log10(calibration.tip_error)
    pylab.plot(x_val, log10([x.tip_error for x in calibration_run]))
    pylab.xlim([x_val[0], x_val[-1]])
    pylab.axhline(tip_error, color='k', linestyle='--')
    pylab.text(parser.num_points - 3,0.08+tip_error,
               '${:0.4f}$ mm'.format(calibration.tip_error),
               bbox={'facecolor':'white',
                     'edgecolor':'white'})
    pylab.xlabel('number of measurements')
    pylab.ylabel('log10 norm error [mm]')
    pylab.title('Tool tip calibration results from automatic measurements')
    pylab.grid()
    if not _debug:
        filename = 'tip_error'
        if _optimize:
            filename += '_opt'
        pylab.savefig(os.path.join(parser.yaml_dir,'{}.png'.format(filename)))
        pylab.clf()
    else:
        pylab.show()

    pylab.plot(x_val, log10([x.wobj_error_x for x in calibration_run]))
    pylab.plot(x_val, log10([x.wobj_error_y for x in calibration_run]))
    pylab.plot(x_val, log10([x.wobj_error_z for x in calibration_run]))
    pylab.xlim([x_val[0], x_val[-1]])
    pylab.legend(['$\hat{x}_{an}$',
                  '$\hat{y}_{an}$',
                  '$\hat{z}_{an}$'])
    pylab.xlabel('number of measurements')
    pylab.ylabel('log10 angle error [deg]')
    pylab.title('Work object calibration results from automatic measurements')
    pylab.grid()
    if not _debug:
        filename = 'wobj_error'
        if _optimize:
            filename += '_opt'
        pylab.savefig(os.path.join(parser.yaml_dir,'{}.png'.format(filename)))
        pylab.clf()
    else:
        pylab.show()

##    # plotting - result
##    for k in range(1):
##        index = k
##        plotting = StPlot(length_unit='mm')
##        plotting.draw_robot(calibration.robot_frames_global[index])
##        plotting.draw_frame(calibration.wobj[0], size=0.4)
##        try:
##            plotting.draw_frame(parser.yaml_wobj, size=0.5)
##        except:
##            pass
##        plotting.draw_tool(calibration.flange[index], calibration.tool)
##    plotting.show()
