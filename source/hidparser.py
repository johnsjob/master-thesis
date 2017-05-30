import hid
import sys
import json
import operator
from collections import OrderedDict
import time
class HIDReport:

    def _fmt_hid_piece(self, piece):
        result = list(piece)
        #result.reverse()
        #return int(reduce(lambda x,y: str(x)+str(y), result), 16)
        if len(piece) == 1:
            import pdb; pdb.set_trace()
        return int(reduce(lambda x,y: str(y)+str(x), result), 16)

    def __init__(self, report=[], hex_report=None):
        result = OrderedDict()

        if hex_report is None:
            hex_rep = map(lambda val: format(val, 'x'), report)
        else:
            hex_rep = hex_report

        result['id']       = int(hex_rep[0], 16)
        result['buttons']  = int(hex_rep[1], 16)
        result['x']        = self._fmt_hid_piece( hex_rep[2:2+4] )
        result['y']        = self._fmt_hid_piece( hex_rep[6:6+4] )
        result['pressure'] = self._fmt_hid_piece( hex_rep[-10:-10+2] )
        result['tiltx']    = self._fmt_hid_piece( hex_rep[-8:-8+2] )
        if result['tiltx'] & 32768 == 32768:
            result['tiltx'] -= 65535
        result['tilty']    = self._fmt_hid_piece( hex_rep[-6:-6+2] )
        if result['tilty'] & 32768 == 32768:
            result['tilty'] -= 65535
        result['twist']    = self._fmt_hid_piece( hex_rep[-4:] )
        result['raw']      = list(hex_rep)
        self.report = OrderedDict(result)

    def __str__(self):
        return json.dumps(self.report, indent=4)

    @property
    def data(self):
        return dict(self.report)

    @property
    def pos(self):
        return self.report['x'], self.report['y']

    @property
    def angles(self):
        return [self.report[key]/100.0 for key in ['tiltx','tilty','twist']]

    @property
    def buttons(self):
        return self.report['buttons']

    @property
    def pressure(self):
        return self.report['pressure'] / 2047.0

class BTPen:
    def _enum(self):
        devices = hid.enumerate()
        if not devices:
            raise Exception('No bluetooth device found!')

        for dev in devices:
            if 'pen' in dev['product_string'].lower():
                return dev

    def vid(self):
        return self._enum()['vendor_id']

    def pid(self):
        return self._enum()['product_id']

    def __init__(self):
        self.dev = hid.device()
        try:
            self.dev.open(self.vid(), self.pid())
        except IOError:
            raise IOError('Could not connect to bluetooth pen!')
        print json.dumps(self._enum(), indent=4)
        self._init = True

    def read(self):
        hid_report = self.dev.read(20)
        print hid_report
        return HIDReport(hid_report)

    def record(self, num_reports=100, num_warmup=50):
        res = []
        while len(res) < (num_reports + num_warmup):
            read_hid = self.read()
            data = read_hid.data
            res.append(data)
            print 'Coord #{}'.format(len(res))
        return res[num_warmup:]

    def __del__(self):
        if self._init:
            self.dev.close()

            del self.dev
            self.dev = None

            del self._init
            self._init = None

from numpy import array as mat
import pylab as p

def list_to_spaced_vals(x):
    res = reduce(lambda x,y: str(x) + ' ' + str(y), x)
    res = 'E: 3.14 20 {}\n'.format(res)
    return res

def filter_func(p):
    if (p[0] < 1000) or (p[1] < 1000):
        return False
    else:
        return True

if __name__ == '__main__':
    pen = BTPen()
    reports = pen.record(num_reports=900, num_warmup=20)
    del pen

    with open('./test.txt','w+') as fp:
        lines = [list_to_spaced_vals(r['raw']) for r in reports]
        fp.writelines(lines)

    positions = mat([[p['x'], p['y'], p['tiltx']/100.0, p['tilty']/100.0, p['twist']/100.0] for p in reports])
    pos = mat(filter(filter_func, positions))
    pos_err = mat(filter(lambda x: not filter_func(x), positions))

    from pylab import *
    ad_to_mm = 0.3
    mm_to_cm = 0.1
    correction_factor = 0.0625
    cm = ad_to_mm * mm_to_cm * correction_factor

    plot( (pos[:,0] - pos[0,0])*cm,
         -(pos[:,1] - pos[0,1])*cm, 'b')
    xlabel('x-axis (centimeters)')
    ylabel('y-axis (centimeters)')
    grid()
    show()

    plot(positions[:,0], -positions[:,1], 'r')
    xlabel('x-axis (anoto distance)')
    ylabel('y-axis (anoto distance)')
    grid()
    show()

    plot( positions[:,0]*cm,
         -positions[:,1]*cm, 'r', linewidth=0.3)
    plot( pos[:,0]*cm,
         -pos[:,1]*cm, 'b')
    xlabel('x-axis (centimeters)')
    ylabel('y-axis (centimeters)')
    grid()
    show()
