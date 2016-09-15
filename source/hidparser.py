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
        old_pos = self.read().data
        skipped = []
        total = []
        while len(res) < (num_reports + num_warmup):
            read_hid = self.read()
            data = read_hid.data
            total.append(data)
##            if (abs(data['x'] - old_pos['x']) > 1000) or (abs(data['y'] - old_pos['y']) > 1000):
##                print 'skipping:'
##                print json.dumps(data, indent=4)
##                skipped.append(data)
##                continue
            old_pos = dict(data)
            res.append(data)
            print 'Coord #{}'.format(len(res))
#            print json.dumps(data, indent=4)
        return res[num_warmup:], total, skipped

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

if __name__ == '__main__':
    pen = BTPen()
    reports, total, skipped = pen.record(num_reports=600, num_warmup=20)
    del pen
    pos = mat([[d['x'], d['y']] for d in reports])
    tx = mat([d['tiltx'] for d in reports])
    ty = mat([d['tilty'] for d in reports])
    tw = mat([d['twist'] for d in reports])

    with open('./test.txt','w+') as fp:
        lines = [list_to_spaced_vals(r['raw']) for r in reports]
        fp.writelines(lines)
