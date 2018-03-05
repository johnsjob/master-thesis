import hid
import sys
import json
import operator
from collections import OrderedDict
import time
import os.path as path

import hidparser
from numpy import array as mat

from pylab import plot, show, xticks, yticks, xlabel, ylabel, grid, savefig, clf

from plotsettings import PlotSettings

basepath = r"C:\Users\***REMOVED***\Dropbox\exjobb\results\hid_capture"

def main():
    lines = None
    with open(path.join(basepath, "raw_hid_data.txt")) as fp:
        lines = fp.readlines()
    reports = [hidparser.HIDReport(hex_report = x.split(" ")[3:-1]).data for x in lines]

    positions = mat([[p['x'], p['y'], p['tiltx']/100.0, p['tilty']/100.0, p['twist']/100.0] for p in reports])
    pos = mat(filter(hidparser.filter_func, positions))
    pos_err = mat(filter(lambda x: not hidparser.filter_func(x), positions))

    ad_to_mm = 0.3
    mm_to_cm = 0.1
    correction_factor = 0.0625
    cm = ad_to_mm * mm_to_cm * correction_factor

    plot( (pos[:,0] - pos[0,0])*cm,
         -(pos[:,1] - pos[0,1])*cm, 'b')
    xlabel('x [cm]', fontsize=PlotSettings.label_size)
    ylabel('y [cm]', fontsize=PlotSettings.label_size)
    xticks(fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    grid()
    savefig(path.join(basepath, "johnnys_local.png"),
        bbox_inches='tight')
    clf()

    plot(positions[:,0], -positions[:,1], 'r')
    xlabel('x [ad]', fontsize=PlotSettings.label_size)
    ylabel('y [ad]', fontsize=PlotSettings.label_size)
    xticks(fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    grid()
    savefig(path.join(basepath, "johnnys_errors_global.png"),
        bbox_inches='tight')
    clf()

    plot( positions[:,0]*cm,
         -positions[:,1]*cm, 'r', linewidth=0.3)
    plot( pos[:,0]*cm,
         -pos[:,1]*cm, 'b')
    xlabel('x [cm]', fontsize=PlotSettings.label_size)
    ylabel('y [cm]', fontsize=PlotSettings.label_size)
    xticks(fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)
    grid()
    savefig(path.join(basepath, "johnnys_all.png"),
        bbox_inches='tight')
    clf()


if __name__ == '__main__':
    main()
