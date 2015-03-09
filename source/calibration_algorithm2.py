from __future__ import division
from helperfunctions import *
import numpy
#----------------------------------------#
num_points = 100
###========================================#
#Part two of calibration algorithm
import denavit_hartenberg as RK
Xtcp = mat([0,0,0,1])

L = 0.1 #10cm
delta_prim = L*get_normalized_vector(rand(), rand(), rand())
delta_prim = mat(list(delta_prim) + [1])

data = {'T':[], 'Xtcp':[], 'Xprim':[]}
for _ in xrange(0,3):
    a = rand_range(-45, 45)
    b = rand_range(-45, 45)
    c = rand_range(-45, 45)
    d = rand_range(-45, 45)
    e = rand_range(-45, 45)
    f = rand_range(-45, 45)
    T44 = RK.calc_tool_IRB120(a,b,c,d,e,f)
    data['T'].append(T44)
    data['Xtcp'].append(T44.dot([0,0,0,1]))
    data['Xprim'].append(data['Xtcp'][-1][0:3] + delta_prim[0:3])
    data['Xprim'][-1] = mat(list(data['Xprim'][-1]) + [1])

