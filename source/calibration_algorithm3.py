from __future__ import division
from helperfunctions import *
import plane_relative as plane_tools
import numpy

from numpy import arccos
#----------------------------------------#
num_points = 30
#========================================#
#print; print "Init plots..."
#ax,_ = init_plot()
###========================================#
##This part tests the theory of the algorithm
##in the case of a flat surface where it works even when varying orientation
##(I did not expect this, but the power of LS-methods are stronger than I thought?
print "\nSecond stage of algorithm - find delta-prim"
import numpy
mat = numpy.array

L = 100

N = num_points
o = mat([rand(), rand(), rand()])
dirx = [1,0,0]
diry = [0,1,0]
do = mat([1,2,3])
do = (do / norm(do))*L
plane = plane_tools.define_plane(o, dirx, diry)

l_xtcp = []
l_xprim = []
l_R = []

def rad_to_ang(v):
    return v*180/pi

def sys(dx, dy, dR):
    S = [dx, dy] + (-dR).reshape(9).tolist()
    return mat(S)

#doesn't work eventhough theoretically the same
def sys2(dx, dy, dR):
    S = [dx,dx,dx,dy,dy,dy] + (-dR).reshape(9).tolist()
    return mat(S)

def vec_diff(v1, v2):
    err = norm(v1 - v2)
    norm_err = abs(norm(v1) - norm(v2))
    angle_err = rad_to_ang(arccos( (v1/norm(v1)).dot((v2/norm(v2))) ))
    return err, norm_err, angle_err

for k in xrange(0,N):
    rot = rand()*360 - 180
    tilt = rand()*360 - 180
    skew = rand()*90 - 45
    px = 1000*rand()
    py = 1000*rand()
    l_xprim.append([px, py])

    Rrel = plane_tools.get_plane_relative_R(plane, rot, tilt, skew)
    l_R.append(Rrel)
    
    d = matmul_series(Rrel, do)
    Xtcp0 = get_plane_point(plane, px, py) - d
    l_xtcp.append(Xtcp0)

l_xtcp = mat(l_xtcp)
l_R = mat(l_R)
l_xprim = mat(l_xprim)

dxtcp = diff(l_xtcp,0)
dR = diff(l_R,0)
dxprim = diff(l_xprim,0)

lhs = []
rhs = []
l_cond = []
from numpy.linalg import solve, det, inv, cond
for i in xrange(0,N):
    A = sys(dxprim[i,0], dxprim[i,1], dR[i])
    b = dxtcp[i]
    lhs.append(A)
    rhs.append(b)
    if i > 0:
        a = mat(lhs[0:-1]); a = a.T.dot(a)
        l_cond.append(cond(a))
    
lhs = mat(lhs)
rhs = mat(rhs)

L = lhs.T.dot(lhs)
R = lhs.T.dot(rhs)

r = solve(L,R)
plot(numpy.log(l_cond)); show()
