from __future__ import division
from helperfunctions import *
import plane_relative as plane_tools
import numpy

from numpy import arccos
#----------------------------------------#
num_points = 26 #num_points = 12 absolute minimum, actually 12+1
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
o = mat([1000*rand(), 1000*rand(), 1000*rand()])
R = rotation_matrix_rot_tilt_skew( rand()*360-180, rand()*90-45, rand()*360-180 )
dirx = R[:,0]
diry = R[:,1]
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

def calibrate_tool0():
    pass

for k in xrange(0,N):
    rot = rand()*360 - 180
    tilt = rand()*90 - 80
    skew = rand()*360 - 180
    px = 300*rand()
    py = 300*rand()
    l_xprim.append([px, py])

    Rrel = plane_tools.get_plane_relative_R(plane, rot, tilt-180, skew)
    l_R.append(Rrel)
    
    d = matmul_series(Rrel, do)
    Xtcp0 = get_plane_point(plane, px, py) - d
    l_xtcp.append(Xtcp0)

l_xtcp = mat(l_xtcp)
l_R = mat(l_R)
l_xprim = mat(l_xprim)

#diffs
dxtcp = diff(l_xtcp, axis=0)
dR = diff(l_R, axis=0)
dxprim = diff(l_xprim, axis=0)

lhs = []
rhs = []
l_cond = []
l_err = []
from numpy.linalg import solve, det, inv, cond
for i in xrange(0,N-1):
    A = sys(dxprim[i,0], dxprim[i,1], dR[i])
    b = dxtcp[i]
    lhs.append(A)
    rhs.append(b)
    if i > 0:
        a = mat(lhs[0:-1]); a = a.T.dot(a)
        l_cond.append(cond(a))
        #length of condition numbers is N-3
        L = mat(lhs).T.dot(lhs)
        R = mat(lhs).T.dot(rhs)

        r = solve(L, R)
        comp = mat([dirx, diry,do, do, do])
        res = mat([r[0,:], r[1,:], r[2:5,0], r[5:8,1], r[8:11,2]])
        err = norm(comp-res)
        l_err.append(err)
    
lhs = mat(lhs)
rhs = mat(rhs)

L = lhs.T.dot(lhs)
R = lhs.T.dot(rhs)

r = solve(L,R)
comp = mat([dirx, diry,do, do, do])
res = mat([r[0,:], r[1,:], r[2:5,0], r[5:8,1], r[8:11,2]])
err = abs(comp-res)
print r[2:5,0] - do
print vec_diff(r[2:5,0],do)
print
print
print r[5:8,1] - do
print vec_diff(r[5:8,1],do)
print
print
print r[8:11,2] - do
print vec_diff(r[8:11,2],do)
print
print
t = range(3,N+1)
logcond = numpy.log10(l_cond)
logerr = numpy.log10(l_err)
plot(t, logcond, label='Condition number',linewidth=2);
plot(t, logerr, label='Error (frobenious norm)',linewidth=2);
hlines(-1, t[0], t[-1], label='Tolerance = 10^-1')
xlim(t[0], t[-1])
xlabel('Number of points collected', fontsize=14)
ylabel('log10', fontsize=14)

index = 12-3
plt.annotate("number of points = 12",
            xy=(t[index]+0.1, logerr[index]+0.2), xycoords='data',
            xytext=(t[index]+0.5, logerr[index]+5), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
grid()
title('Calibration algorithm verification')
legend()
show()
