from __future__ import division
from helperfunctions import *
import plane_relative as plane_tools
import numpy

from numpy import arccos
#----------------------------------------#
num_points = 120 #num_points = 12 absolute minimum, actually 12+1
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
R = rotation_matrix_rot_tilt_skew( -(rand()*360-180), rand()*90-60, rand()*360-180 )
dirx = R[:,0]
diry = R[:,1]
do = mat([1,2,3])
do = (do / norm(do))*L
plane = plane_tools.define_plane(o, dirx, diry)

l_xtcp = []
l_anoto2D = []
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
#----------------------------------------
def generate_random_Anoto_Point(L):
    px = L*rand()
    py = L*rand()
    return px, py
#----------------------------------------
def append_to_points2D(x):
    l_anoto2D.append(x)
def append_to_relative_plane_orientation(x):
    l_R.append(x)
def append_to_Xtcp_o(x):
    l_xtcp.append(x)
#----------------------------------------
def convert_to_matrix(x):
    return mat(x)
#----------------------------------------
for k in xrange(0,N):
    rot = rand()*360-180
    tilt = rand()*90-60
    skew = rand()*360-180
    px,py = generate_random_Anoto_Point(1)
    append_to_points2D([px, py])

    Rrel = plane_tools.get_plane_relative_R(plane, rot, tilt, skew)

    append_to_relative_plane_orientation(Rrel)
    
    #this is technically correct but sloppy
    #compared to the derivation
    d = matmul_series(Rrel, do)
    Xtcp0 = get_plane_point(plane, px, py) - d

    append_to_Xtcp_o(Xtcp0)

#convert the lists to ndarrays
l_xtcp = convert_to_matrix(l_xtcp)
l_R = convert_to_matrix(l_R)
l_anoto2D = convert_to_matrix(l_anoto2D)

#diffs
dxtcp = diff(l_xtcp, axis=0)
dR = diff(l_R, axis=0)
danoto2D = diff(l_anoto2D, axis=0)

import time
start_time = time.time()
lhs = []
rhs = []
l_cond = []
l_err = []
from numpy.linalg import solve, det, inv, cond
for i in xrange(0,N-1):
    A = sys(danoto2D[i,0], danoto2D[i,1], dR[i])
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
        comp = mat([dirx, diry, do, do, do])
        res = mat([r[0,:], r[1,:], r[2:5,0], r[5:8,1], r[8:11,2]])
        err = norm(comp-res)
        l_err.append(err)
    
lhs = mat(lhs)
rhs = mat(rhs)

L = lhs.T.dot(lhs)
R = lhs.T.dot(rhs)

r = solve(L,R)
stop_time = time.time()
time_spent = stop_time - start_time
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
print 'Time spent: ' + str(time_spent)
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
            xy=(t[index]+0.01, logerr[index]+0.2), xycoords='data',
            xytext=(t[index]+0.5, logerr[index]+5), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
grid()
title('Calibration algorithm verification')
legend()
show()
