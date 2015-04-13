from __future__ import division
from helperfunctions import *
import plane_relative as plane_tools
import numpy

from numpy import arccos
#----------------------------------------#
#num_points = 12 absolute minimum, actually 12+1
num_points = 120
#========================================#
#print; print "Init plots..."
#ax,_ = init_plot()
###========================================#
import numpy
mat = numpy.array

#num points we want to sample
N = num_points
#placing paper origin
o = mat([1000*rand(), 1000*rand(), 1000*rand()])
#defining paper orientation
R = rotation_matrix_rot_tilt_skew( -(rand()*360-180), rand()*90-60, rand()*360-180 )
#extracting basis vectors
dirx = R[:,0]
diry = R[:,1]
#delta vector which we want to find in tool-space
L = 100
do = mat([1,2,3])
do = (do / norm(do))*L #length L
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
from numpy.linalg import solve, det, inv, cond
def solve_tool0_tip(array_forward_kinematics_T44, array_anoto2D):
    try:
        N,m,n = array_forward_kinematics_T44.shape
    except Exception as e:
        print 'solve_tool0_tip:\n\tWrong shape or type for input parameter: array_forward_kinematics_T44'
    try:
        m,n = array_anoto2D.shape
    except Exception as e:
        print 'solve_tool0_tip:\n\tWrong shape or type for input parameter: array_anoto2D'

    l_xtcp = array_forward_kinematics_T44[:,0:3,3]
    l_R = array_forward_kinematics_T44[:,0:3,0:3]

    dxtcp = diff(l_xtcp, axis=0)
    dR = diff(l_R, axis=0)
    danoto2D = diff(array_anoto2D, axis=0)

    lhs = []
    rhs = []
    l_cond = []
    l_err = []

    for i in xrange(0,N-1): #one less after forward-differences....
        A = sys(danoto2D[i,0], danoto2D[i,1], dR[i])
        b = dxtcp[i]
        lhs.append(A)
        rhs.append(b)

    lhs = mat(lhs)
    rhs = mat(rhs)

    L = lhs.T.dot(lhs)
    R = lhs.T.dot(rhs)

    c = cond(L)
    r = solve(L,R)
    return r,c
#----------------------------------------
def extract_solutions(sol):
    dirx = sol[0,:]
    diry = sol[1,:]
    do1 = sol[2:5,0]
    do2 = sol[5:8,1]
    do3 = sol[8:11,2]
    return mat([dirx, diry, do1, do2, do3])
#----------------------------------------
if __name__ == '__main__':
    print "Sampling points..."
    #generating points and "forward-kinematics"
    for k in xrange(0,N):
        rot = 30 + rand()*0.001 - 0.0005
        tilt = 40 + rand()*0.001 - 0.0005
        skew = 50 + rand()*0.001 - 0.0005

        px,py = generate_random_Anoto_Point(0.001)
        append_to_points2D([px, py])

        Rrel = plane_tools.get_plane_relative_R(plane, rot, tilt, skew)
        append_to_relative_plane_orientation(Rrel)
        
        #this is technically correct but sloppy
        #compared to the derivation
        #d = matmul_series(Rrel, do)
        #Xtcp0 = get_plane_point(plane, px, py) - d

        Xtcp0 = plane_tools.get_plane_relative_skew_point(plane, px, py, rot, tilt, skew, do)
        append_to_Xtcp_o(Xtcp0)

    #convert the lists to ndarrays
    l_xtcp = convert_to_matrix(l_xtcp)
    l_R = convert_to_matrix(l_R)
    l_anoto2D = convert_to_matrix(l_anoto2D)

    #turn geometry into forward-kinematics
    T44 = zeros((num_points,4,4))
    T44[:,0:3,0:3] = l_R
    T44[:,0:3,3] = l_xtcp
    T44[:,3,:] = [0,0,0,1]

    print "Solving ..."
    import time
    start_time = time.clock()

    r, cond_num = solve_tool0_tip(T44, l_anoto2D)

    stop_time = time.clock()
    time_spent = stop_time - start_time
    print
    print 'Time spent solving '+str(N)+' points: ' + str(time_spent) +' seconds.'

    print
    print 'Preparing plots...'
    comp = mat([dirx, diry,do, do, do])
    l_cond = []
    l_err = []
    for k in xrange(3,N):
        r,cond_num = solve_tool0_tip(T44[0:k,:,:], l_anoto2D[0:k])
        l_cond.append(cond_num)

        res = mat([r[0,:], r[1,:], r[2:5,0], r[5:8,1], r[8:11,2]])
        err = abs(comp-res)
        l_err.append(norm(err))
    print r[2:5,0] - do
    print 'err, norm_err, angle_err = ' + str(vec_diff(r[2:5,0],do))
    print
    print
    print r[5:8,1] - do
    print 'err, norm_err, angle_err = ' + str(vec_diff(r[5:8,1],do))
    print
    print
    print r[8:11,2] - do
    print 'err, norm_err, angle_err = ' + str(vec_diff(r[8:11,2],do))


    t = range(3,N)
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
    title('Calibration algorithm verification using simulated geometry')
    legend()
    show()
