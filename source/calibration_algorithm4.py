from __future__ import division
#----------------------------------------#

#custom imports
from helperfunctions_math import *
from helperfunctions_plot import *
import plane_relative as plane_tools
#----------------------------------------#
def convert_to_matrix(x):
    return mat(x)
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
#################################################
# define delta vector which
# we want to find (in tool-space / local space)
L = 100
do = mat([1,2,3,0])
do = (do / norm(do))*L #length L

# define the paper-orientation in global (robot-system) directions
plane = plane_tools.define_plane_from_directions(o, dirx, diry)

# Orientation of the tool in tool (local) coordinate system 
local_tool_orientation = rotation_matrix_rot_tilt_skew(-10, 20, 30)
#----------------------------------------
def merge_dicts(*list_of_dicts):
    ret = {}
    keys = []
    for d in list_of_dicts:
        keys += d.keys()
    keys = set().union(keys)
    for k in keys:
        if not ret.has_key(k):
            ret[k] = []
        for d in list_of_dicts:
            if d.has_key(k):
                empty = False
                try:
                    empty = len(d[k]) == 0
                except:
                    pass
                if not empty:
                    ret[k].append( d[k] )
    for k in keys:
        if len(ret[k]) == 1:
            ret[k] = ret[k][0]
        elif len(ret[k]) == 0:
            del ret[k]
        else:
            ret[k] = mat(ret[k])
    return ret
#----------------------------------------
def rad_to_ang(v):
    return v*180/pi
#----------------------------------------
def sys2(dx, dy, dR):
    r11,r12,r13,r21,r22,r23,r31,r32,r33 = (-dR).reshape(9)
    S1 = [dx, dy, r11, 0,   0, r12, 0, 0,   r13]
    S2 = [0,   0, r21, dx, dy, r22, 0, 0,   r23]
    S3 = [0,   0, r31, 0,   0, r32, dx, dy, r33]
    row_value = 3
    col_value = 9
    return mat([S1, S2, S3]),  row_value, col_value
#----------------------------------------
def vec_diff(v1, v2):
    err = norm(v1 - v2)
    norm_err = abs(norm(v1) - norm(v2))
    angle_err = rad_to_ang(acos( (v1/norm(v1)).dot((v2/norm(v2))) ))
    return err, norm_err, angle_err
#----------------------------------------
def generate_random_Anoto_Point(L):
    px = L*rand()-L/2.0
    py = L*rand()-L/2.0
    return px, py
#----------------------------------------
from numpy.linalg import solve, det, inv, cond
#----------------------------------------
def solve_tool0_tip_alt(array_forward_kinematics_T44, array_anoto2D, array_lhs_sys_eq = None):
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
        A, row_value, col_value = array_lhs_sys_eq(danoto2D[i,0], danoto2D[i,1], dR[i])
        b = dxtcp[i]
        lhs.append(A)
        rhs.append(b)
    lhs = mat(lhs).reshape(((N-1) * row_value, col_value))

    #shape the rhs depending on shape-info from lhs
    if row_value != 1:
        rhs = mat(rhs).reshape((N-1) *  row_value)
    else:
        rhs = mat(rhs)
        
    L = lhs.T.dot(lhs)
    R = lhs.T.dot(rhs)

    c = cond(L)
    r = solve(L,R)
    return r,c
#----------------------------------------
def generate_Xflange_orientation(plane,rot, tilt, skew):
    # Generate Xtcp-orientation in world coordinates (Rrel)
    # Planes are homoenous matrices, if we want the orientation
    # we need the 0:3,0:3 submatrix
    return plane_tools.define_plane_relative_from_angles(plane, (0,0,0), rot, tilt, skew)[:3,:3]
#----------------------------------------
if __name__ == '__main__':
    print "Sampling points..."
    point_spread = 300
    
    collected_info = {}
    collected_info['plane'] = plane
    collected_data = []
    #generating points and "forward-kinematics"
    l = []
    for k in xrange(0,N):
        info = {}
        info['angles'] = \
        {
            'rot':  rand_range(-180,180),
            'tilt': rand_range(-60, 60),
            'skew': rand_range(-180,180) 
        }
        # Xtcp (flange) orientation in global space, generated relative to the paper plane
        info['Xflange_orientation_relative_to_paper_plane'] = \
                                                generate_Xflange_orientation(collected_info['plane'],**info['angles'])
        #generate pen-tip position in Anoto2d in mm
        px,py = generate_random_Anoto_Point(point_spread)
        info['pentip_2d'] = [px,py]

        # generate global Xtcp position in mm
        info['Xtcp0'] = (plane_tools.get_plane_point(plane, px, py)[:3] - \
                 info['Xflange_orientation_relative_to_paper_plane'].dot(do[:3]))

        #generate relative-tool-orientation in world coordinates
        info['global_tool_orientation'] = matmul(info['Xflange_orientation_relative_to_paper_plane'], local_tool_orientation)
        info['forward_kinematics'] = homogenous_matrix(info['Xflange_orientation_relative_to_paper_plane'],
                                                       info['Xtcp0'])
        collected_data.append(info)
    collected_info['data'] = merge_dicts(*collected_data)

    print "Solving for dirx, diry, do..."
    sys_of_eq = sys2
    import time
    start_time = time.clock()

    #r, cond_num = solve_tool0_tip_alt(T44, l_anoto2D, sys_of_eq)
    r, cond_num = solve_tool0_tip_alt(collected_info['data']['forward_kinematics'],
                                      collected_info['data']['pentip_2d'],
                                      sys_of_eq)
    stop_time = time.clock()
    time_spent = stop_time - start_time

    #solve for orientation s which should be same as local_tool_orientation
    #s = linalg.solve(l_R.reshape(360,3).T.dot(l_R.reshape(360,3)), l_R.reshape(360,3).T.dot(l_Rd_rel.reshape(360,3)))
    s = linalg.solve(collected_info['data']['Xflange_orientation_relative_to_paper_plane'].reshape(360,3).T.dot(collected_info['data']['Xflange_orientation_relative_to_paper_plane'].reshape(360,3)),
                     collected_info['data']['Xflange_orientation_relative_to_paper_plane'].reshape(360,3).T.dot(collected_info['data']['global_tool_orientation'].reshape(360,3)))

    
    print
    print 'Time spent solving '+str(N)+' points: ' + str(time_spent) +' seconds.'


    print
    print 'Preparing plots...'
    comp = mat([dirx, diry, do[:3]]).T.reshape(9)
    l_cond = []
    l_err = []
    for k in xrange(3,N):
        #r,cond_num = solve_tool0_tip_alt(T44[0:k,:,:], l_anoto2D[0:k], sys_of_eq)
        r, cond_num = solve_tool0_tip_alt(collected_info['data']['forward_kinematics'][0:k,:,:],
                                      collected_info['data']['pentip_2d'][0:k],
                                      sys_of_eq)
        l_cond.append(cond_num)

        res = mat(r)
        err = abs(comp-res)
        l_err.append(norm(err))    
    print 'solution error = ' + str(res.reshape((3,3))[:,2] - do[:3])
    print 'err, norm_err, angle_err = ' + str(vec_diff(res.reshape((3,3))[:,2],do[:3]))

    t = range(3,N)
    logcond = numpy.log10(l_cond)
    logerr = numpy.log10(l_err)
    plot(t, logcond, label='Condition number',linewidth=2);
    plot(t, logerr, label='Error (frobenious norm)',linewidth=2);
    hlines(-1, t[0], t[-1], label='Tolerance = 10^-1')
    xlim(t[0], t[-1])
    xlabel('Number of points collected', fontsize=14)
    ylabel('log10', fontsize=14)

    index = 4-3
    plt.annotate("number of points = 4",
                xy=(t[index]+0.01, logerr[index]+0.2), xycoords='data',
                xytext=(t[index]+0.8, logerr[index]+4.5), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"),
                )
    grid()
    title('Calibration algorithm verification using simulated geometry')
    legend()
    show()
