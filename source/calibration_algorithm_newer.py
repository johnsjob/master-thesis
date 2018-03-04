# -*- coding: cp1252 -*-
from __future__ import division
#----------------------------------------#

from plotsettings import PlotSettings

import numpy
import numpy as n
import utils

from pylab import xlim, ylim, xticks, yticks, savefig, clf
import pylab

from numpy.linalg import solve, lstsq, det, inv, cond, svd
from numpy import array as mat, log10, diag
#----------------------------------------
# custom imports
from helperfunctions_math import *
from helperfunctions_plot import *
import plane_relative as plane_tools

from orientation_verification import fit_to_ori, solve_ori

from numpy.random import set_state
import os.path as path

rand_state = ('MT19937', n.array([4272543851, 2552595568,  431289734,  494160517,   15621524,
       3365871479, 3691534276,  705774780, 1590143843, 3193439880,
       1039175963, 1244054132, 1422283452, 2857425769, 1167033376,
       1816213813,   70429566,  360206226, 2171151328, 1097054948,
        505766291,  770900782, 3287631580, 1576770853,  830622809,
       2749681270, 2891043295, 2465256331, 2254136045, 2070768961,
       3687692167, 2655748804, 3242291735, 2270477560, 1984549220,
       1039808922, 3546118402, 1780548169, 2206121502, 1023426340,
       3295624999, 3484651512, 2014309578,  157984313, 4251147539,
       4005212747, 4140390252, 2815973799, 2352408759, 2635162117,
       1998696550, 3365582203, 3251800172,  750425027, 3935112254,
         94389170, 1067726855, 1580075728, 1103897548,  837728437,
       3525864943, 3199603848, 3485782908, 3926313547,  225002741,
        296702674, 3330584280, 2004615401, 2190257401, 4034219265,
        952995947, 2910526444, 1175371011,  775006296, 2264902488,
       1746018010, 2346434199, 1329250084, 1735925444, 1565739309,
       3163213225, 2810472833,  278486153, 2003253365,  767782526,
       2576694361, 3190253534, 4007734294,  963137327, 2216314375,
       1171287175, 3191849216,   90662231, 1113854804, 2018893147,
       2160543852, 1966471807, 2665628543, 1478092281, 1356191861,
       2377413350, 1685215606, 3200434710,  215713147, 1149009508,
        750234647,   51285890, 3415901350, 1908805302, 2564823793,
       4086943532, 1154465625, 2812582711, 1313698744, 1464542074,
        628575564, 3467442113,  262538024, 1539349489,  105778280,
       1340439291, 3193774558, 3820626456, 3970067007,  460380309,
       2007024001, 1203128358, 2396521974, 3213948721, 4204815225,
       2991889344, 3428319561,  980196822, 1192137901, 2095491766,
         24163295, 1768641900, 2466640856,  478847655, 3516082741,
       2583191275, 1499773476,  520213283, 2904960824, 3121130600,
        138830350,  999288673, 1469756142, 2819112499,  372182893,
       3358203788, 1919810198, 1700817935, 1268711542, 2735058923,
       3698166596, 3592470875, 1368613993, 4259997437,  895950648,
        375199795,  774762978, 2763549185, 3563597156,  572695698,
       2611945041, 2297354438, 3989387496, 4012355517, 1069018611,
        569258509,  422178708, 3368811776,  871572870,  781555363,
       3800004077,  655658139, 1740519675, 3977836839, 3162217913,
       2325641008, 3827099172, 2472667687, 3800504732,  359395637,
       4035050265, 2702378682, 3524593550, 1639473156, 2102392035,
       2870463327, 3462917578,  531643722, 3249546712, 2908486439,
       3447916333, 3226549081, 1417659903, 3037320711, 2616149503,
       2174820072, 3715952665, 1945003336, 2148002419,  814185234,
       2533198722, 3772534311, 2998061779, 2646442629, 2097775577,
       3380846124,  505612553, 2042176558,  618795428, 4154865382,
       2189204612, 1764893667, 2247003868, 3933985315, 2462735047,
       1362223686, 2887465562,  238939174, 3943027662, 2038460028,
            80381,  286698100, 3247832378,  826183289, 3647839909,
       3178952381, 3996773947, 2963317493, 1500400204, 4171599576,
        730836410, 1911450816, 3240989624, 2330932312,  290607160,
       1088140051,  977262652,  278473212, 2033692573, 3684321306,
       4191372969,  401225524, 3169088221,  789882877, 2237934883,
        723884282, 2459285154, 1712087956, 1313904242, 2232644616,
       1430538682,  902379449,  254374405, 1265773558, 1735109162,
        320129303, 2409348358, 2034038923, 2164002746, 2641350580,
       1588030904, 3842648347, 3959928884, 1265010665,  664970386,
        350297868,   81402069, 2367789874, 2133391915, 3070152174,
       2375128869, 1489053106, 2688467162, 2539861783, 3750888536,
       3512641095, 1007274389,  833590728, 2557639837, 3891008268,
        525677895, 2200581025, 2252148149, 3443165090, 2266534990,
       2868283426, 2868512463, 4288041962, 1075859169, 3548892157,
        561900143, 1526758428,  607374607, 3383668439, 2970490756,
        211972893, 1832196058, 2115095330, 2501754866, 1926057408,
       3549121423, 1589457972, 2187205360, 1719871194, 2334643159,
       2962401949, 1059952969, 3865105532,  672994583, 1887369648,
       2781265323,  996538778, 1189130349, 3056710030,  590287236,
        200389202,  325875651, 3254216101, 4251766720, 4186165781,
       1256180691, 2330800059, 2601484348, 2517662995, 2749179712,
       2975626613, 2178140129, 2186592243, 1008130172,   42074057,
        411585630, 3330274418,  531862994,  673007415, 4083379938,
       4146394531, 2286532701,  600240366, 1920359267,  787679106,
       3891476599, 2960407034, 2589408893, 3420743896, 3927480023,
       2372244998, 4064616511,  902876861, 2548467162,  422027529,
       4040117777, 2738260440, 1054135964,   50005755, 3955338863,
       1885532058, 2405087531,  940946427,  828823819, 1614765256,
       1980783942, 2524125451, 3922525276, 1711421885, 4120533276,
       4156607440, 2088031088,  150940318,  681414452, 1471174076,
       3461639241, 2954913587,  457049925, 3033703346, 3036134707,
       1808980674, 2164115068, 3663971676,  481190969, 1850684964,
         85139049, 3976883833, 1446669662, 3814644778,  945663411,
       1880847787, 3676028404, 1424478057, 2645476837, 4102400879,
       3655124415, 3564104975, 1803425705, 1892102482, 1986327607,
       4193602623, 1695185914, 2938324753, 3651738432, 1866705792,
        756226745, 2969863054,  369803500, 3141222629, 2047565623,
        102323486, 2550219910,  188644475, 1837541557, 1417433452,
       3307837063, 2953587459,  277649125, 1551150828,  258332801,
       2238630914, 3410505080, 4180353887, 4076337939,  389146891,
        172129886, 3802627675,  160511488, 3371335186, 1819822241,
       3274026632, 3468441638, 2163142493, 1157119576, 3925046268,
       1509720544, 3181077116,  991943105, 4080655136,  961783984,
       3363850309, 3167951146, 2307041196,  570395480, 3125369893,
       2796815924, 2067159292,  790142040,  544833357, 2966770282,
       2034912314,  866398400, 3853906509, 4072598412, 2891900314,
       1155405070, 3959339346, 3760873061, 2188216760, 2121804981,
       1040864842, 3189370037, 1310915580, 2188368098, 1759115324,
       3760527494, 3007512096, 2564986159, 3577093393, 2959465442,
       1504684329, 3570793130, 3854430525, 3497147094,  669597867,
       3192879565,  323992815,    8935255,  479604783, 1663801453,
        680511234, 3590484256, 3948981571, 1411204847, 2707696879,
       3601963786, 2115429800, 2084455003, 1310873991, 2001513117,
        962867883, 4282299421,  469743232, 1492545227, 3543745808,
       4187728760, 3300987542, 3933372998, 1232840460, 3234050320,
       3258662897, 3671779643, 1208267876, 3867653262, 2830053242,
       3027278826,  195677582, 1137090720,  782207323, 4274308232,
       1270316589, 3265089858,  483628850, 4137739191, 3426891494,
       3526514702,  765589747, 3693642945,  264615155, 2720837542,
       4161552190, 2692607383, 1649762610, 2953446300, 2955299046,
        577744663, 1979434295, 4058532746,  555838069, 2988490402,
       3457918508, 4091451558,  354040455, 3488088385, 2385291494,
       2254273928, 3608399029, 3189261434, 3800007612,  513934922,
        971810313, 4104923056, 2311487010,  313115168,  646530070,
       2551787773, 3031517883, 2006876419,  610697920, 1335017057,
       1232294310, 1437206872, 3943964307, 4167387920, 1257089831,
       3746835411, 3977817655, 2230232675, 1183451324,  547832960,
       3084663401, 4237722783, 3662128372, 2089534800, 4143617969,
        729843910, 2448390574, 3511176101, 2085825118, 1706421885,
       3839138579, 3953141230,  682126295, 4010472807, 3891311201,
       2335765342, 3715598444, 3678594326, 2596063472, 1767343251,
       1149299825, 2506779807, 3062062067, 3033092796, 3111884820,
       2277969505, 1813050701,  538924208, 2169063902, 4165001438,
       3005211200, 1864198589,  966125046, 4058317216, 3757960583,
       2811244655, 2656576553, 1485285916, 2209225635, 2685539193,
       3570828645, 3945713571, 1590746505, 2535445560, 2264455438,
        401303369, 2007599862, 2898370386,  307080526, 1124921745,
       3663695897,  893946809, 1639603531, 3105583879, 3520078117,
       2976634573, 2595753638,  837585180,  116712083, 1912073772,
       3936174162, 2493230111, 3725515935, 3170110637], dtype=n.uint32), 624, 0, 0.0)

set_state(rand_state)
#----------------------------------------#

# num_points = 12 absolute minimum, actually 12+1
num_points = 500
#========================================#
# placing paper origin
o = mat([1000*rand(), 1000*rand(), 1000*rand()])

### defining paper orientation
r,t,s = -rand_range(-180,180), rand_range(-180,180), rand_range(-180, 180)

# define the paper-orientation in global (robot-system) directions
# however, the paper uses "euclidean" orientation ('local' = euclidean)
plane = plane_tools.define_plane_from_angles(o, r, t, s, 'local')

#################################################
# define the units so less mistakes are made
chosen_unit = 'mm'

unit_descriptor = {
    'm' : 1.0,
    'mm': 1000.0
    }

unit = unit_descriptor[chosen_unit]
metre = unit
millimetre = unit / 1000.0
#----------------------------------------
# define delta vector which we want to find (in tool-space / local space)
# in unit lengths
L = 233 * millimetre
local_delta_vector = mat([1,2,3])
local_delta_vector = (local_delta_vector / norm(local_delta_vector))*L #length L

# Orientation of the tool in tool (local) coordinate system
local_tool_orientation = rotation_matrix_rot_tilt_skew(-10, 20, 30)
#----------------------------------------
# define the anoto point spread in unit lengths
#plane_point_spread = 47 * millimetre
plane_point_spread = 200 * millimetre

# define max-tilt of pen
pen_max_tilt = 40
#pen_max_tilt = 10

#----------------------------------------
def merge_dicts(*list_of_dicts):
    # init
    ret = {}
    keys = []

    # get all unique keys
    for d in list_of_dicts:
        keys += d.keys()
    keys = set().union(keys)

    # for all keys ...
    for k in keys:
        # prepare a k:th-list if none exists
        if not ret.has_key(k):
            ret[k] = []
        # for all dicts ...
        for d in list_of_dicts:
            # if dict has key ...
            if d.has_key(k):
                # check so that the key is not an empty list ...
                empty = False
                try:
                    empty = len(d[k]) == 0
                except:
                    # not a list/array-type, equivalent to non-empty list
                    pass
                # append item or non-empty list
                if not empty:
                    ret[k].append( d[k] )
    # for all keys ...
    for k in keys:
        # if we only got one item for this key from al the dicts ...
        if len(ret[k]) == 1:
            # un-list it
            ret[k] = ret[k][0]
        # remove empy lists if any manage to get here
        elif len(ret[k]) == 0:
            del ret[k]
        else:
            # turn remaining lists into numpy-arrays
            ret[k] = mat(ret[k])
    return ret
#----------------------------------------
def rad_to_ang(v):
    return v*180/pi
#----------------------------------------
def vec_ang(v,w):
    res = matmul(v,w) / (norm(v) * norm(w))
    return numpy.arccos(res) * 180 / numpy.pi
#----------------------------------------
def vec_diff(v1, v2):
    err = norm(v1 - v2)
    norm_err = abs(norm(v1) - norm(v2))
    angle_err = rad_to_ang(acos( (v1/norm(v1)).dot((v2/norm(v2))) ))
    return err, norm_err, angle_err
#----------------------------------------
def problem_formulation(dx, dy, dR):
    r11,r12,r13,r21,r22,r23,r31,r32,r33 = (-dR).reshape(9)
    S1 = [dx, dy, r11, 0,   0, r12, 0, 0,   r13]
    S2 = [0,   0, r21, dx, dy, r22, 0, 0,   r23]
    S3 = [0,   0, r31, 0,   0, r32, dx, dy, r33]
    row_value = 3
    col_value = 9
    return mat([S1, S2, S3]),  row_value, col_value
#----------------------------------------
def generate_random_Anoto_Point(L):
    px = L*rand()
    py = L*rand()
    return px, py
#----------------------------------------
def solve_tool0_tip_alt(array_forward_kinematics_T44,
                        array_anoto2D,
                        array_lhs_sys_eq = None):
    try:
        num_points,m,n = array_forward_kinematics_T44.shape
    except Exception as e:
        print 'solve_tool0_tip:\n\tWrong shape or type for input parameter: array_forward_kinematics_T44'
    try:
        m,n = array_anoto2D.shape
    except Exception as e:
        print 'solve_tool0_tip:\n\tWrong shape or type for input parameter: array_anoto2D'

    l_xtcp = array_forward_kinematics_T44[:, 0:3, 3]
    l_R = array_forward_kinematics_T44[:, 0:3, 0:3]

    dxtcp = diff(l_xtcp, axis=0)
    dR = diff(l_R, axis=0)
    danoto2D = diff(array_anoto2D, axis=0)

    lhs = []
    rhs = []
    l_cond = []
    l_err = []

    for i in xrange(0, num_points-1): #one less after forward-differences....
        A, row_value, col_value = array_lhs_sys_eq(danoto2D[i,0], danoto2D[i,1], dR[i])
        b = dxtcp[i]
        lhs.append(A)
        rhs.append(b)
    lhs = mat(lhs).reshape(((num_points-1) * row_value, col_value))

    #shape the rhs depending on shape-info from lhs
    if row_value != 1:
        rhs = mat(rhs).reshape((num_points-1) *  row_value)
    else:
        rhs = mat(rhs)

    L = lhs.T.dot(lhs)
    R = lhs.T.dot(rhs)

    X, Y, D = solve(L, R).reshape(3,3).T
    #X, Y, D = lstsq(lhs,rhs)[0].reshape(3,3).T
    X = X / numpy.linalg.norm(X)
    Y = Y / numpy.linalg.norm(Y)
    Z = cross(X, Y)
#    Z = Z / numpy.linalg.norm(Z)
    result = mat([X,Y,Z,D]).T
    condition = cond(L)

    return result, condition
#----------------------------------------
def generate_Xflange_orientation(plane,rot, tilt, skew):
    """
        Generate Xtcp-orientation in world coordinates, using Anoto-paper
        orientation formulation
        Planes are homoenous matrices, if we want the orientation
        we need the 0:3,0:3 submatrix.
    """
    return plane_tools.define_plane_relative_from_angles(plane, (0,0,0),
                                                         rot, tilt, skew,'global')[:3,:3]
#----------------------------------------
def setup_geometry(current_plane, point_spread, num_points, perturbations=None):
    global local_delta_vector, local_tool_orientation,\
           millimetre, pen_max_tilt

    geometry_info = {'plane':current_plane}
    geometry_info['local_tool_orientation'] = local_tool_orientation
    geometry_info['local_delta_vector'] = local_delta_vector

    geometry_info['correct_solution_geometry'] = mat(list(geometry_info['plane'][:3,:3].T.flatten()) +
               list(geometry_info['local_delta_vector'])).reshape(4,3).T


    #generating points and "forward-kinematics"
    collected_data = []
    for k in xrange(0,num_points):
        info = {}
        info['angles'] = \
        {
            'rot':  rand_range(-180,180),
            'tilt': rand_range(0, pen_max_tilt),
            'skew': rand_range(-90, 180)
        }
        # Xtcp (flange) orientation in global space, generated relative to the paper plane
        info['Xflange_orientation_relative_to_paper_plane'] = \
                                                generate_Xflange_orientation(geometry_info['plane'],**info['angles'])
        #generate pen-tip position in Anoto2d in mm
        px,py = generate_random_Anoto_Point(point_spread)
        info['pentip_2d'] = [px,py]

        # generate global Xtcp position in mm
        info['Xtcp0'] = (plane_tools.get_plane_point(geometry_info['plane'], px, py)[:3] - \
                 matmul( info['Xflange_orientation_relative_to_paper_plane'], local_delta_vector[:3]) )
        # ^OK

        # generate relative-tool-orientation in world coordinates
        info['global_tool_orientation'] = matmul( info['Xflange_orientation_relative_to_paper_plane'],
                                                  local_tool_orientation )
        # ^OK

        info['forward_kinematics'] = homogenous_matrix( info['Xflange_orientation_relative_to_paper_plane'],
                                                        info['Xtcp0'] )
        # perturbations
        if type(perturbations) in [list, tuple]:
            if 'tip' in perturbations:
                tilt = info['angles']['tilt']
                r = lambda: (0.6*rand()-0.3) * millimetre
                if tilt >= 0:
                    r2 = lambda: (0.02*rand()-0.01)*(1-abs(info['angles']['tilt'])/numpy.max(abs(info['angles']['tilt']))) * millimetre
                    info['pentip_2d'] = [px + r() + r2(), py + r() + r2()]
                else:
                    r2 = lambda: (0.02*rand()-0.01) * millimetre
                    info['pentip_2d'] = [px + r() + r2(), py + r() + r2()]
        # ^OK
        collected_data.append(info)
    geometry_info['data'] = merge_dicts(*collected_data)
    geometry_info['data']['angles'] = merge_dicts(*geometry_info['data']['angles'])

    return geometry_info
#----------------------------------------
def find_solution_pen_tip(geometry_info, included_solutions_from_start = -1):
    result, cond_num = solve_tool0_tip_alt(geometry_info['data']['forward_kinematics'][:included_solutions_from_start],
                                      geometry_info['data']['pentip_2d'][:included_solutions_from_start],
                                      problem_formulation)
    wobj2, _ = fit_to_ori(result[:,:3])
    result[:,:3] = wobj2
    return result, cond_num

def find_solution_pen_ori(geometry_info, included_solutions_from_start = -1):
    # solve for orientation s which should be same as local_tool_orientation
    #l,m,n = geometry_info['data']['Xflange_orientation_relative_to_paper_plane'].shape

    flange_orientation = geometry_info['data']['forward_kinematics'][:,:3,:3]
    pen_orientation = geometry_info['data']['global_tool_orientation'][:,:3,:3]

    solved_tool_orientation, _ = solve_ori(flange_orientation, pen_orientation)

    return solved_tool_orientation, 0.0
#----------------------------------------
def perform_solution_run(geometry_info):
    interval = range(3,num_points)
    list_of_solving = []
    for k in interval:
        solve_info = {}
        tip_wobj_res = find_solution_pen_tip(geometry_info, k)

##        solve_info['point_spread_x'] = numpy.std(geometry_info['data']['pentip_2d'][:k], axis=0)[0]
##        solve_info['point_spread_y'] = numpy.std(geometry_info['data']['pentip_2d'][:k], axis=0)[1]

        solve_info['tipwobj-result']   = tip_wobj_res[0]
        solve_info['tip-result']       = tip_wobj_res[0][:,3]
        solve_info['wobj-result']      = tip_wobj_res[0][:,:3]
        solve_info['tip-cond_num']     = tip_wobj_res[1]

        solve_info['orientation-result'], solve_info['orientation-cond_num'] = find_solution_pen_ori(geometry_info, k)
        sol2 = solve_info['orientation-result']
        solve_info['orientation-result'] = sol2

        solve_info['err-tipwobj'] = abs(geometry_info['correct_solution_geometry'] - solve_info['tipwobj-result'])
        solve_info['err-tip']     = numpy.linalg.norm(solve_info['err-tipwobj'][:,3])
        solve_info['err-wobj']    = numpy.linalg.norm(solve_info['err-tipwobj'][:,:3])

        solve_info['err-ori']     = numpy.linalg.norm(geometry_info['local_tool_orientation'] - sol2)
        list_of_solving.append(solve_info)

    solving_data = merge_dicts(*list_of_solving)
    solving_data['interval'] = interval
    print 'solution max  error tip(1) = {}\n'.format( numpy.max( abs( solving_data['err-tip'][1:]) ))

    print 'solution max  error tip(20) = {}'.format( numpy.max( abs( solving_data['err-tip'][20:]) ))
    print 'solution mean error tip(20) = {}\n'.format( numpy.mean( abs( solving_data['err-tip'][1:21]) ))

    print 'solution max  error tip(40) = {}'.format( numpy.max( abs( solving_data['err-tip'][40:]) ))
    print 'solution mean error tip(40) = {}\n'.format( numpy.mean( abs( solving_data['err-tip'][1:41]) ))
    print 'solution error ori = {}'.format( numpy.max( abs( solving_data['err-ori'][1:]) ))
    return solving_data

def make_plots(solving_data):
    global chosen_unit

    logcond = log10( solving_data['tip-cond_num'] )
    plot(solving_data['interval'], logcond,
                       'b--',label='Condition number tip/wobj',
                       linewidth=2)

    logerr  = log10( solving_data['err-tip'] )
    plot(solving_data['interval'], logerr,
                      'b',label='Error tip (frobenious norm)',
                      linewidth=2);

    logerr  = log10( solving_data['err-wobj'] )
    plot(solving_data['interval'], logerr,
                      'g',label='Error wobj (frobenious norm)',
                      linewidth=2);

    if chosen_unit == 'mm':
        tol = -1
        hlines(tol, solving_data['interval'][0],
                   solving_data['interval'][-1],
                   label='Tolerance = 10^{}'.format(tol))
    else:
        tol = -4
    hlines(-4, solving_data['interval'][0],
               solving_data['interval'][-1])

    xlim(solving_data['interval'][0], solving_data['interval'][-1])
    xlabel('Number of points collected', fontsize=PlotSettings.label_size)
    ylabel('log10'.format(chosen_unit), fontsize=PlotSettings.label_size)
    xticks(fontsize=PlotSettings.tick_size)
    yticks(fontsize=PlotSettings.tick_size)

    index = 4-3
    plt.annotate("number of points = 4",
                xy=(solving_data['interval'][index]+0.01, logerr[index]+0.2), xycoords='data',
                xytext=(solving_data['interval'][index]+0.8, logerr[index]+4.5), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3"),
                )
    grid()
    title('Calibration algorithm verification using simulated geometry', fontsize=PlotSettings.title_size)
    legend(fontsize=PlotSettings.legend_size)
    show()

#----- SCRIPT START ---------

if __name__ == '__main__':
  # generating data
    res = []
    ks = range(90) # default value is 90
    for k in ks:
        print 'Run {} % complete!'.format(100* k / len(ks))
        with utils.timing.Timer() as timer:
            try:
                print "Sampling points..."
                geometry_info = setup_geometry(plane, plane_point_spread,
                                               num_points, perturbations=['tip'])
                print 'Collecting solving information...'
                solving_data = perform_solution_run(geometry_info)

                res.append(solving_data)
                print
                print 'Preparing plots...'
            except Exception as e:
                print "Run #{} of {}: failed:".format(k+1, len(ks))
                print str(e)
#                raise
                continue

#--------------- generate plots ---------------------------------

##    make_plots(solving_data) # unused in report

    if chosen_unit == 'mm':
        length_tol = -1 # 1/10th millimetre
    else:
        length_tol = -4 # 1/10th millimetre (in meters)

    for key,tol, un in zip(['err-tip', 'err-tip-zoom', 'err-wobj'],
                           [length_tol, length_tol, -4],
                           [chosen_unit, chosen_unit, '']):
        val = mat([
                    x[key.replace('-zoom','')] for x in res
                    ])
        maxval = numpy.max(val, axis=0)
        minval = numpy.min(val, axis=0)
        meanval = numpy.mean(val, axis=0)
        if un:
            unit_str = '[{}]'.format(un)
        else:
            unit_str = ''
        clf() # clear figure
        key_label = key.replace('-', ' ').replace('zoom', '')
        plot(log10(maxval),'b', label = 'max {}'.format(key_label, unit_str))
        plot(log10(meanval),'g', label = 'mean {}'.format(key_label, unit_str))
        plot(log10(minval),'r', label = 'min {}'.format(key_label, unit_str))
        hlines(tol, solving_data['interval'][0],
                    solving_data['interval'][-1],
                    label='Tolerance = 10^{}'.format(tol, unit_str))
        legend(fontsize=PlotSettings.legend_size)
        if not key == 'err-tip-zoom':
          xlim(solving_data['interval'][0], solving_data['interval'][-1])
        else:
          xlim(solving_data['interval'][0], 101)
        xlabel('Number of measured points', fontsize=PlotSettings.label_size)
        ylabel('Log$_{{10}}$ error {}'.format(unit_str),
               fontsize=PlotSettings.label_size)
        xticks(fontsize=PlotSettings.tick_size)
        yticks(fontsize=PlotSettings.tick_size)
        if 'tip' in key:
            ylim(-4, 4)
        elif 'wobj' in key:
            ylim(-6, 1)
        title('Calibration algorithm verification with repetition '+\
              'using simulated geometry', fontsize=PlotSettings.title_size)
        grid()

        # Qt4 backend - maximize plots
        pylab.get_current_fig_manager().window.showMaximized()

        # save figures
        basepath = r"C:\Users\***REMOVED***\Dropbox\exjobb\results\calibrationalg_res"

        if key == 'err-tip':
          figpath = path.join(basepath, "tool_tip")
          filename = "tool_tip_calibration_tilt40_spread200_len233_2"
          filetype = ".png"
          savefig(path.join(
            figpath, filename + filetype), bbox_inches='tight')
        elif key == 'err-tip-zoom':
          figpath = path.join(basepath, "tool_tip")
          filename = "tool_tip_calibration_tilt40_spread200_len233_22"
          filetype = ".png"
          savefig(path.join(
            figpath, filename + filetype), bbox_inches='tight')
        elif key == 'err-wobj':
          figpath = path.join(basepath, "wobj")
          filename = "wobj_calibration_tilt40_spread200_len233_2"
          filetype = ".png"
          savefig(path.join(
            figpath, filename + filetype), bbox_inches='tight')
