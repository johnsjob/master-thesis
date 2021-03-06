import random

import numpy, numpy as n
from numpy import pi, linspace, meshgrid as mesh, zeros,\
     arccos as acos, log10, array, array as mat, uint32

from numpy.linalg import norm, det, inv
from numpy.random import uniform, seed, set_state, get_state
from helperfunctions_plot import *
import calibration_algorithm as cal
from helperfunctions_math import rotation_matrix_rot_tilt_skew as ori,\
                                 rotation_matrix_x_y_z as xyz,\
                                 homogenous_matrices, nzip, nmap,\
                                 quat_to_rot, rot_to_quat,\
                                 rot_tilt_skew as rts
from pylab import axhline
from plane_relative import generate_symmetric_curve,\
                           get_transformed_points, attach_to_base_frame,\
                           create_circle_curve, place_curve, attach_frames,\
                           orientation_frames

from denavit_hartenberg140 import forward_kinematics, DH_TABLE as dh_params,\
     calc_valid_invkin_irb140 as invkin

from denavit_hartenberg import homogenous_matrix as hom
from jacobian import jacobian_from_joints

from standardplot import StPlot
import pylab as plt

import itertools as it

import sys
import time
import json

import utils
import cPickle as pickle

from collections import OrderedDict

# settings
numpy.set_printoptions(precision=4)
numpy.set_printoptions(suppress=True)

rand_state = ('MT19937', array([4272543851, 2552595568,  431289734,  494160517,   15621524,
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
       3936174162, 2493230111, 3725515935, 3170110637], dtype=uint32), 624, 0, 0.0)

# functions
normalize = lambda x: x / norm(x)

def print_dict(d):
    print json.dumps(d, indent=4)

def targets(xmin,xmax,ymin,ymax, num):
    '''
    Generates a ('num'X'num') rectangle grid of targets in the range of
    [xmin, xmax], [ymin, ymax]. Depending on application these can either be
    relative to a work boject or absolute.

    return:
        ndarray : (num**2, 4, 4)
    '''
    x = linspace(xmin,xmax,num)
    y = linspace(ymin,ymax,num)
    xy = mat(mesh(x,y)).T
    res = nmap(lambda row: nmap(lambda pair: hom(0,0,0,[pair[0],pair[1],0]) ,row), xy)
    a,b,c,d = res.shape
    res = res.reshape((a*b,c,d))
    return res

def ik(pose):
    return invkin(pose)[0]

def flange(pose):
    '''
    Given a pose, calulate the forward-kinematics frame obtained
    from reaching the point by inverse-kinematics.
    '''
    return forward_kinematics(*ik(pose), **dh_params)['flange']

def pen_pos(pose, dim=2):
    return pose[:dim,3]

def pen_ori(pose):
    return pose[:3,:3]


def err(mag):
    mag = mag * 1e-3 # millimeters to meters
    return uniform(-mag, mag)

def err_pen_pos(mag, dim=2):
    mag = mag * 1e-3 # millimeters to meters
    ret = mat([uniform(-mag, mag) for _ in range(dim)])
    ret[2:] = 0
    return

def err_pen_ori(mag):
    return mat([uniform(-mag, mag) for _ in range(9)]).reshape(3,3)

def err_flange(mag):
    '''
        noise function used to generate a homogenous matrix with(-mag, mag)
        noise on the translation part.

        unit of mag is in milimeters
    '''
    mag = mag*1e-3 # millimeters to meters
    pos = mat([uniform(-mag, mag) for x in range(3)])
    res = zeros((4,4))
    res[:3,3] = pos
    return res

def solution_run(geometry_info):
    num_meas = len(geometry_info['data']['forward_kinematics'])
    res = [ cal.find_solution_pen_tip(geometry_info, k)
            for k in range(3, num_meas) ]
    return mat(res), [k for k in range(3, num_meas)]

def solution_errors(solution, comparison):
    c = solution[1]
    solution = solution[0]
    wobj = comparison[:,:3]
    tool = comparison[:,3]

    tool_err = abs(tool - solution[:,3]) * 1000
    wobj_err = abs(wobj - solution[:,:3])
    errors = {
        'wobj': {
            'norm': norm(wobj_err),
            'x_ang': acos(wobj[:,0].dot(normalize(solution[:,0])))*180.0/pi,
            'y_ang': acos(wobj[:,1].dot(normalize(solution[:,1])))*180.0/pi,
            'n_ang': acos(wobj[:,2].dot(normalize(solution[:,2])))*180.0/pi,
            'unit': 'deg'
            },
        'tool': {
            'norm': norm(tool_err),
            'x': tool_err[0],
            'y': tool_err[1],
            'z': tool_err[2],
            'unit': 'mm'
            },
        'cond': c
        }
    return errors

if __name__ == '__main__':
    # init - random generator
    set_state(rand_state)

    # init - misc
    num_points_side = 20
    # rot, tilt, skew
    # debug purposes
    rot  = zeros(n.round((num_points_side**2)**(1/3.0)))
    tilt = linspace(-20, 20, n.round((num_points_side**2)**(1/3.0)))
    skew = linspace(0, 45, n.round((num_points_side**2)**(1/3.0)))
    # best for now
    rot  = zeros(n.round((num_points_side**2)**(1/3.0)))
    tilt = linspace(-20, 20, n.round((num_points_side**2)**(1/3.0)))
    skew = linspace(-90, 90, n.round((num_points_side**2)**(1/3.0)))
    angles = mat(list(it.product(rot, tilt, skew)))

    # init - tool
    joints = (0,0,0,0,0,0)
    dh_params['tool'] = hom(0,0,0,
                           [0.005,0.005,0.233])
    dh_params['tool'][:3,:3] = ori(50,-45,0).dot(dh_params['tool'][:3,:3])
    dh_params['tool'][:3,:3] = dh_params['tool'][:3,:3].dot(ori(0,180,90)).dot(xyz(5,-10,0))
    tool = dh_params['tool']

    # init - wobj
    wobj = hom(180-45,180,0,[0.2,-0.2,-0.2])
    wobj[:3,:3] = wobj[:3,:3].dot(ori(0,180,90))

    # init reference
    ref = zeros((3,4))
    ref[:3,:3] = wobj[:3,:3]
    ref[:,3] = tool[:3,3]

    # init - paper targets
    tars = targets(0, 100e-3,
                   0, 200e-3, num=num_points_side)
    #tars = nmap(lambda x: x.dot(hom(ori(0,20,uniform(-180,180)))), tars)
    tars = nmap(lambda x: x[0].dot(hom(ori(x[1][0], x[1][1], x[1][2]))), zip(tars,angles))

    # init - global targets
    global_tars = nmap(wobj.dot, tars)

    # calculate flanges from global targets
    flanges = nmap(flange, global_tars)
    flanges = nmap(lambda x: x + err_flange(0.03), flanges)

    # obtain relative pen positions from local targets
    # unit: meters
    penxy0 = nmap(pen_pos, tars)

    # convert to 'anoto coords' with projection error
    # unit: meters
    penxy = nmap(lambda x: pen_pos(x) + mat((-(x[2,0]/0.70710678118654757)*0.8e-3*2,
                                             -(x[2,1]/0.70710678118654757)*0.8e-3*2)), tars)
    #penxy = nmap(lambda x: pen_pos(x)), tars)

    # convert to "anoto coords" in ad "anoto distance"
    penxyanoto = penxy*1000/0.3 + 1e9

    # obtain relative pen orietnation from local targets
    penori0 = global_tars
    penori = nmap(lambda x: hom(x[:3,:3] + err_pen_ori(2e-3)), tars)

    # supply information to calibration algorithm
    geometry_info = {
            'data':
                {
                    'forward_kinematics': flanges,
                    'pentip_2d': penxyanoto,
                    'global_tool_orientation': penori
                }
        }

    # save data to pickle
    pickle_data = []
    for k in xrange(len(flanges)):
        tmp = {
            'pen_pos': geometry_info['data']['pentip_2d'][k],
            'pen_ori': geometry_info['data']['global_tool_orientation'][k][:3,:3],
            'pen_fsr': numpy.NAN,
            'pen_fsr_adc': numpy.NAN,
            'robot_joints': ik(global_tars[k]),
            'robot_tcp': forward_kinematics(*global_tars[0], **dh_params)['tcp'],
            'robot_flange': forward_kinematics(*global_tars[0], **dh_params)['flange'],
            }
        pickle_data.append(tmp)

    with open('calibration_verifitcation.pickle', 'wb+') as fp:
        pickle.dump(pickle_data, fp)

    with open('tool_wobj.yaml', 'w+') as fp:
        json.dump({
            'tool': {
                'pos': list(dh_params['tool'][:3,3]),
                'q': list(rot_to_quat(dh_params['tool'][:3,:3]))
                },
            'wobj': {
                'pos': list(wobj[:3,3]),
                'q': list(rot_to_quat(dh_params['tool'][:3,:3]))
                }
            }, fp, indent=4)



    # perform calibration
    solution = cal.find_solution_pen_tip(geometry_info)
    d = solution_errors(solution, ref)
    print_dict(d)
    run, xvals = solution_run(geometry_info)
    errors = map(lambda x: solution_errors(x, ref), run)
    # perform calibration tool orientation
    (sol1,sol2,tensor,(u,s,v)),_ = cal.find_solution_pen_ori(geometry_info)
    print sol2[2] - tool[:3,:3]
    solved_tool_ori = sol2[2]


    # present results (tip, and wobj)
    plot(xvals,
         log10([e['tool']['norm'] for e in errors]))
    xlim(xvals[0],xvals[-1])
    axhline(-1,color='k')
    grid()
    show()

    plot(penxy0[:,0],penxy0[:,1],'k.')
    plot(penxy[:,0],penxy[:,1],'r.')
    grid()
    show()

#    plot(penxy0[:num_points_side,0],penxy0[:num_points_side,1],'k.')
#    plot(penxy[:num_points_side,0],penxy[:num_points_side,1],'r.')
#    grid()
#    show()

    # plot setup
    pl = StPlot()
    robot_info = forward_kinematics(*joints, **dh_params)
    pl.draw_robot(robot_info['robot_geometry_global'])
    pl.draw_tool(robot_info['flange'],
                 dh_params['tool'])
    pl.draw_frame(wobj, size=0.1)
    map(lambda x: pl.draw_frame(x, size=0.01), global_tars)
    pl.show()

