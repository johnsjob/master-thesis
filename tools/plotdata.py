# -*- coding: utf-8 -*-
from __future__ import division
#--------------------------------------------------------------------------#
from parsingtools import *
import pylab as P
import numpy as N
from fixedpoint import fpsg32 as q
from fixedpoint import fpGetMax,fpGetMin
#--------------------------------------------------------------------------#
def extract(ls,tol):
    ret = []
    index = []
    for i,x in enumerate(ls):
        if abs(x) > tol:
            ret.append(x)
        else:
            ret.append(0)
        index.append(i)
    return index,ret
#--------------------------------------------------------------------------#
def get_nonzero_ind(valerr):
    ret = []
    for i,x in enumerate(valerr):
        if abs(x) > 0:
            ret.append(i)
    return ret
#--------------------------------------------------------------------------#
def bar(data,keyanal,keymeas,valerr):
    P.bar(*get_diff_nonzero(data,keyanal,keymeas,valerr))
#--------------------------------------------------------------------------#
def plotall(data,keyanal,keymeas):
    P.plot(*get_diff_all(data,keyanal,keymeas))
    P.title('[fixpoint representation]')
    P.xlabel('Matrix-list index')
    P.ylabel('fixpoint value')
#--------------------------------------------------------------------------#
def plotlogall(data,keyanal,keymeas):
    P.plot(*get_logdiff_all(data,keyanal,keymeas))
#--------------------------------------------------------------------------#
def get_diff_nonzero(data,keyanal,keymeas,valerr):
    return get_nonzero_ind(valerr),abs(data[keyanal][get_nonzero_ind(valerr)]-data[keymeas][get_nonzero_ind(valerr)])
#--------------------------------------------------------------------------#
def get_diff_all(data,keyanal,keymeas):
    return range(0,len(data[keyanal])), abs(data[keyanal]-data[keymeas])
#--------------------------------------------------------------------------#
def get_logdiff_all(data,keyanal,keymeas):
    #return range(0,len(data[keyanal])), N.log10(abs(data[keyanal]-data[keymeas])/2**24)
    return range(0,len(data[keyanal])), N.log10(abs(data[keyanal]-data[keymeas]))
#--------------------------------------------------------------------------#
def get_var(data,key):
    return range(0,len(data[key])), data[key]
#--------------------------------------------------------------------------#
def plot_var(data,key):
    P.plot(*get_var(data,key))
#--------------------------------------------------------------------------#
data = {}
#--------------------------------------------------------------------------#
filepath = get_input_path("case0__theta_is_largerThan_45__skew_is_inRange_0_to_minus_90__rot_is_near_0")
filepath = get_input_path("case1__theta_is_largerThan_45__skew_is_inRange_0_to_90__rot_is_near_0")
filepath = get_input_path("case2__theta_is_largerThan_45__skew_is_near_0__rot_is_inRange_0_to_minus_90")
filepath = get_input_path("case3__theta_is_largerThan_45__skew_is_near_0__rot_is_inRange_0_to_90")
filepath = get_input_path("case4__theta_is_inRange_45_to_minus_45__skew_is_near_0__rot_is_near_0")
filepath = get_input_path("case5__theta_is_near_0__skew_is_near_0__rot_is_near_0")

filepath = get_input_path("case5__theta_is_near_0__skew_is_near_0__rot_is_near_0")

lines = process_file(filepath)
data['anal_K'] = N.array( get_float_from_field(lines,"K")   )
data['anal_P'] = N.array( get_float_from_field(lines,"P")   )
data['anal_Q'] = N.array( get_float_from_field(lines,"Q")   )
data['anal_RR'] = N.array( get_float_from_field(lines,"RR") )

data['anal_S'] = N.array( get_float_from_field(lines,"S")  )
data['anal_U'] = N.array( get_float_from_field(lines,"_U") )
data['anal_V'] = N.array( get_float_from_field(lines,"V")  )

data['anal_W'] = N.array( get_float_from_field(lines,"W") )
data['anal_X'] = N.array( get_float_from_field(lines,"X") )
data['anal_Y'] = N.array( get_float_from_field(lines,"Y") )
data['anal_Z'] = N.array( get_float_from_field(lines,"Z") )

data['anal_s'] = N.array( get_float_from_field(lines,"_s") )
data['anal_A'] = N.array( get_float_from_field(lines,"_A") )
data['anal_D'] = N.array( get_float_from_field(lines,"D")  )

res = N.array( get_float_from_field(lines,"xpt")  )
data['anal_xpt'] = res[:,0]
data['anal_ypt'] = res[:,1]

data['anal_Mt'] = N.array( get_float_from_field(lines,"Mt") )
data['anal_t'] = N.array( get_float_from_field(lines,"_t") )

data['anal_R'] = N.array( get_all_R(lines) )
data['anal_T'] = N.array( get_all_T(lines) )

t11 = data['anal_T'][:,0,0]
t12 = data['anal_T'][:,0,1]
t13 = data['anal_T'][:,0,2]
t21 = data['anal_T'][:,1,0]
t22 = data['anal_T'][:,1,1]
t23 = data['anal_T'][:,1,2]
t31 = data['anal_T'][:,2,0]
t32 = data['anal_T'][:,2,1]
t33 = data['anal_T'][:,2,2]
#--------------------------------------------------------------------------#
filepath = "meas.txt"
lines = process_file(filepath)

data['meas_K'] = N.array( get_float_from_field(lines,"K") )
data['meas_P'] = N.array( get_float_from_field(lines,"P") )
data['meas_Q'] = N.array( get_float_from_field(lines,"Q") )
data['meas_RR'] = N.array( get_float_from_field(lines,"_RR") )

data['meas_S'] = N.array( get_float_from_field(lines,"S") )
data['meas_U'] = N.array( get_float_from_field(lines,"U") )
data['meas_V'] = N.array( get_float_from_field(lines,"V") )

data['meas_W'] = N.array( get_float_from_field(lines,"W") )
data['meas_X'] = N.array( get_float_from_field(lines,"X") )
data['meas_Y'] = N.array( get_float_from_field(lines,"Y") )
data['meas_Z'] = N.array( get_float_from_field(lines,"Z") )

data['meas_s'] = N.array( get_float_from_field(lines,"_s") )
data['meas_staljare'] = N.array( get_float_from_field(lines,"staljare") )
data['meas_snamnare'] = N.array( get_float_from_field(lines,"snamnare") )
data['meas_A'] = N.array( get_float_from_field(lines,"_A") )
data['meas_D'] = N.array( get_float_from_field(lines,"D") )
data['meas_xpt'] = N.array( get_float_from_field(lines,"_xpt") )
data['meas_ypt'] = N.array( get_float_from_field(lines,"_ypt") )

data['meas_Mt'] = N.array( get_float_from_field(lines,"_Mt") )
data['meas_t'] = N.array( get_float_from_field(lines,"_t") )
data['meas_Mttaljare'] = N.array( get_float_from_field(lines,"Mttaljare") )
data['meas_Mtnamnare'] = N.array( get_float_from_field(lines,"Mtnamnare") )
data['meas_R'] = N.array(N.array( get_all_R_int(lines) ) / 2**24.0, N.float32)
######################################################################
data['tol'] = []
data['tol9'] = []
data['R_err'] = []
data['origin'] = []
data['err'] = []
#--------------------------------------------------------------------------#
tol = 2e-3
tol9 = 4.5*2e-3
for i,e in enumerate(data['anal_R']):
    R_meas = data['meas_R']
    #data['R_err'].append(N.linalg.norm( e.reshape(9) - R_meas[i].reshape(9) ))
    data['R_err'].append(N.sum(N.abs( e.reshape(9) - R_meas[i].reshape(9)) ))
    data['err'].append(N.abs( e.reshape(9) - R_meas[i].reshape(9)))
    #data['tol'].append( q(N.sqrt(2e-3**2+2e-3**2+2e-3**2)) )
    data['tol'].append( tol )
    data['tol9'].append( tol9 )
    data['origin'].append(0)
#ind,valerr = extract(data['R_err'],2000)
ind,valerr = extract(data['R_err'],tol9)
tol = N.array(data['tol'])
#terr = abs(data['anal_t']-data['meas_t'])
#--------------------------------------------------------------------------#
err = abs(data['meas_R'] - data['anal_R'])
x,y,z = err.shape
err2 = err.reshape((x,y*z)).reshape(x*y*z)
errsum = [sum(sum(k)) for k in err]
tol = [tol[0] for k in range(0,x*y*z)]
r11 = err[:,0,0]
r12 = err[:,0,1]
r13 = err[:,0,2]
r21 = err[:,1,0]
r22 = err[:,1,1]
r23 = err[:,1,2]
r31 = err[:,2,0]
r32 = err[:,2,1]
r33 = err[:,2,2]
###--------------------------------------------------------------------------#
x0 = 50
x1 = 106
###--------------------------------------------------------------------------#
##P.plot(data['origin'],"k--",linewidth=3)
##P.plot(data['tol9'],"k--")
##P.plot(data['R_err'])
##P.plot(ind,valerr,'r')
##P.plot(terr,'g--')
##P.title('Elementwise L1 norm of the element errors - TOL = %i ( %0.5f... )' % (tol9, tol9 / 2**24.0))
##P.xlabel('Matrix-list index')
##P.ylabel('fixpoint value')
##P.plt.xlim(0,x)
###P.plt.xlim(x0,x1)
##P.show()
###--------------------------------------------------------------------------#
P.plot(err2)
P.plot(tol,"k--")
P.title('Element errors for all elements for all matrices')
P.xlabel('Matrix element index')
P.ylabel('float value')
P.plt.xlim(0,x*y*z)
#P.plt.xlim(x0*y*z,x1*y*z)
P.show()
###--------------------------------------------------------------------------#
######print "Tolerance: %i" % tol[0]
######var = "snamnare"
##plotall(data,"anal_Mt","meas_Mt")
######plotall(data,"anal_s","meas_s")
#######plot_var(data,"meas_s")
######plot_var(data,"meas_Mttaljare")
######plot_var(data,"meas_Mtnamnare")
##P.plot(tol,"k--")
##P.legend(["Mt error",'tol'])
##P.plt.xlim(0,x)
##P.plt.xlim(x0,x1)
##P.show()
###--------------------------------------------------------------------------#
##plotall(data,"anal_K","meas_K")
##plotall(data,"anal_P","meas_P")
##plotall(data,"anal_Q","meas_Q")
##plotall(data,"anal_RR","meas_RR")
##P.plot(tol,"k--")
##P.legend(['K','P','Q','RR','tol'])
##P.plt.xlim(0,x)
##P.plt.xlim(x0,x1)
##P.show()
##
##plotall(data,"anal_S","meas_S")
##plotall(data,"anal_U","meas_U")
##plotall(data,"anal_V","meas_V")
##P.plot(tol,"k--")
##P.legend(['S','U','V','tol'])
##P.plt.xlim(0,x)
##P.plt.xlim(x0,x1)
##P.show()
##
##plotall(data,"anal_W","meas_W")
##plotall(data,"anal_X","meas_X")
##plotall(data,"anal_Y","meas_Y")
##plotall(data,"anal_Z","meas_Z")
##P.plot(tol,"k--")
##P.legend(['W','X','Y','Z','tol'])
##P.plt.xlim(0,x)
##P.plt.xlim(x0,x1)
##P.show()
##
##plotall(data,"anal_s","meas_s")
##plotall(data,"anal_A","meas_A")
##plotall(data,"anal_D","meas_D")
##P.plot(tol,"k--")
##P.legend(['s','A','D','tol'])
##P.plt.xlim(0,x)
##P.plt.xlim(x0,x1)
##P.show()
#--------------------------------------------------------------------------#
##Ranal = data['anal_R'][50:110]
##Rmeas = data['meas_R'][50:110]
##Rerr = abs(Ranal - Rmeas)
##leg = []
##log = numpy.log
##t = []
##r = []
##tee = []
##for i,x in enumerate(Rerr):
##    if (x > tol[0]).any():
##        print "ERROR INDEX: %i" % i
##        d = data['anal_T'][i].reshape(9)
##        #d -= -45788159
##        d = numpy.array(d,dtype='float32')
##        d /= 16777216
##        d *= 89929
##        P.plot(x.reshape(9),'g')
##        P.plot(d,'-.')
##        P.plt.xlim(0,8)
##        leg.append(str(i))
##        
##        tee.append(terr[i])
##        t.append(data['anal_T'][i].reshape(9))
##        r.append(x.reshape(9))
##P.plot(tol,"k--")
##P.plot(data['origin'],'k--',linewidth=2)
##leg.append('tol')
##P.legend(leg)
##P.ylabel("Fixpoint value")
##P.xlabel("Matrix element index")
##P.title("Case 1 - Matrix indexes with most errors")
##P.show()
