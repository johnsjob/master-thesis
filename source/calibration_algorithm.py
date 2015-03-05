from __future__ import division
#----------------------------------------#
from scipy import *
from numpy import *
from pylab import *
from numpy.linalg import *
from numpy.random import rand
#----------------------------------------#
def rad(x):
    return x * pi / 180
#----------------------------------------#
def Rz(angle):
    angle = rad(angle)
    return array([[ cos(angle), -sin(angle), 0 ],
                        [ sin(angle),  cos(angle), 0 ],
                        [  0   ,    0  , 1 ]])
#----------------------------------------#
def Rx(angle):
    angle = rad(angle)
    return array([[ 1,     0     ,           0 ],
                  [ 0, cos(angle), -sin(angle) ],
                  [ 0, sin(angle),  cos(angle)]])
#----------------------------------------#
'''
    creates a rotation mapping from subspace to worldspace
    using euler angles in degrees
'''
def R(rot, tilt, skew):
    return dot(Rx(rot), dot(Rz(tilt), Rx(skew)))
#----------------------------------------#
'''
    creates a rotation mapping from subspace to worldspace
    using vector directions
'''
def Rv(dirz, diry):
    dirz = dirz / norm(dirz)
    diry = diry / norm(diry)
    diry = grahm_schm(diry, dirz)
    dirx = cross(diry, dirz)
    return array([dirx, diry, dirz]).T
#----------------------------------------#
'''
    Orthogonal projection of vector a onto vector b
'''
def P(a, b):
    c = b / norm(b)
    return a.dot(c) * c
#----------------------------------------#
'''
    Calculates the orthogonal component of A relative to the
    parallell direction vector B
'''
def grahm_schm(a, b):
    return a - P(a, b)
#----------------------------------------#
def define_plane(origin, dirx, diry):
    origin = array(origin); dirx = array(dirx); diry = array(diry)
    normal = cross(dirx, diry)
    normal = normal / norm(normal)
    diry = grahm_schm(diry, dirx)
    diry = diry / norm(diry)
    dirx = dirx / norm(dirx)
    return [origin, dirx, diry, normal]
#----------------------------------------#
def init_plot():
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca(projection='3d') 
    return ax, fig
#----------------------------------------#
def plot_plane(ax, plane):
    dx = array([get_plane_point(plane,0,0), get_plane_point(plane,1,0)])
    dy = array([get_plane_point(plane,0,0), get_plane_point(plane,0,1)])
    dz = array([get_plane_point(plane,0,0), get_plane_point(plane,0,0) + plane[3]])
    ax.plot(dx[:,0],dx[:,1],dx[:,2],'b')
    ax.plot(dy[:,0],dy[:,1],dy[:,2],'g')
    ax.plot(dz[:,0],dz[:,1],dz[:,2],'r')
#----------------------------------------#
def get_plane_point(plane,x,y):
    origin, dx, dy, n = plane
    pos = origin + x*dx + y*dy
    return pos
#----------------------------------------#
def get_relative_point(plane,x,y,Dx,Dy,L):
    origin, dx, dy, n = plane
    pos = get_plane_point(plane,x,y)
    vec = n + Dx*dx + Dy*dy
    vec = vec / norm(vec)
    rel = pos + L*vec
    dir = pos - rel
    dir = dir / norm(dir)
    return rel, dir, pos,(dx*x,dy*y)
#========================================#
print; print "Init plots..."
ax,_ = init_plot()
#========================================#
print; print "Define planes..."
pl = define_plane([1,1,1],[2,-2,-1],[1,3,0])
o, dx, dy, n = pl

paper = define_plane([0,0,0],[1,0,0],[0,1,0])
op, dxp, dyp, np = paper
#----------------------------------------#
print; print "Placing points in coordinate systems.."
p = []
p_paper = []
num_points = 4
for i in range(0, num_points):
    x = rand(); y = rand()
    p.append(get_plane_point(pl,x,y))
    p_paper.append(get_plane_point(paper,x,y))
p = array(p)
p_paper = array(p_paper)
###========================================#
#picks directions in both plane coordinates
p_dir = p[0] - p[1], p[1] - p[2]
p_paper_dir = p_paper[0] - p_paper[1], p_paper[1] - p_paper[2]
p_dir = array(p_dir)
p_paper_dir = array(p_paper_dir)
p_paper_dir = p_paper_dir[:,0:2] #no z-coord
#----------------------------------------#
print; print "solve for basis in world coordinates from our sampled directions:"
if num_points == 3: #exact solution possible, (num_points - 1) * 3 direction eqs, (num_points * 2) unknowns
    plane_basis = linalg.solve(p_paper_dir, p_dir)
elif num_points > 3: #overdetermined system, ((num_points - 1) - 2) * 3 direction eqs extra and (num_points - 3) * 2 unknowns extra
    plane_basis = linalg.solve(p_paper_dir.T.dot(p_paper_dir), p_paper_dir.T.dot(p_dir))
else: #Underdetermines system
    raise "Underdetermined system - too few points."
print plane_basis[0], plane_basis[1]
print dx, dy
#----------------------------------------#
print; print "Solve for the paper origin in world coordinates:"
p_pos = p[0], p[1], p[2]
p_paper_pos = p_paper[0], p_paper[1], p_paper[2]
o0 = p_pos[0] - plane_basis[0]*p_paper_pos[0][0] - plane_basis[1]*p_paper_pos[0][1]
o1 = p_pos[1] - plane_basis[0]*p_paper_pos[1][0] - plane_basis[1]*p_paper_pos[1][1]
o2 = p_pos[2] - plane_basis[0]*p_paper_pos[2][0] - plane_basis[1]*p_paper_pos[2][1]
print o
print o0
print o1
print o2
###========================================#
print "Experimental solution? I think this is what I have in the report"
import numpy
mat = numpy.array


Rk = []
dprim = mat([0.1, 0, 5])
xprimk = []
Xtcpk = []
dtk = []
Tk = []

N = 24
for k in range(0,N):
    Rk.append( R(rand()*45, rand()*45, rand()*45) )
    xprimk.append( mat([rand()*100, rand()*100, 0, 1]) )

    Xtcp = (xprimk[-1][0:3] + Rk[-1].dot(dprim)).tolist() + [1]
    Xtcpk.append( Xtcp )

    dtk.append( xprimk[-1] )

    T = zeros((4,4))
    T[0:3,0:3] = Rk[-1]
    T[:,  3] = dtk[-1]
    T[3, :] = [0, 0, 0, 1]
    Tk.append( T )    
Rk = mat(Rk)
xprimk = mat(xprimk)
Xtcpk = mat(Xtcpk)
dtk = mat(dtk)
Tk = mat(Tk)
###
dRk = Rk[range(1, N)] - Rk[range(0, N-1)]

dxprimk = xprimk[range(1, N)] - xprimk[range(0, N-1)]
##dxprimk[:,3] = 1

dXtcpk = Xtcpk[range(1, N)] - Xtcpk[range(0, N-1)]
##dXtcpk[:,3] = 1

ddtk = dtk[range(1, N)] - dtk[range(0, N-1)]

dTk = Tk[range(1, N)] - Tk[range(0, N-1)]
##dTk[:,3,3] = 1

tmp = dTk.reshape(((N-1)*4,4))
A = tmp
A = tmp.T.dot(A)
#A /= A[3,3]

B = dXtcpk - dxprimk
#B[:,3] = 1
B = B.reshape((N-1)*4)
B = tmp.T.dot(B)
#B /= B[3]

s = numpy.linalg.solve(A,B)
#s /= s[-1]
print s
print dprim
print norm(s)

print "All in all I remembr thinking these were all correct, but I will have to recheck my derivations and report."


###========================================#
##ax.scatter(p[:,0], p[:,1], p[:,2])
##ax.scatter(p_paper[:,0], p_paper[:,1], p_paper[:,2])
##plot_plane(ax, pl)
##plot_plane(ax, paper)
##show()
