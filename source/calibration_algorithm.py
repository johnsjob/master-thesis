from __future__ import division
from helperfunctions import *
import numpy
#----------------------------------------#
num_points = 100
#----------------------------------------#
def get_plane_point(plane,x,y):
    """
        From a subspace plane coordinate system with local coordinates (x,y)
        return the coordinates in corresponding world-space.
    """
    origin, basis_x, basis_y, _ = plane
    pos = origin + x*basis_x+ y*basis_y
    return pos
#========================================#
print; print "Init plots..."
ax,_ = init_plot()
#----------------------------------------#
print; print "Define planes..."
transformed_paper = define_plane([1,1,1],[2,-2,-1],[1,3,0])
(origin_transformed_paper, basis_x_transformed_paper,
 basis_y_transformed_paper, normal_transformed_paper) = transformed_paper

untransformed_paper = define_plane([0,0,0],[1,0,0],[0,1,0])
(origin_paper, basis_x_paper, basis_y_paper, normal_paper) = untransformed_paper
#----------------------------------------#
print; print "Placing points in coordinate systems.."
transformed_point = []
untransformed_point = []
for i in xrange(0, num_points):
    x, y = (rand(), rand())
    transformed_point.append( get_plane_point(transformed_paper,x,y) )
    untransformed_point.append( get_plane_point(untransformed_paper,x,y) )

transformed_point = array(transformed_point)
untransformed_point = array(untransformed_point)
###========================================#
##  picks directions in both plane coordinates,
##  depecated code only tested 3 points
##index = numpy.int32(rand()*(num_points-3))
##print "index picked: " + str(index)
##transformed_directions = transformed_point[index] - transformed_point[index+1], transformed_point[index+1] - transformed_point[index+2]
##untransformed_directions = untransformed_point[index] - untransformed_point[index+1], untransformed_point[index+1] - untransformed_point[index+2]
###========================================#
transformed_directions = -diff(transformed_point.T).T
untransformed_directions = -diff(untransformed_point.T).T
untransformed_directions = untransformed_directions[:,0:2] #no z-coord
#----------------------------------------#
print; print "solve for basis in world coordinates from our sampled directions:"
if num_points == 3: #exact solution possible, (num_points - 1) * 3 direction eqs, (num_points * 2) unknowns
    print "-> solving exact..."
    plane_basis = linalg.solve(untransformed_directions, transformed_directions)
elif num_points > 3: #overdetermined system, ((num_points - 1) - 2) * 3 direction eqs extra and (num_points - 3) * 2 unknowns extra
    print "-> solving with LMS..."
    plane_basis = linalg.solve(untransformed_directions.T.dot(untransformed_directions), untransformed_directions.T.dot(transformed_directions))
else: #Underdetermines system
    raise Exception("Underdetermined system - too few points.")

plane_basis_error = numpy.abs(plane_basis[0] - basis_x_transformed_paper, plane_basis[1] - basis_y_transformed_paper)
print "plane basis vectors error = " + str(plane_basis_error)
#----------------------------------------#
print; print "Solve for the transformed_paper origin in world coordinates:"
p_pos = transformed_point[0], transformed_point[1], transformed_point[2]
p_paper_pos = untransformed_point[0], untransformed_point[1], untransformed_point[2]
o0 = p_pos[0] - plane_basis[0]*p_paper_pos[0][0] - plane_basis[1]*p_paper_pos[0][1]
o1 = p_pos[1] - plane_basis[0]*p_paper_pos[1][0] - plane_basis[1]*p_paper_pos[1][1]
o2 = p_pos[2] - plane_basis[0]*p_paper_pos[2][0] - plane_basis[1]*p_paper_pos[2][1]
print origin_transformed_paper
print o0
print o1
print o2
###========================================#
##This part tests the theory of the algorithm
##in the case of a flat surface where it works even when varying orientation
##(I did not expect this, but the power of LS-methods are stronger than I thought?
print "\nSecond stage of algorithm - find delta-prim"
import numpy
mat = numpy.array


Rk = []
dprim = mat([rand(), rand(), rand()])
dprim = 100*rand()*dprim / norm(dprim)
xprimk = []
Xtcpk = []
dtk = []
Tk = []

N = num_points
for k in xrange(0,N):
    Rk.append( rotation_matrix(rand()*360-180, rand()*180-90, rand()*360-180) )
    xprimk.append( mat([rand()*1000, rand()*1000, 0, 1]) )

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
dRk = Rk[xrange(1, N)] - Rk[xrange(0, N-1)]

dxprimk = xprimk[xrange(1, N)] - xprimk[xrange(0, N-1)]
##dxprimk[:,3] = 1

dXtcpk = Xtcpk[xrange(1, N)] - Xtcpk[xrange(0, N-1)]
##dXtcpk[:,3] = 1

ddtk = dtk[xrange(1, N)] - dtk[xrange(0, N-1)]

dTk = Tk[xrange(1, N)] - Tk[xrange(0, N-1)]
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

print "\normal_transformed_paper=== Result ==="
print '\tdelta-prim = ' + str(dprim)
print '\tnorm = ' + str(norm(dprim))
print '\normal_transformed_paper\ts     = ' + str(s)
print '\tnorm = ' + str(norm(s))
print '\normal_transformed_paper\terror = ' + str(numpy.abs(s[0:3]-dprim))


###========================================#
ax.scatter(transformed_point[:,0], transformed_point[:,1], transformed_point[:,2])
ax.scatter(untransformed_point[:,0], untransformed_point[:,1], untransformed_point[:,2])
plot_plane(ax, transformed_paper)
plot_plane(ax, untransformed_paper)
show()
