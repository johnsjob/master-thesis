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
##print; print "Init plots..."
##ax,_ = init_plot()
#----------------------------------------#
print; print "Define planes..."
transformed_paper = define_plane([1,1,1],[2,-2,-1],[1,3,0])
(origin_transformed_paper, basis_x_transformed_paper,
 basis_y_transformed_paper, normal_transformed_paper) = transformed_paper

untransformed_paper = define_plane([0,0,0],[1,0,0],[0,1,0])
(origin_paper, basis_x_paper, basis_y_paper, normal_paper) = untransformed_paper
#----------------------------------------#
print; print "Placing points in coordinate systems ("+str(num_points)+")..."
transformed_point = []
untransformed_point = []
for i in xrange(0, num_points):
    x, y = (rand(), rand())
    transformed_point.append( get_plane_point(transformed_paper,x,y) )
    untransformed_point.append( get_plane_point(untransformed_paper,x,y) )

transformed_point = array(transformed_point)
untransformed_point = array(untransformed_point)

transformed_directions = -diff(transformed_point.T).T
untransformed_directions = -diff(untransformed_point.T).T
untransformed_directions = untransformed_directions[:,0:2] #no z-coord
#----------------------------------------#
print; print "solve for basis in world coordinates from our sampled directions:"
if num_points == 3: #exact solution possible, (num_points - 1) * 3 direction eqs, (num_points * 2) unknowns
    print "-> solving exact..."
    plane_basis = linalg.solve(untransformed_directions,
                               transformed_directions)
elif num_points > 3: #overdetermined system, ((num_points - 1) - 2) * 3 direction eqs extra and (num_points - 3) * 2 unknowns extra
    print "-> solving with LMS..."
    plane_basis = linalg.solve(untransformed_directions.T.dot(untransformed_directions),
                               untransformed_directions.T.dot(transformed_directions))
else: #Underdetermines system
    raise Exception("Underdetermined system - too few points.")

plane_basis_error = numpy.abs(plane_basis[0] - basis_x_transformed_paper, plane_basis[1] - basis_y_transformed_paper)
print "transformed plane basis vectors error = " + str(plane_basis_error)
#----------------------------------------#
print; print "Solve for the transformed_paper origin in world coordinates 100 times"
print "using different points to make sure what point we choose doesn't matter:"
err = []
for _ in xrange(0,100):
    index = numpy.int32( rand()*(num_points-1) )
    calc_origin = transformed_point[index] - plane_basis[0]*untransformed_point[index][0] - plane_basis[1]*untransformed_point[index][1]
    err.append(origin_transformed_paper - calc_origin)
print "origin error: " + str(numpy.max(numpy.abs(err),0))
#----------------------------------------#
def get_normalized_vector(*components):
    return mat(components) / norm(components)
###========================================#
#Part two of calibration algorithm
import denavit_hartenberg as RK
Xtcp = mat([0,0,0,1])

L = 0.1 #10cm
delta_prim = L*get_normalized_vector(rand(), rand(), rand())
delta_prim = mat(list(delta_prim) + [1])

data = {'T':[], 'Xtcp':[], 'Xprim':[]}
for _ in xrange(0,3):
    a = rand_range(-45, 45)
    b = rand_range(-45, 45)
    c = rand_range(-45, 45)
    d = rand_range(-45, 45)
    e = rand_range(-45, 45)
    f = rand_range(-45, 45)
    T44 = RK.calc_tool_IRB120(a,b,c,d,e,f)
    data['T'].append(T44)
    data['Xtcp'].append(T44.dot([0,0,0,1]))
    data['Xprim'].append(data['Xtcp'][-1][0:3] + delta_prim[0:3])
    data['Xprim'][-1] = mat(list(data['Xprim'][-1]) + [1])

