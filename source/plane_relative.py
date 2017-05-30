from __future__ import division
#--------------------------#
import sys
#--------------------------#
import numpy as n
from numpy import cos, sin, pi, array
from helperfunctions_math import mat, nmap, nzip,\
     rotation_matrices as orientation_frames,\
     homogenous_matrices, matmul, rotation_matrix_rot_tilt_skew,\
     expand_matrix

from helperfunctions import convert_to_mat
#=====================================================#
#sys.path.append("../int/misc-tools/")
#--------------------------#
def get_relative_R(R0, R01):
    M = matmul(R0, R01)
    return M
#--------------------------#
def attach_to_base_frame(base_frame, *frames):
    return nmap(lambda x: matmul(base_frame, x), frames)
#----------------------------------------------------------------------------------------------------------#
def define_plane_from_directions(origin, dirx, diry, system_type = 'local'):
    """
        This function returns a coordinate system frame (NxN matrix),
        for a plane defined by a matrix containing a point of origin (t)
        and a set of basis vectors (R) all of dimension N: [[R, t,]
                                                            [0, 1]]
    """
    if not system_type in ['local','global']:
        raise ArithmeticError('Only local or global coordinate system descriptions allowed')
    origin, dirx, diry = convert_to_mat(origin, dirx, diry)
    diry = gram_schmith_step(diry, dirx)
    diry = diry / norm(diry)
    dirx = dirx / norm(dirx)
    normal = mat([rand(), rand(), rand()])
    normal = gram_schmith_step(normal, dirx)
    normal = gram_schmith_step(normal, diry)

    plane_transform = expand_matrix( mat([dirx, diry, normal,origin]).T,1,0)
    if system_type == 'global':
        return mat_flip(plane_transform)
    else:
        return plane_transform
#----------------------------------------#
def define_plane_from_angles(origin, r,t,s, system_type = 'local'):
    """
        This function returns a coordinate system frame (NxN matrix),
        for a plane defined by a matrix containing a point of origin (t)
        and a set of basis vectors (R) all of dimension N: [[R, t,]
                                                            [0, 1]]
    """
    if not system_type in ['local','global']:
        raise ArithmeticError('Only local or global coordinate system descriptions allowed')
    
    R = rotation_matrix_rot_tilt_skew(r,t,s)
    plane_transform = expand_matrix(R)
    plane_transform[:3,3] = origin
    plane_transform = expand_matrix( mat([R[:,0], R[:,1], R[:,2],origin]).T,1,0)

    if system_type == 'global':
        return mat_flip(plane_transform)
    else:
        return plane_transform
#----------------------------------------#
def define_plane_relative_from_angles(plane_transform, rel_origin, r,t,s, system_type='local'):
    """
    Calculates the plane orientation and position relative another plane,
    from local position and orientation.

    The reuslting plane inherits the coordinate system type from plane_transform
    and is either 'local' or 'global'
    """
    rel_plane = define_plane_from_angles(rel_origin, r, t, s, system_type)
    absolute_plane_transform = plane_transform.dot(rel_plane)
    return absolute_plane_transform
#----------------------------------------#
def define_plane_relative_from_plane(plane_transform, rel_plane):
    """
    Calculates the plane orientation and position relative another plane,
    from local position and orientation.

    The reuslting plane inherits the coordinate system type from plane_transform
    and is either 'local' or 'global'
    """
    absolute_plane_transform = plane_transform.dot(rel_plane)
    return absolute_plane_transform
#----------------------------------------#
def get_plane_point(plane_transform,x,y):
    """
        From a subspace plane coordinate system with local coordinates (x,y)
        return the coordinates in corresponding world-space.
    """
    pos = plane_transform.dot([x,y,0,1])
    return pos
#----------------------------------------#
def attach_frames(*args):
    return homogenous_matrices(zip(*args))
    
def place_curve(point_matrix, frame):
    point_matrix_tf = get_transformed_points(frame, point_matrix)
    return point_matrix_tf

def create_circle_curve(diameter=0.3, num_p=50):
    point_matrix = generate_symmetric_curve(num_points=num_p,
                                            ampl_factor=diameter)
    return point_matrix
#----------------------------------------#
def generate_symmetric_curve(t=None, x_func=n.cos, y_func=n.sin,
                   curve_factor = 2*n.pi, freq = 1, ampl_factor=0.1, num_points=50, offset=0):
    """
    Genrates a Homogenous point-curve with z=0 in the following matrix form:
    [[x1, y1, 0, 1},
         .....     ,
     [xn, yn, 0, 1}]
     in local (unransformed) plane-coordinates in metres.

    parameters:
                curve_factor : ...
                freq         : ...
                ampl_factor  : 'diameter' of the curve
                num_points   : ...
    """
    if t is None:
        t = n.linspace(0, 1, num_points)

    le = len(t)
    y = y_func(t/t[-1]*curve_factor*freq)
    y = y / n.max( n.abs(y) )

    x = x_func(t/t[-1]*curve_factor*freq)
    x = x / n.max( n.abs(x) )

    ampl_factor = ampl_factor/2.0
    point_matrix = mat(zip(x*ampl_factor, y*ampl_factor,n.zeros(le), n.ones(le)))
    return point_matrix
#----------------------------------------#
def generate_curve(xmin=-0.25, xmax=0.25, x_func=None, y_func=n.sin,
                   curve_factor = n.pi, freq = 1, ampl_factor=0.25, num_points=50, offset=0):
    """
    Genrates a Homogenous point-curve with z=0 in the following matrix form:
    [[x1, y1, 0, 1},
         .....     ,
     [xn, yn, 0, 1}]
     in local (unransformed) plane-coordinates in metres.
    """
    if x_func is None:
        if xmin > xmax:
            xmin, xmax = xmax, xmin
        x_func = n.linspace(xmin, xmax, num_points)

    le = len(x_func)
    y = y_func(x_func / n.abs(xmax) * curve_factor * freq)
    
    if n.max( n.abs(y) ) > 0:
        y = y / n.max( n.abs(y) )
    
    point_matrix = mat(zip(x_func + offset,
                           y * ampl_factor,
                           n.zeros(le),
                           n.ones(le)))
    return point_matrix
#----------------------------------------#
def get_transformed_points(plane_transform, point_matrix):
    """
    Transforms the local plane coordinates to global cordinates
    by performing <[m.T 1], [[R.T t],[0 1]]> where m is the untransformed point-matrix.
    """
    transf_points = point_matrix.dot( plane_transform.T )
    return transf_points

#--------------------------#
def mat_flip(M):
    transf_flip = mat([[1, 0, 0,  0],
                       [0, -1, 0, 0],
                       [0, 0, -1, 0],
                       [0, 0, 0,  1]])
    return transf_flip.dot(M)
#--------------------------#
def mat_flip3(M):
    transf_flip = mat([[1, 0, 0],
                       [0, -1, 0],
                       [0, 0, -1]])
    return transf_flip.dot(M)
#--------------------------#
