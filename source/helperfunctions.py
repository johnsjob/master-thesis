from __future__ import division
#----------------------------------------#
from scipy import *
from numpy import *
from pylab import *
from numpy.linalg import *
from numpy.random import rand
from numpy import array as mat
import sympy as S
#----------------------------------------#
def matmul_series(*matrix_factors):
    '''
        Takes a list of matrices as arguments and perform
        a series of matrix multiplications from left to right
        in the order given.

        The parameters may contain vectors and scalars as well as long as
        the matrix-vector dimensions are the same.

        Note:
        The sequence can not start with a scalar, it must start
        with either a vector or a matrix of type numpy.ndarray.
    '''
    return reduce(dot, matrix_factors, 1)
#----------------------------------------#
def matmul_series_symbolic(*matrix_factors):
    '''
        Takes a list of matrices as arguments and perform
        a series of matrix multiplications from left to right
        in the order given.

        The parameters may contain vectors and scalars as well as long as
        the matrix-vector dimensions are the same.

        Note:
        The sequence can not start with a scalar, it must start
        with either a vector or a matrix of type numpy.ndarray.
    '''
    return reduce(S.multiply, matrix_factors, 1)
#----------------------------------------#
def deg_to_rad(x):
    return x * pi / 180
#----------------------------------------#
def rotation_matrix_z(angle):
    '''
        creates a rotation Z-mapping from subspace to worldspace
        using euler angles in degrees
    '''
    angle = deg_to_rad(angle)
    return array([[ cos(angle), -sin(angle), 0 ],
                  [ sin(angle),  cos(angle), 0 ],
                  [  0   ,    0  , 1 ]])
#----------------------------------------#
def rotation_matrix_x(angle):
    '''
        creates a rotation X-mapping from subspace to worldspace
        using euler angles in degrees
    '''
    angle = deg_to_rad(angle)
    return array([[ 1,     0     ,           0 ],
                  [ 0, cos(angle), -sin(angle) ],
                  [ 0, sin(angle),  cos(angle)]])
#----------------------------------------#
def rotation_matrix_y(angle):
    '''
        creates a rotation Y-mapping from subspace to worldspace
        using euler angles in degrees
    '''
    angle = deg_to_rad(angle)
    return array([[ cos(angle),     0     ,   sin(angle) ],
                  [ 0,              1     ,   0 ],
                  [ -sin(angle),    0     ,   cos(angle)]])
#----------------------------------------#
def rotation_matrix_rot_tilt_skew(rot, tilt, skew):
    '''
        creates a rotation ZXZ-mapping from subspace to worldspace
        using euler angles in degrees (extrinsic mapping)
    '''
    return matmul_series(rotation_matrix_z(-rot), rotation_matrix_x(tilt), rotation_matrix_z(skew))
#----------------------------------------#
def rotation_matrix_skew_tilt_rot(rot, tilt, skew):
    '''
        creates a rotation ZXZ-mapping from subspace to worldspace
        using euler angles in degrees (intrinsic mapping)
    '''
    return matmul_series(rotation_matrix_z(skew), rotation_matrix_x(tilt), rotation_matrix_z(-rot))
#----------------------------------------#
def orthogonal_projection_vectors(a, b):
    '''
        Orthogonal projection of vector a onto vector b
    '''
    c = b / norm(b)
    return a.dot(c) * c
#----------------------------------------#
def gram_schmith_step(a, b):
    '''
        Calculates the orthogonal component of A relative to the
        parallell direction vector B
    '''
    return a - orthogonal_projection_vectors(a, b)
#----------------------------------------#
def coordinate_system_from_two_directions(dirz, diry):
    '''
        Generates an orthogonal coordinate system from two vectors in world space,
        the vectors themselves need not be orthogonal.

        Alternative interpretation:
        creates a (general) rotation mapping from subspace to worldspace
        using vector directions
    '''
    dirz = dirz / norm(dirz)
    diry = diry / norm(diry)
    diry = gram_schmith_step(diry, dirz)
    dirx = cross(diry, dirz)
    return array([dirx, diry, dirz]).T
#----------------------------------------#
def convert_to_mat(*args):
    """
        Converts a list of arguments to numpy.ndarray,
        the arguments are expected to be tuples or lists,
        scalar values are also accepted but it's not the functions
        intended use.
    """
    return map(mat, args)
#----------------------------------------#
def define_plane(origin, dirx, diry):
    """
        This function returns a coordinate system,
        defined as a set containing a point of origin
        and a set of basis vectors all of dimension N.
        {O,{basis_x, basis_y, ...}}
    """
    origin, dirx, diry = convert_to_mat(origin, dirx, diry)
    diry = gram_schmith_step(diry, dirx)
    diry = diry / norm(diry)
    dirx = dirx / norm(dirx)
    normal = cross(dirx, diry)
    return origin, dirx, diry, normal
#----------------------------------------#
def get_normalized_vector(*components):
    return mat(components) / norm(components)
#----------------------------------------#
def get_plane_point(plane,x,y):
    """
        From a subspace plane coordinate system with local coordinates (x,y)
        return the coordinates in corresponding world-space.
    """
    origin, dx, dy, n = plane
    pos = origin + x*dx + y*dy
    return pos
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
def plot_equal_perspective(ax, xList, yList, zList):
    xmin = xList[0]; xmax = xList[1]
    ymin = yList[0]; ymax = yList[1]
    zmin = zList[0]; zmax = zList[1]
    Xb = 0.5*abs(xmax-xmin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(xmax + xmin)
    Yb = 0.5*abs(ymax-ymin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(ymax + ymin)
    Zb = 0.5*abs(zmax-zmin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(zmax + zmin)
    # Comment or uncomment following both lines to test the fake bounding box:
##    for xb, yb, zb in zip(Xb, Yb, Zb):
##       ax.plot([xb], [yb], [zb], 'w')    
#----------------------------------------#
def _swap(a,b):
        tmp = b
        b = a
        a = tmp
        return a, b
#----------------------------------------#
def rand_range(low, high):
        if low > high:
                low, high = _swap(low, high)
        return low + rand()*(high-low)
#----------------------------------------#
