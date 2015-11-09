#numpy imports
from numpy import linalg
from numpy import zeros, diff
from numpy import pi, cos as _cos, sin as _sin,\
                  arccos as _acos, arcsin as _asin,\
                  arctan2 as _atan2, arctan as _atan

from numpy import array as mat, dot, cross, inf
from numpy.linalg import norm
import numpy
import numpy as n

#random import
from random import random as rand
#----------------------------------------------------------------------------------------------------------#
def deg_to_rad(x):
    return x * pi / 180.0
#----------------------------------------#
def rad_to_deg(x):
    return x / pi * 180.0
#----------------------------------------#
def atan(x, unit='rad'):
    if unit == 'rad':
        return _atan(x)
    elif unit == 'deg':
        return rad_to_deg(_atan(x))
#----------------------------------------#
def atan2(y, x, unit='rad'):
    if unit == 'rad':
        return _atan2(y,x)
    elif unit == 'deg':
        return rad_to_deg(_atan2(y,x))
#----------------------------------------#
def cos(x, unit='rad'):
    if unit == 'rad':
        return _cos(x)
    elif unit == 'deg':
        return _cos( deg_to_rad(x))
#----------------------------------------------------------------------------------------------------------#
def sin(x, unit='rad'):
    if unit == 'rad':
        return _sin(x)
    elif unit == 'deg':
        return _sin( deg_to_rad(x))
#----------------------------------------#
def acos(x, unit='rad'):
    if unit == 'rad':
        return _acos(x)
    elif unit == 'deg':
        return rad_to_deg( _acos(x) )
#----------------------------------------------------------------------------------------------------------#
def asin(x, unit='rad'):
    if unit == 'rad':
        return _asin(x)
    elif unit == 'deg':
        return rad_to_deg( _asin(x) )
#----------------------------------------------------------------------------------------------------------#
def matmul(*matrix_factors):
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
#----------------------------------------------------------------------------------------------------------#
def matmul_series(*matrix_factors):
    '''
        Takes a list of matrices as arguments and perform
        a series of matrix multiplications from left to right
        in the order given, but contrary to matmul
        it returns each partial result in a list.

        The parameters may contain vectors and scalars as well as long as
        the matrix-vector dimensions are the same.

        Note:
        The sequence can not start with a scalar, it must start
        with either a vector or a matrix of type numpy.ndarray.
    '''
    res = []
    res.append(matrix_factors[0])
    for i in xrange(1,len(matrix_factors)):
        res.append(reduce(dot, matrix_factors[:i+1], 1))
    return res
#----------------------------------------------------------------------------------------------------------#
def homogenous_matrix(*args):
    """
    Creates a homogenous matrix ( 4x4 matrix of type [[R,t],[0,1]] ),
    allows input on the forms:
    1: (rot, tilt, skew, x,y,z)
    2: (rot, tilt, skew, t)
    3: (angles, x,y,z)
    4: (R, t), where R is list or numpy.array
    """
    l_args = len(args)

    if  l_args == 2:
        R,t = args
    elif l_args == 4: 
        if type(args[0]) in (tuple, list, numpy.array()):
            (rot, tilt, skew), t = args[0], args[1:]
        else:
            (rot, tilt, skew), t = args[:3], args[3:][0]
    elif l_args == 6:
        (rot, tilt, skew), t = args[:3], args[3:]        

    if not type(t) in [list, numpy.array]:
        raise ArithmeticError('translation part can not be a scalar.')
    elif len(t) < 3:
        raise ArithmeticError('translation part must be of dimension 3.')
    

    R = rotation_matrix_rot_tilt_skew(rot, tilt, skew)
    T = zeros((4,4))
    T[:3, :3] = R
    T[:3, 3] = t
    T[3, :] = [0, 0, 0, 1]
    return T
#----------------------------------------------------------------------------------------------------------#
def homogenous_translation_x( tx ):
    return mat([[1,        0,      0,     tx],
                [0,        1,      0,      0],
                [0,        0,      1,      0],
                [0,        0,      0,      1]])
#----------------------------------------------------------------------------------------------------------#
def homogenous_translation_z( tz ):
    return mat([[1,     0,     0,      0],
                [0,     1,     0,      0],
                [0,     0,     1,      tz],
                [0,     0,     0,      1]])
#----------------------------------------------------------------------------------------------------------#
def homogenous_rotation_x( ang ):
    #convert to radians
    c = cos(ang * pi / 180)
    s = sin(ang * pi / 180)
    return mat([[1,     0,      0,      0],
                [0,     c,     -s,     0],
                [0,     s,     c,      0],
                [0,     0,     0,      1]])
#----------------------------------------------------------------------------------------------------------#
def homogenous_rotation_z( ang ):
    #convert to radians
    c = cos(ang * pi / 180)
    s = sin(ang * pi / 180)
    return mat([[c,     -s,    0,     0],
                [s,     c,     0,     0],
                [0,     0,     1,     0],
                [0,     0,     0,     1]])
#----------------------------------------#
def expand_matrix(matrix, k=1,l=1):
    if k < 0 or l < 0:
        raise ArithmeticError('Only positive values allowed for changing the dimension of the matrix.')
    m,n = matrix.shape
    res = numpy.eye(m+k,n+l)
    res[0:m, 0:n] = matrix
    return res
#----------------------------------------#    
def rotation_matrix_z(angle):
    '''
        creates a rotation Z-mapping from subspace to worldspace
        using euler angles in degrees
    '''
    angle = deg_to_rad(angle)
    return mat([[ cos(angle), -sin(angle), 0 ],
                  [ sin(angle),  cos(angle), 0 ],
                  [  0   ,    0  , 1 ]])
#----------------------------------------#
def rotation_matrix_x(angle):
    '''
        creates a rotation X-mapping from subspace to worldspace
        using euler angles in degrees
    '''
    angle = deg_to_rad(angle)
    return mat([[ 1,     0     ,           0 ],
                  [ 0, cos(angle), -sin(angle) ],
                  [ 0, sin(angle),  cos(angle)]])
#----------------------------------------#
def rotation_matrix_y(angle):
    '''
        creates a rotation Y-mapping from subspace to worldspace
        using euler angles in degrees
    '''
    angle = deg_to_rad(angle)
    return mat([[ cos(angle),     0     ,   sin(angle) ],
                  [ 0,              1     ,   0 ],
                  [ -sin(angle),    0     ,   cos(angle)]])
#----------------------------------------#
def rotation_matrix_rot_tilt_skew(rot, tilt, skew):
    '''
        creates a rotation ZXZ-mapping from subspace to worldspace
        using euler angles in degrees (extrinsic mapping)

        (!) Note - rot is rotated counter-clockwise compared to ordinary euler z-rotation.
    '''
    return matmul(rotation_matrix_z(-rot), rotation_matrix_x(tilt), rotation_matrix_z(skew))
#----------------------------------------#
def rotation_matrix_skew_tilt_rot(rot, tilt, skew):
    '''
        creates a rotation ZXZ-mapping from subspace to worldspace
        using euler angles in degrees (intrinsic mapping)

        (!) Note - rot is rotated counter-clockwise compared to ordinary euler z-rotation.
    '''
    return matmul(rotation_matrix_z(skew), rotation_matrix_x(tilt), rotation_matrix_z(-rot))
#----------------------------------------#
def coordinate_system_from_two_directions(dirz, diry):
    '''
        Generates an orthogonal coordinate system from two vectors in world space,
        the vectors themselves need not be orthogonal.

        Alternative interpretation:
        creates a (general) rotation mapping from subspace to worldspace
        using vector directions in world-space
    '''
    dirz = dirz / norm(dirz)
    diry = diry / norm(diry)
    diry = gram_schmith_step(diry, dirz)
    dirx = cross(diry, dirz)
    return mat([dirx, diry, dirz]).T
#----------------------------------------#
def get_normalized_vector(*components):
    return mat(components) / norm(components)
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
        Calculates with the help of vector A the orthogonal component relative to the
        parallell direction vector B
    '''
    return a - orthogonal_projection_vectors(a, b)
#----------------------------------------#
##########################################
#----------------------------------------#
def rand_range(low, high):
        if low > high:
                low, high = high, low
        return low + rand()*(high-low)
