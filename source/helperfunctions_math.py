#numpy imports
import numpy

from numpy import zeros, diff, sqrt
from numpy import pi, cos as _cos, sin as _sin,\
                  arccos as _acos, arcsin as _asin,\
                  arctan2 as _atan2, arctan as _atan

from numpy import sign
from numpy import array as mat, dot, cross, inf
from numpy import linalg
from numpy.linalg import norm

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
    0:    R                                       (1x)
    1:    (rot, tilt, skew, x,y,z)                (3x, 3x)
    2:    (rot, tilt, skew, t)                    (3x, 1x)
    3:    (angles, x,y,z)                         (1x, 3x)
    3.1:  (R, x,y,z)                              (1x, 3x)
    4:    (R, t), where R is list or numpy.array  (1x, 1x)
    """
    l_args = len(args)
    valid_types  = [tuple, list, numpy.ndarray]
    
    ## TODO: separate into subcase-functions (!!!)
    if  l_args == 1: #case 0
        if type(args[0]) in [list, tuple]:
            R = mat(args[0])
        else:
            R = args[0]
        m,n = R.shape
        if m == n:
            R = expand_matrix(R)
            return R
        else:
            return R
    elif  l_args == 2: #case 4
        R,t = args
    elif l_args == 4:
        if type(args[0]) in valid_types:
            if type(args[0][0] ) in valid_types: #case 3.1
                R,t = args[0], args[1:]
            else:
                (rot, tilt, skew), t = args[0], args[1:] #case 3
        else: #case 2
            (rot, tilt, skew), t = args[:3], args[3:][0]
    elif l_args == 6: #case 1
        (rot, tilt, skew), t = args[:3], args[3:]
    
    if not type(t) in valid_types:
        raise ArithmeticError('translation part can not be a scalar.')
    elif len(t) < 3:
        raise ArithmeticError('translation part must be of dimension 3.')

    if not 'R' in locals():
        R = rotation_matrix_rot_tilt_skew(rot, tilt, skew)

    T = zeros((4,4))
    m,n = R.shape

    T[:m, :n] = R
    T[:3, 3] = t[:3]
    T[3, :] = [0, 0, 0, 1]
    return T
#----------------------------------------------------------------------------------------------------------#
def homogenous_matrices(frames):
    """
    Creates a matrix of homogenous matrices( 4x4 matrix of type [[R,t],[0,1]] ),
    each frame in the input 'frames' can have one of following descriptions:
    0:    R                                       (1x)
    1:    (rot, tilt, skew, x,y,z)                (3x, 3x)
    2:    (rot, tilt, skew, t)                    (3x, 1x)
    3:    (angles, x,y,z)                         (1x, 3x)
    3.1:  (R, x,y,z)                              (1x, 3x)
    4:    (R, t), where R is list or numpy.array  (1x, 1x)
    """
    return mat( map(lambda x: homogenous_matrix(*x), frames) )
#----------------------------------------------------------------------------------------------------------#
def rot_to_quat(R):
    """
    """

    r11, r12, r13,\
    r21, r22, r23,\
    r31, r32, r33 = R.reshape(9) 

    q0 = ( r11 + r22 + r33 + 1.0) / 4.0
    q1 = ( r11 - r22 - r33 + 1.0) / 4.0
    q2 = (-r11 + r22 - r33 + 1.0) / 4.0
    q3 = (-r11 - r22 + r33 + 1.0) / 4.0
    if(q0 < 0.0):
        q0 = 0.0
    if(q1 < 0.0):
        q1 = 0.0
    if(q2 < 0.0):
        q2 = 0.0
    if(q3 < 0.0):
        q3 = 0.0
    q0 = sqrt(q0)
    q1 = sqrt(q1)
    q2 = sqrt(q2)
    q3 = sqrt(q3)
    if((q0 >= q1) and (q0 >= q2) and (q0 >= q3)):
        q0 *= +1.0
        q1 *= sign(r32 - r23)
        q2 *= sign(r13 - r31)
        q3 *= sign(r21 - r12)
    elif((q1 >= q0) and (q1 >= q2) and (q1 >= q3)):
        q0 *= sign(r32 - r23)
        q1 *= +1.0
        q2 *= sign(r21 + r12)
        q3 *= sign(r13 + r31)
    elif((q2 >= q0) and (q2 >= q1) and (q2 >= q3)):
        q0 *= sign(r13 - r31);
        q1 *= sign(r21 + r12);
        q2 *= +1.0
        q3 *= sign(r32 + r23);
    elif(q3 >= q0 and q3 >= q1 and q3 >= q2):
        q0 *= sign(r21 - r12);
        q1 *= sign(r31 + r13);
        q2 *= sign(r32 + r23);
        q3 *= +1.0
    else:
        print("coding error\n")

    q = mat( (q0, q1, q2, q3) )
    q = q / norm(q)

    return q
#----------------------------------------------------------------------------------------------------------#
def quat_to_rot(q):
    """
    """
    q0, q1, q2, q3 = q

    r1 = [
        q0**2 + q1**2 - q2**2 - q3**2,
        2*(q1*q2 + q0*q3),
        2*(q1*q3 - q0*q2)
        ]

    r2 = [
        2*(q1*q2 - q0*q3),
        q0**2 - q1**2 + q2**2 - q3**2,
        2*(q2*q3 + q0*q1)
        ]

    r3 = [
        2*(q1*q3 + q0*q2),
        2*(q2*q3 - q0*q1),
        q0**2 - q1**2 - q2**2 + q3**2
        ]

    R = mat(zip(r1,r2,r3))
    return R
#----------------------------------------------------------------------------------------------------------#
def apply_transform_on_frames(T44, frames):
    return mat( map(lambda x: matmul(T44, x), frames) )
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
        raise ArithmeticError('Only positive values (or zero) allowed for changing the dimension of the matrix.')
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
def rotation_matrices(rts):
    '''
        creates a rotation ZXZ-mapping from subspace to worldspace
        using euler angles in degrees (intrinsic mapping)

        (!) Note - rot is rotated counter-clockwise compared to ordinary euler z-rotation.
    '''
    M = lambda x: homogenous_matrix(rotation_matrix_rot_tilt_skew(*x))
    return mat( map(M, rts) )
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
#----------------------------------------#
def rot_tilt_skew(R):
    zs_vector = R[2,:] #matmul(R.T, numpy.array([0, 0,  1]) )
    zp_vector = -R[:,2] #matmul(R , numpy.array([0, 0, -1]) )

    #angle between normal and zp_vector direction
    #rot  = atan2( zp_vector[0], zp_vector[1] ) * 180 / pi
    rot  = atan2( zp_vector[0], zp_vector[1]) * 180 / pi
    tilt = -atan( norm(zp_vector[0:2]) / zp_vector[2] ) * 180 / pi
    skew = -atan2( zs_vector[0], zs_vector[1] ) * 180 / pi

    if tilt < 0.01:
        x_vector = R[:,0]
        z_vector = R[:,2]
        rot = atan2(x_vector[0], x_vector[1]) * 180 / pi
    return (rot, tilt, skew)
