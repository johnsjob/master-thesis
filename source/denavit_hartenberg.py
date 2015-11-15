from helperfunctions_math import *
#----------------------------------------------------------------------------------------------------------#
rad = lambda x: x * pi / 180.0
deg = lambda x: x * 180.0 / pi
cos2 = lambda x: n.cos(rad(x))
sin2 = lambda x: n.sin(rad(x))

cos_sats = lambda a,b,th: a**2 + b**2 - 2*a*b*cos(rad(th)); #ok
ang_sats = lambda c,a,b: deg(acos((c**2 - a**2 - b**2)/(-2*a*b))); #ok
ang_sats2 = lambda c,a,b: deg(acos((c**2 - a**2 - b**2)/(2*a*b))); #ok
round = lambda x: custom_round(x)
atan = lambda x: deg(n.arctan(x))
atan2 = lambda y,x: deg(n.arctan2(y,x))

up_to = lambda i: custom_round(matmul(*[debug[x] for x in range(i)]))
#----------------------------------------------------------------------------------------------------------#
def transform_to_next(A, alpha, D, theta, convention='standard'):
    """
    Calculats transform from frame J = I-1 to I
    in order AI = Rzj( theta )Tzj( D )Txi( A )Rxi( alpha )
    and parameters are given in order A, alpha, D, theta, 
    according to  the standard convention:

            matrix =   mat([[ct,    -st*ca,  st*sa, A*ct],
                            [st,     ct*ca, -ct*sa, A*st],
                            [ 0,        sa,     ca,   D ],
                            [ 0,         0,      0,   1 ]]);

    If convention is changed to modified then the transformation
    is remapped according to the equivalent modified Denivit-hartenberg,
    and performs the mapping from I to I+1.
    """
    Rz_J = homogenous_rotation_z(theta)
    Tz_J = homogenous_translation_z(D)
    Tx_I = homogenous_translation_x(A)
    Rx_I = homogenous_rotation_x(alpha)
    matrix = matmul(Rz_J, Tz_J, Tx_I, Rx_I)
    
    if convention.lower() == 'standard':
        return matrix
    elif convention.lower() == 'modified':
        modified = mat([ matrix[0,0],  -matrix[1,0],   matrix[2,0],  A,
                        -matrix[0,1],   matrix[1,1],  -matrix[2,1], -matrix[2,1]*D,
                         matrix[0,2],  -matrix[1,2],   matrix[2,2],  matrix[2,2]*D,
                               0    ,         0    ,         0    ,  1,]).reshape((4,4))
        return modified
    else:
        raise ArithmeticError("As of writing this function only two conventions are allowed: 'standard' or 'modified'.")
#----------------------------------------------------------------------------------------------------------#
def forward_kinematics(*joint_values,**kwargs):
    """
    Performs the denavit-hartenberg algorithm
    and calculates T44 = A0 * A1 * A2 * ... * An multiplication.

    The parameters are entered in the order A, alpha, D, theta repeatedly
    if no order is specified. However, is another order specified then
    the values must be entered in that order.
    """

    # check so that only supported parameters
    # exist in kwargs and abort if any other parameter
    # had been created by mistake due to typos or similar
    input_param_diff = set(kwargs) - set(['order','convention','unit','table'])
    if len(input_param_diff) > 0:
        raise Exception('Unsupported parameters: ' + str(*input_param_diff))
    # handle which unit is being used
    if not kwargs.has_key('unit'):
        kwargs['unit'] = 'metre'
    elif kwargs['unit'] == 'm':
        kwargs['unit'] = 'metre'
    elif kwargs['unit'] == 'mm':
        kwargs['unit'] = 'millimetre'
        
    # supply a standard order
    if not kwargs.has_key('order'):
        kwargs['order'] = ['A','alpha','D','theta']
    # supply the standard denivit-hartenberg if no convention given
    if not kwargs.has_key('convention'):
        kwargs['convention'] = 'standard'

    row_length = 5
    nbr_of_joints = len(joint_values)
    if len(kwargs['table']) == 1 and type(kwargs['table'][0]) in [list, tuple]:
        raise ArithmeticError("Function does not use lists or tuples, please unpack using * operator.")
    elif not (len(kwargs['table']) % row_length == 0):
        raise ArithmeticError("Invalid number of Denavit-Hartenberg parameters - you also need to supply type of joint")

    matrices = []
    for k in xrange(0, nbr_of_joints):
        # Performing the operation
        # A, alpha, D, theta = params['A'], params['alpha'], params['D'], params['theta']
        # in a very general way
        var_names, var_values, joint_type = kwargs['order'],\
                                            kwargs['table'][row_length*k : row_length*k + row_length-1],\
                                            kwargs['table'][row_length*k + row_length-1]
        for i in xrange(row_length-1):
            exec('%s = %0.16f' % (var_names[i], var_values[i]))
            
        if kwargs['unit'] == 'mm' or kwargs['unit'] == 'millimetre':
            A = A * 1e-3 #converts mm to meters
            D = D * 1e-3 #convers mm to meters
        elif kwargs['unit'] == 'm' or kwargs['unit'] == 'metre':
            pass
        else:
            raise ArithmeticError("Unknown unit of length, only meters('m' or 'metre')"+\
            +" or millimeters ('mm' or 'millimetre') allowed.")
        if joint_type == 'R':
            matrices.append( transform_to_next(A, alpha, D, theta+joint_values[k]) )
        elif joint_type == 'P':
            matrices.append( transform_to_next(A, alpha, D+joint_values[k], theta) )

    # collect information about the Denivit-Hartenberg table
    dh_table = {
    'table' :    kwargs['table'],
    'unit':      kwargs['unit'],
    'convention':kwargs['convention'],
    'order':     kwargs['order'],
    }

    global_geometry = matmul_series(*matrices)
    global_geometry.insert(0, homogenous_matrix(0,0,0,0,0,0))
    
    result = {
        'T44': matmul(*matrices),
        'robot_geometry_local': matrices,
        'robot_geometry_global': global_geometry,
        'dh_table': dh_table
        }

    # perform matrix chain-multiplication / serial-multiplication with matrix product
    return result
#----------------------------------------------------------------------------------------------------------#
def calc_wcp(T44, L=None):
    return (T44[:,3] - T44[:,2]*L)[0:3]
#----------------------------------------------------------------------------------------------------------#
def inverse_kinematics_spherical_wrist(dh_table, j1, j2, j3, T44):
    #Calculate last angles
    R = T44[0:3,0:3]
    robot_info = forward_kinematics(j1, j2, j3, **dh_table)
    R3 = robot_info['T44'][:3,:3]

    R36 = R3.T.dot(R)
    X = R36[:,0]
    Y = R36[:,1]
    Z = R36[:,2]
    # for order of parameters check numpy.info(numpy.arctan2)
    j4 = atan2(Z[1],Z[0])
    j5 = atan2(norm(Z[0:2]), Z[2])
    j6 = atan2(X[2], Y[2]) + 90
    R36 = R36.T
    X = R36[:,0]
    Y = R36[:,1]
    Z = R36[:,2]
    # for order of parameters check numpy.info(numpy.arctan2)
    j41 = -(atan2(X[2], Y[2]) + 90)
    j51 = -atan2(norm(Z[0:2]), Z[2])
    j61 = -(atan2(Z[1],Z[0]))

    j42 = j4 + n.sign(j4-j41)*180
    j52 = j51
    j62 = j6 + n.sign(j6-j61)*180

    ##import pdb; pdb.set_trace()
    
    return j4, j5, j6, j41, j51, j61, j42, j52, j62
#----------------------------------------------------------------------------------------------------------#

