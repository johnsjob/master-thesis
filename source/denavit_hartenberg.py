from helperfunctions_math import *
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
def DH_params(*joint_values,**kwargs):
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
    'table' : kwargs['table'],
    'unit': kwargs['unit'],
    'convention':kwargs['convention'],
    'order': kwargs['order'],
    }    
    # perform matrix chain-multiplication / serial-multiplication with matrix product
    return matmul(*matrices), matrices, dh_table
#----------------------------------------------------------------------------------------------------------#
def calc_tool(_DH, *joint_values):
    DH = _DH.copy()
    joint_types = DH['table'][4 :: 5]
    revolute_index = DH['order'].index('theta')
    prismatic_index = DH['order'].index('D')
    DH['table'] = list( DH['table'] )
    
    for i,k in enumerate(joint_values):
        if joint_types[i] == 'R':
            DH['table'][i*5+revolute_index] += k
        elif k == 'P':
            DH['table'][i*5+prismatic_index] += k

    tool0, Ai, _  = DH_params(*joint_values, **DH)
    return tool0, Ai
#----------------------------------------------------------------------------------------------------------#
def calc_wcp(T44, L=None):
    return (T44[:,3] - T44[:,2]*L)[0:3]
