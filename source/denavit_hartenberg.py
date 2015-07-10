from helperfunctions_math import *
#----------------------------------------------------------------------------------------------------------#
def transform_to_next(A, alpha, D, theta):
    """
    Calculats transform from frame J = I-1 to I
    in order AI = Rzj( theta )Tzj( D )Txi( A )Rxi( alpha )
    and parameters are given in order A, alpha, D, theta.
    """
    Rz_J = homogenous_rotation_z(theta)
    Tz_J = homogenous_translation_z(D)
    Tx_I = homogenous_translation_x(A)
    Rx_I = homogenous_rotation_x(alpha)
    return matmul(Rz_J, Tz_J, Tx_I, Rx_I)
#----------------------------------------------------------------------------------------------------------#
def DH_params(*params, **kwargs):
    """
    Performs the denavit-hartenberg algorithm
    and calculates T44 = A0 * A1 * A2 * ... * An multiplication.

    The parameters are entered in the order A, alpha, D, theta repeatedly.
    """
    # handle which unit is being used
    if not kwargs.has_key('unit'):
        kwargs['unit'] = 'm'
    unit = kwargs['unit']

    nbr_of_sections = int(len(params) / 4)
    if len(params) == 1 and type(params[0]) in [list, tuple]:
        raise ArithmeticError("Function does not use lists or tuples, please unpack using *.")
    elif not (len(params) % 4 == 0):
        raise ArithmeticError("Invalid number of Denavit-Hartenberg parameters.")

    matrices = []
    for k in xrange(0, nbr_of_sections):
        A, alpha, D, theta = params[4*k:4*k+4]
        if unit == 'mm':
            A = A * 1e-3
            D = D * 1e-3
        elif unit == 'm':
            pass
        else:
            raise ArithmeticError("Unknown unit of length, only meters(\'m\') or millimeters meters(\'mm\') allowed.")
            
        matrices.append( transform_to_next(A, alpha, D, theta) )
    return matmul(*matrices), mat(matrices)
#----------------------------------------------------------------------------------------------------------#
def calc_wcp(T44, L=None):
    return (T44[:,3] - T44[:,2]*L)[0:3]
