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
def matmul_symbolic(*matrix_factors):
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
def convert_to_mat(*args):
    """
        Converts a list of arguments to numpy.ndarray,
        the arguments are expected to be tuples or lists,
        scalar values are also accepted but it's not the functions
        intended use.
    """
    return map(mat, args)
#----------------------------------------#
