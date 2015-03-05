from __future__ import division
from numpy import dot

def matmul_series(matrix_factors):
	return reduce(dot, matrix_factors, 1)
