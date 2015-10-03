import numpy

from numpy import array as mat
from random import random

from time import time

s = [random() for x in xrange(100)]
s = mat(s).reshape(10,10)

start = time()
s = s.T
end = time()
print 'time: %0.32f' % (end-start)
