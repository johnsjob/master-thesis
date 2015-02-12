import parsingtools as parse

import numpy as n
import matplotlib.pylab as pl

from numpy import array as mat
import time

from functools import reduce # Valid in Python 2.6+, required in Python 3
import operator

def cumprod(L):
        return reduce(operator.mul, L, 1)

def do_sort_param_x(data, param_x, param_y):
        return do_sort_x(data[param_x], data[param_y])

def do_sort_x(px, py):
        out = zip(px, py)
        out_s = sorted(out, key=lambda x: x[0], reverse=False)
        return zip(*out_s)

def do_plot_params(data, param_x, param_y, style = '+', title='Parameter plot', tol=None):
        x,y = do_sort_param_x(data, param_x, param_y)
        pl.title(title)
        pl.xlabel(param_x)
        pl.ylabel(param_y)
        pl.plot(x, y, style)

def do_plot(data, param_x, param_y, style = '+', title='Standard plot', tol=None):
        x,y = do_sort_param_x(data, param_x, param_y, label=param_y)
        pl.title(title)
        pl.xlabel(param_x)
        pl.plot(x, y, style)

def do_show():
        pl.legend()
        pl.show()

def do_parse(filepath):
        # parse the file
        start_time = time.time()
        data = parse.parse_file(filepath)
        stop_time = time.time()
        print "\ntime: "+str(stop_time - start_time)

        print n.sort(data.keys())
        return data

def do_err_floor(L, tol=2e-3):
        tol = tol
        return [n.floor(x/tol)*tol for x in L]

def get_submat_maxerr(E, lab='E', tol=2e-3, condition='largerthan'):
        """
HELP:
======================================================================
Function to retrieve the list of elements in matrix form where
the max-error is larger than a given tolerance.

input: error matrix, e.g abs(R - Rt)
======================================================================
        """
        if type(E) in [list, tuple]:
                help(get_submat_maxerr)
                raise Exception("Wrong type, also make sure it's an error matrix - see help above.")
        sh = E.shape
        if len(sh) <= 1:
                raise Exception("Wrong shape of matrix.");
        elif len(sh) == 2:
                  if not (n.mod( n.sum(sh), 2) == 0):
                          raise Exception("Matrix is not square.");
        elif len(sh) == 3:
                  if not (n.mod( n.sum(sh[1:]), 2) == 0):
                          raise Exception("List of matrices are not square.");
        else:
                raise Exception("This typbe of matrix / list is not supported yet.");
        if len(sh) == 2:
                size = cumprod(sh)
                rowlen = sh[0]
        elif len(sh) == 3:
                size = cumprod(sh[1:])
                rowlen = sh[1]
                listlen = sh[0]
        #get indexes
        if condition == 'largerthan':
                I = n.max(E,0) > tol
        elif condition == 'smallerthan':
                I = n.max(E,0) < tol
        elif condition == 'smallerthaninclusive':
                I = n.max(E,0) <= tol
        else:
                raise Exception("Unsupported condition."+str("\n%79s" % "Implemented are: 'largerthan', 'smallerthan', 'smallerthaninclusive'"))

        #get submatrix
        if len(sh) == 3:
                subE = E[:,I]
        elif len(sh) == 2:
                subE = E[I]

        #generate labels
        labels = [lab + str(i / rowlen + 1) + str( n.mod(i, rowlen) + 1) for i in xrange(0, size)]
        labels = mat(labels).reshape(rowlen, rowlen)
        labels = labels[I]
        return subE, labels

def do_plot_submat(E, labels, dat=None, param_x='Index', param_y='error', style='-', tol=2e-3):
        if dat is not None:
                x,y = do_sort_x(data[param_x],E)
                pl.plot(x, y, style);
                pl.xlabel(param_x)
        else:
                pl.plot(E, style);
                pl.xlabel('Index')
        pl.ylabel(param_y)

        pl.axhline(tol, color='k', linewidth=2)
        pl.legend(labels)
        #pl.ylim(0, 13 * tol)
        pl.show()

import os, sys
if __name__ == '__main__':
        var = raw_input("DID YOU REMEMBER TO MAKE? (Y/N): ")

        if 'n' in var.lower():
                sys.exit(0)

        # path to the outputed data from the test-case
        path = '/home/johnnys/git/orientationmath/source/'
        fil = 'test.txt'
        filepath = path + fil

        data = do_parse(filepath)
        err = do_err_floor(data['#err'])
        data['max_err'] = do_err_floor(data['max_err'])
        
        err, labels = get_submat_maxerr( n.reshape(err, (20000, 3, 3)), 'e')
        do_plot_submat(err, labels, data, param_x='skew', param_y='max-error')

        #check xro against yro - overflow
        #check t against _s - overflow
