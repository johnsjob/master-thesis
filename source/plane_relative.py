from __future__ import division
#--------------------------#
import sys
#--------------------------#
import numpy as n
from numpy import cos, sin, pi
from numpy import array as mat
#=====================================================#
from helperfunctions import *
#sys.path.append("../int/misc-tools/")
#--------------------------#
num_points = 3
#--------------------------#
def get_plane_relative_point(plane, px, py, rot, tilt, skew, L):
    orig, basis_x, basis_y, normal = plane
    M = get_plane_relative_R(plane, rot, tilt, skew)
    return get_plane_point(plane, px, py) + L*M.dot(mat([0,0,1]))
#--------------------------#
def get_plane_relative_R(plane, rot, tilt, skew, flipped = True):
    orig, basis_x, basis_y, normal = plane
    transf_plane = mat([basis_x, basis_y, normal]).T

    rot_mat = lambda: rotation_matrix_rot_tilt_skew(rot, -tilt, skew)
    if not flipped:
        R = rot_mat()
    else:
        
        R = rot_mat().dot(mat_flip(1))
    M = matmul_series(transf_plane,R,transf_plane.T)
    return M
#--------------------------#
def mat_flip(M):
    transf_flip = mat([[1, 0, 0],
                       [0, -1, 0],
                       [0, 0, -1]])
    return transf_flip.dot(M)
#--------------------------#
if __name__ == '__main__':
    print; print "Init plots..."
    ax,_ = init_plot()
    #----------------------------------------#
    print; print "Define planes..."
    #untransformed_paper = define_plane([0,0,0],[1,0,0],[0,0,1])
    untransformed_paper = define_plane([0,0,0],[1,0,0],[0,1,0])
    (origin_paper, basis_x_paper, basis_y_paper, normal_paper) = untransformed_paper

    R = get_plane_relative_R(untransformed_paper,0,0,0,flipped=False)
    untransformed_paper = define_plane([0,0,0],R.dot([1,0,0]),R.dot([0,1,0]))
    (origin_paper, basis_x_paper, basis_y_paper, normal_paper) = untransformed_paper
    tmp = R

    #transformed_paper = define_plane([0,0,1],R.dot(basis_x_paper), R.dot(basis_y_paper))
    R = get_plane_relative_R(untransformed_paper,45,45,45)
    transformed_paper = define_plane(R.dot([0,0,-1]),R.dot(basis_x_paper), R.dot(basis_y_paper))
    (origin_transformed_paper, basis_x_transformed_paper,
     basis_y_transformed_paper, normal_transformed_paper) = transformed_paper
    #----------------------------------------#
    print; print "Placing points in coordinate systems ("+str(num_points)+")..."
    transformed_point = []
    untransformed_point = []
    rel_point = []
    for i in xrange(0, num_points):
        x, y = (rand(), rand())
        transformed_point.append( get_plane_point(transformed_paper,x,y) )
        untransformed_point.append( get_plane_point(untransformed_paper,x,y) )
    for i in xrange(0, 10):
        rot, tilt, skew = rand_range(-180,180), rand_range(-45,45), rand_range(-180,180)
        rel_point.append( get_plane_relative_point(untransformed_paper, untransformed_point[0][0], untransformed_point[0][1], rot, tilt, skew, 0.2) )
        rel_point.append( get_plane_relative_point(transformed_paper, untransformed_point[0][0], untransformed_point[0][1], rot, tilt, skew, 0.2) )

    transformed_point = array(transformed_point)
    untransformed_point = array(untransformed_point)
    rel_point = array(rel_point)

    transformed_directions = -diff(transformed_point.T).T
    untransformed_directions = -diff(untransformed_point.T).T
    untransformed_directions = untransformed_directions[:,0:2] #no z-coord
    ###========================================#
    ax.scatter(transformed_point[:,0], transformed_point[:,1], transformed_point[:,2])
    ax.scatter(untransformed_point[:,0], untransformed_point[:,1], untransformed_point[:,2])
    ax.scatter(rel_point[:,0], rel_point[:,1], rel_point[:,2],color='r')
    plot_plane(ax, transformed_paper)
    plot_plane(ax, untransformed_paper)
    show()

