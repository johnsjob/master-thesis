from helperfunctions_plot import *
from plane_relative import *
from denavit_hartenberg140 import *

import itertools as it

def work_it(M, func=n.diff, axis=1):
    return np.apply_along_axis(func, axis, arr=M)

if __name__ == '__main__':
    ax, fig = init_plot(1.0)

    plane0 = define_plane_from_angles([0,0,0],
                                      30, 0, 0)
    plane_rel = define_plane_relative_from_angles(plane0, [0.1,0,0],
                                      30,45,0)


    ######
    plot_plane(ax, plane0, scale_factor=0.1)
    plot_plane(ax, plane_rel, scale_factor=0.1)

    #ax.scatter(point_matrix_tf[:,0],point_matrix_tf[:,1],point_matrix_tf[:,2])
    s = 0.3
    plot_equal_perspective(ax, [-s, s],
                               [-s, s],
                               [-s, s])
    show()

### note (orientations): local * global = global with relative local transformation
