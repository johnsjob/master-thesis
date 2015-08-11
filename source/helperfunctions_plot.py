#plot imports
from matplotlib.pylab import plot, hlines, xlim,\
xlabel, ylabel, plt, grid, title, legend, show

from plane_relative import *
from numpy import array
import numpy as np

def init_plot():
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca(projection='3d') 
    return ax, fig
#----------------------------------------#
def plot_plane(ax, plane,style='-'):
    s =mat([plane[:3,3], plane[:3,3] + plane[:3,0]])
    ax.plot(s[:,0],s[:,1],s[:,2],'b'+style)

    s =mat([plane[:3,3], plane[:3,3] + plane[:3,1]])
    ax.plot(s[:,0],s[:,1],s[:,2],'g'+style)

    s =mat([plane[:3,3], plane[:3,3] + plane[:3,2]])
    ax.plot(s[:,0],s[:,1],s[:,2],'r'+style)
#----------------------------------------#
def plot_equal_perspective(ax, xList, yList, zList):
    xmin = xList[0]; xmax = xList[1]
    ymin = yList[0]; ymax = yList[1]
    zmin = zList[0]; zmax = zList[1]
    Xb = 0.5*abs(xmax-xmin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(xmax + xmin)
    Yb = 0.5*abs(ymax-ymin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(ymax + ymin)
    Zb = 0.5*abs(zmax-zmin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(zmax + zmin)


    ax.set_xlim3d(xmin, xmax)
    ax.set_ylim3d(ymin, ymax)
    ax.set_zlim3d(zmin, zmax)
##    for xb, yb, zb in zip(Xb, Yb, Zb):
##       ax.plot([xb], [yb], [zb], 'w')    
#----------------------------------------#
if __name__ == '__main__':
    print; print "Init plots..."
    ax,_ = init_plot()
    #----------------------------------------#
    ### flipped = True means that we are plotting the coordinate system locally,
    ### X,Y,Z axes align with Euclidean axes
    ###
    ### flipped = False means that we are plotting the coordinate system
    ### globally in the robot world system
    flipped = True
    r, t, s = 0,0,0
    print; print "Define planes..."
    plane1 = define_plane_from_angles([1,0,0],45,45,180)
    plane2 = define_plane_from_angles([0,0,0],0,10,0)
    plane3 = define_plane_relative_from_plane(plane1, plane2)
    plane4 = define_plane_relative_from_angles(plane3, [0,0,0], 0,10,0)

    func = lambda x: n.sin(x)
    s = generate_curve(ampl_factor=1, y_func = func, freq=1, offset=-0.0)
    ### <[v.T 1], [[R.T t],[0 1]]>
    q = get_transformed_points(plane1, s)
    ax.scatter(q[:,0], q[:,1], q[:,2])
    
    plot_equal_perspective(ax, [-2,2], [-2,2], [-2,2])
    plot_plane(ax, plane1)
    plot_plane(ax, plane3, '-.')
    plot_plane(ax, plane4, '--')

    show()
