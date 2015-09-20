#plot imports
from matplotlib.pylab import plot, hlines, xlim,\
xlabel, ylabel, plt, grid, title, legend, show
from numpy import array as mat
import numpy as np

def init_plot(aspect=1/4.0):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=plt.figaspect(aspect))
    ax = fig.gca(projection='3d') 
    return ax, fig
#----------------------------------------#
def plot_plane(ax, plane,style='-',scale_factor=1):
    s =mat([plane[:3,3], plane[:3,3] + plane[:3,0]*scale_factor])
    ax.plot(s[:,0],s[:,1],s[:,2],'b'+style)

    s =mat([plane[:3,3], plane[:3,3] + plane[:3,1]*scale_factor])
    ax.plot(s[:,0],s[:,1],s[:,2],'g'+style)

    s =mat([plane[:3,3], plane[:3,3] + plane[:3,2]*scale_factor])
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
    ax, _ = init_plot()
    
