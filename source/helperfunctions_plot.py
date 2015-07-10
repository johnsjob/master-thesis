#plot imports
from matplotlib.pylab import plot, hlines, xlim,\
xlabel, ylabel, plt, grid, title, legend, show

def init_plot():
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.gca(projection='3d') 
    return ax, fig
#----------------------------------------#
def plot_plane(ax, plane):
    dx = array([get_plane_point(plane,0,0), get_plane_point(plane,1,0)])
    dy = array([get_plane_point(plane,0,0), get_plane_point(plane,0,1)])
    dz = array([get_plane_point(plane,0,0), get_plane_point(plane,0,0) + plane[3]])
    ax.plot(dx[:,0],dx[:,1],dx[:,2],'b')
    ax.plot(dy[:,0],dy[:,1],dy[:,2],'g')
    ax.plot(dz[:,0],dz[:,1],dz[:,2],'r')
#----------------------------------------#
def plot_equal_perspective(ax, xList, yList, zList):
    xmin = xList[0]; xmax = xList[1]
    ymin = yList[0]; ymax = yList[1]
    zmin = zList[0]; zmax = zList[1]
    Xb = 0.5*abs(xmax-xmin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(xmax + xmin)
    Yb = 0.5*abs(ymax-ymin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(ymax + ymin)
    Zb = 0.5*abs(zmax-zmin)*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(zmax + zmin)
    # Comment or uncomment following both lines to test the fake bounding box:
##    for xb, yb, zb in zip(Xb, Yb, Zb):
##       ax.plot([xb], [yb], [zb], 'w')    
#----------------------------------------#
