from pylab import plot, show, grid, title, clf
import cPickle as pickle
from numpy import mean, array as mat
with open(r"C:\Users\***REMOVED***\Dropbox\exjobb\results\measurements\measurements_automatic\boxerpen_upright_manual.pickle") as fp:
   p1 = pickle.load(fp)
with open(r"C:\Users\***REMOVED***\Dropbox\exjobb\results\measurements\measurements_automatic\boxerpen_20deg_manual.pickle") as fp:
   p2 = pickle.load(fp)
p = mat( list(p1) + list(p2) )
mx = mean(p[:,0])
my = mean(p[:,1])
print len(p)

title("Boxer pen")
plot(p[:,0]-mx, p[:,1]-my, '.b');
grid(); show()

# - - -

with open(r"C:\Users\***REMOVED***\Dropbox\exjobb\results\measurements\measurements_manual\pen_infos\ED117_pen_tip.pickle") as fp:
    p = pickle.load(fp)

xy = mat(p['coords'])
mxy = mean(xy, axis=0)
xy = (xy - mxy) * 0.3
print len(xy)

title("ED pen")
plot(xy[:,0], xy[:,1], '.b');
grid(); show()
