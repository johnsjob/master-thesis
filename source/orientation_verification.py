from helperfunctions_math import rotation_matrices, \
homogenous_matrix as hom, rand_range as rnd, pi,\
rotation_matrix_rot_tilt_skew as rts_mat, mat, matmul,\
rot_tilt_skew, nmap

from numpy import set_printoptions
set_printoptions(precision=4)

from numpy import sign, diag, inf, mean
from numpy.linalg import solve, svd, det, norm

from pylab import plot, show

def fit_to_ori(noisy_ori):

    U,s,V = svd(noisy_ori)
    S = diag(s)

    R1 = U.dot(V)
    R2 = U.dot( diag((1,1,sign(det(U.dot(V))))) ).dot(V)

    return R2, (U,S,V)

def solve_ori(rots1, rots2):
    a,b,c = rots1.shape
    d,e,f = rots2.shape

    Y = rots1.reshape(a*b, c)
    Z = rots2.reshape(d*e, f)

    B = Y.T.dot(Y)
    C = Y.T.dot(Z)
    A = solve(B,C)
    return fit_to_ori(A)

# ranges:
# rot = (-pi, pi)
# tilt = (-pi/2, pi/2)
# skew = (-pi, pi)

e = lambda mag: mat([rnd(-1, 1)*(10**mag) for _ in xrange(9)]).reshape(3,3)

def ang(R):
    return mat(rot_tilt_skew(R))

def apply_noise(R, mag=-4):
    return R + e(mag)

if __name__ == "__main__":
    res = {"solerr1":[],
           "solerr2":[],
           "solerr32":[],
           "angerr1":[],
           "angerr2":[],
           "angerr32":[],
           }
    for i in xrange(1000):
        if i % 100 == 0:
            print(i)

        rts = ((rnd(-pi, pi),
            rnd(-pi/2.0, pi/2.0),        rnd(-pi, pi)) for _ in range(100))

        hom_oris = rotation_matrices(rts)

        x = ref = rts_mat(10,20,30)
        y = oris = hom_oris[:,:3,:3]
        z = nmap(ref.dot, oris)

        y = nmap(apply_noise, y)
        z = nmap(apply_noise, z)

        Y = y.reshape(300,3)
        Z = z.reshape(300,3)

        B = Y.T.dot(Y)
        C = Y.T.dot(Z)

        A = solve(B,C)

        U,s,V = svd(C)
        S = diag(s)

        R1 = U.dot(V)
        R2 = U.dot( diag((1,1,sign(det(U.dot(V)))))).dot(V)
        R32, _ = fit_to_ori(A)

        res["solerr1"].append(abs(R1 - x).max())
        res["solerr2"].append(abs(R2 - x).max())
        res["solerr32"].append(abs(R32 - x).max())
        res["angerr1"].append(norm(ang(R1) - ang(x), inf))
        res["angerr2"].append(norm(ang(R2) - ang(x), inf))
        res["angerr32"].append(norm(ang(R32) - ang(x), inf))

    print "MAX:"
    print max(res['solerr1'])
    print max(res['solerr2'])
    print max(res['solerr32'])
    print max(res['angerr1'])
    print max(res['angerr2'])
    print max(res['angerr32'])
    print "MEAN:"
    print mean(res['solerr1'])
    print mean(res['solerr2'])
    print mean(res['solerr32'])
    print mean(res['angerr1'])
    print mean(res['angerr2'])
    print mean(res['angerr32'])
    print "MIN:"
    print min(res['solerr1'])
    print min(res['solerr2'])
    print min(res['solerr32'])
    print min(res['angerr1'])
    print min(res['angerr2'])
    print min(res['angerr32'])
