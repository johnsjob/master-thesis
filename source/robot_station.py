from helperfunctions_plot import *
from plane_relative import *
from denavit_hartenberg140 import *

import itertools as it

def work_it(M, func=n.diff, axis=1):
    return np.apply_along_axis(func, axis, arr=M)

def get_closest_solutions_pair(s0, s1):
##    diff_list = []
##    index_list0 = []
##    index_list1 = []
##    for i0, k in enumerate(s0):
##        for i1, l in enumerate(s1):
##            diff_list.append(k-l)
##            index_list0.append(i0)
##            index_list1.append(i1)
##    index_list0 = mat(index_list0)
##    index_list1 = mat(index_list1)
##    diff_list = mat(diff_list)
##    norm_list = mat(map(norm, diff_list))
##    t = (norm_list - min(norm_list)) == 0.0
##    index0 = index_list0[t][0]
##    index1 = index_list1[t][0]
##    return mat((s0[index0], s1[index1]))
    data = []
    for i, s0i in enumerate(s0):
        for j, s1j in enumerate(s1):
            data.append([norm(s0i - s1j, ord = inf), i, j])
    data = mat(data)

    ret = []
    solution_col_row_pairs = n.argwhere(data == data.min(axis = 0)[0])
    solution_indices = solution_col_row_pairs[:,0]
    for solution_data in data[solution_indices]:
        norm_value, i, j = solution_data
        pair = mat([s0[i], s1[j]])
        return pair

def get_closest_solution(s0, s):
    diff_list = []
    index_list1 = []
    for i1, l in enumerate(s):
        diff_list.append(s0-l)
        index_list1.append(i1)
    index_list1 = mat(index_list1)
    diff_list = mat(diff_list)
    norm_list = mat(map(norm, diff_list))
    t = (norm_list - min(norm_list)) == 0.0
    index1 = index_list1[t][0]
    return s[index1]

def add_solutions(solutions, solution_value, index=5):
    for s in solutions.T:
        tmp1 = s.copy()
        tmp2 = s.copy()
        old_val = s[index]
        tmp1[index] = old_val + solution_value
        yield tmp1
        tmp2[index] = old_val - solution_value
        yield tmp2

def traverse_solutions(*args):
    for solutions in args:
        for s in solutions.T:
            yield s

def make_array(list_of):
    return mat(list_of).T

		
if __name__ == '__main__':
    for count in n.linspace(-180,180,10):
        ax, fig = init_plot()
        fig.clear()
        j1 =  180 #rand_range(-180, 180)
        j2 =  0#rand_range(-90, 110)
        j3 =  0#rand_range(-230, 50)
        j4 =  0#rand_range(-200, 200)
        j5 =  0#rand_range(-115, 115)
        j6 =  0#rand_range(-400, 400)
        joint_values = j1,j2,j3,j4,j5,j6

        T44, debug = forward_kinematics(*joint_values, **DH_TABLE)
        sol = inverse_kinematics_irb140(DH_TABLE, T44)

        plane0 = define_plane_from_angles([0,0,0],0, 0, 0)

        global_robot = matmul_series(*debug)
        global_robot.insert(0, debug[0])
        global_robot.insert(0, plane0)
        global_robot = mat(global_robot)
        global_robot_points = global_robot[:,:3,3]

        point_matrix = generate_symmetric_curve()
        point_matrix_tf = get_transformed_points(T44, point_matrix)

    ######
        ax = fig.add_subplot(1,2,1, projection='3d')
        for p in global_robot:
            plot_plane(ax, p, '--',scale_factor=0.1)

        ax.scatter(point_matrix_tf[:,0],point_matrix_tf[:,1],point_matrix_tf[:,2])
        ax.plot(global_robot_points[:,0], global_robot_points[:,1], global_robot_points[:,2],'k',linewidth=2)
        plot_equal_perspective(ax, [-0.5,0.5],[-0.5,0.5],[0,1])
        #show()
    ######
        plane = global_robot[-1]
        s = point_matrix_tf
        all_solutions = []
        for p in s:
            T44 = n.zeros((4,4))
            T44[:,3] = p
            T44[:3,:3] = plane[:3,:3]
            solutions = inverse_kinematics_irb140(DH_TABLE, T44)
            solutions = filter_solutions(solutions)
            print solutions.T.shape
            all_solutions.append(solutions.T)
        a = mat(all_solutions)

        import time
        start = time.time()
####        l = []
####        for i in xrange(len(a)-1):
####            l.append(get_closest_solutions_pair(a[i], a[i+1]))
####        l = mat(l)

        sol = []
        pair = get_closest_solutions_pair(a[0],a[1])
        sol.append(pair[0])
        for i in xrange(1,len(a)):
            sol.append(get_closest_solution(sol[i-1],a[i]))
        sol = mat(sol)
##        s = list(l[:,0,:])
##        s.append(l[-1,1,:])
##        s = mat(s)
        
        print 'stop: %0.2f' % (time.time() - start)
        r = work_it(work_it(sol, func=diff, axis=0),func=norm, axis=1)
        #r = n.max(n.abs(n.diff(sol,axis=0)),axis=1)
##        if (r >= 180.0).any():
##            print r
##            print n.round(n.max(n.abs(work_it(sol, func=diff, axis=0)),0))
##            import pdb; pdb.set_trace()

        ax0 = fig.add_subplot(1,2,2)
        ax0.plot(n.linspace(0,360,49),r);
        xlabel('curve angle')
        ylabel('solution distance')
        show()
        break
    print n.round(n.max(n.abs(work_it(sol, func=diff, axis=0)),0))
    #show()
    #plot(n.max(abs(s-sol), axis=1)); show()
