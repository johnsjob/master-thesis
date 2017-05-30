class TestIRB140(unittest.TestCase):

    def test_elbow_down(self):
        print 'test_elbow_down'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = 40
            j3 = -100
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['flange'], robot_info['robot_geometry_local']

            s = inverse_kinematics_elbow_down(DH_TABLE, A)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_up(DH_TABLE, A)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_fl, delbow_down_fl
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(1, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)

    def test_elbow_up(self):
        print '\ntest_elbow_up'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = 40
            j3 = -30
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['flange'], robot_info['robot_geometry_local']

            s = inverse_kinematics_elbow_up(DH_TABLE, A)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_down(DH_TABLE, A)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_fl, delbow_down_fl
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(0, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)

    def test_elbow_down_backward_facing(self):
        print '\ntest_elbow_down_backward'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = -90
            j3 = -30
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['flange'], robot_info['robot_geometry_local']

            s = inverse_kinematics_elbow_down(DH_TABLE, A, flipped = True)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_up(DH_TABLE, A, flipped = True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_fl, delbow_down_fl
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(3, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)

    def test_elbow_up_backward_facing(self):
        print '\ntest_elbow_up_backward'
        for _ in xrange(100):
            j1 =  rand_range(-180, 180)
            j2 = -40
            j3 = -100
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            robot_info = forward_kinematics(j1,j2,j3,j4,j5,j6,**DH_TABLE)
            A, debug = robot_info['flange'], robot_info['robot_geometry_local']

            sol = mat([j1,j2,j3,j4,j5,j6])
            s = inverse_kinematics_elbow_up(DH_TABLE, A, flipped = True)
            self.assertNotEqual(n.isnan(n.sum(s)), True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertAlmostEqual(b, j2)
            self.assertAlmostEqual(c, j3)

            s = inverse_kinematics_elbow_down(DH_TABLE, A, flipped = True)
            a,b,c = s[0][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)
            a,b,c = s[1][0:3]
            self.assertAlmostEqual(a, j1)
            self.assertNotAlmostEqual(b, j2)
            self.assertNotAlmostEqual(c, j3)

            # check if they are stored in the right order i.e.
            # elbow_up, elbow_down, elbow_up_backward, delbow_down_backward
            sols = inverse_kinematics_irb140(DH_TABLE, A)
            for i in xrange(2, len(sols.T), 4):
                a,b,c = sols[:,i][:3]
                self.assertAlmostEqual(a, j1)
                self.assertAlmostEqual(b, j2)
                self.assertAlmostEqual(c, j3)

    def test_non_reach_config(self):
            print '\ntest_non_reach_configs'
            for _ in xrange(0,100):
                j1 = rand_range(-180, 180)
                j2 = 90
                j3 = -89
                j4 = rand_range(-200, 200)
                j5 = rand_range(-115, 115)
                j6 = rand_range(-400, 400)

                s0 = j1,j2,j3,j4,j5,j6

                robot_info = forward_kinematics(*s0, **DH_TABLE)
                T44, debug1  = robot_info['flange'], robot_info['robot_geometry_local']
                sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
                sol = sol.T

                for i,s in enumerate(sol):
                    robot_info = forward_kinematics(*s, **DH_TABLE)
                    A, debug2  = robot_info['flange'], robot_info['robot_geometry_local']
                    if i in [l+m*8 for m,_ in enumerate(range(0, len(sol), 8)) for l in [0,1,4,5]]: #all non-flipped solutions only
                        self.assertAlmostEqual(norm(A-T44), 0)
                    else:
                        self.assertTrue(n.isnan(norm(A-T44)))

    def test_just_barely_reach_flipped_configs(self):
            print '\ntest_just_barely_reach_flipped_configs'
            for _ in xrange(0,100):
                j1 = rand_range(-180, 180)
                j2 = -90
                j3 = -89
                j4 = rand_range(-200, 200)
                j5 = rand_range(-115, 115)
                j6 = rand_range(-400, 400)

                s0 = j1,j2,j3,j4,j5,j6

                robot_info = forward_kinematics(*s0, **DH_TABLE)
                T44, debug1  = robot_info['flange'], robot_info['robot_geometry_local']
                sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
                for s in sol.T:
                    robot_info = forward_kinematics(*s, **DH_TABLE)
                    A, debug2  = robot_info['flange'], robot_info['robot_geometry_local']
                    self.assertAlmostEqual(norm(A-T44), 0)

    def test_forward_kinematics_general(self):
        print '\ntest_forward_kinematics_general'

        for counter in xrange(10000):
            fcounter = (counter / 10000.0)*100
            if fcounter % 1.0 == 0.0:
                print fcounter
            j1 = rand_range(-180, 180)
            j2 = rand_range(-90, 110)
            j3 = rand_range(-230, 50)
            j4 = rand_range(-200, 200)
            j5 = rand_range(-115, 115)
            j6 = rand_range(-400, 400)

            # makes sure we never end up at a singular point

            while (abs(j3) - 90) < 1e-7:
                j3 = rand_range(-230, 50)

            s0 = j1,j2,j3,j4,j5,j6
            robot_info = forward_kinematics(j1, j2, j3, j4, j5, j6, **DH_TABLE)
            T44, debug1 = robot_info['flange'], robot_info['robot_geometry_local']

            while norm(calc_wcp(T44,L=0.065)[:2]) < 1e-7:
                j2 = rand_range(-90, 110)
                T44, debug1  = forward_kinematics(j1, j2, j3, j4, j5, j6, **DH_TABLE)

            sol = mat( inverse_kinematics_irb140(DH_TABLE, T44) )
            num_valid_solutions = 0
            for s in sol.T:
                robot_info = forward_kinematics(*s, **DH_TABLE)
                A, debug2  = robot_info['flange'], robot_info['robot_geometry_local']
                num_valid_solutions += check_solution(*s)
                error = norm(A - T44)
                if not n.isnan(error):
                    self.assertAlmostEqual(error, 0)
            self.assertGreaterEqual(num_valid_solutions, 1)
            self.assertEqual(num_valid_solutions, calc_valid_raw_invkin_irb140(T44).shape[0])

            L = []
            for s in iterdim(sol,1):
                if check_solution(*s) == True:
                    L.append(s)
            L = mat(L)
            self.assertTrue(norm(calc_valid_raw_invkin_irb140(T44) - L) == 0.0)
#----------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':
    unittest.main()
