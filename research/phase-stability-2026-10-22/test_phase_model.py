import unittest,math
from phase_model import rhs,jacobian,equilibrium,critical_detuning,simulate

class PhaseChecks(unittest.TestCase):
    def test_threshold_satisfies_cubic_and_decreases_with_activity_penalty(self):
        last=3.
        for alpha in [0,.25,1,4]:
            q=critical_detuning(2,alpha,1)
            self.assertAlmostEqual(q*(1+alpha*q*q),2,places=12);self.assertLess(q,last);last=q
    def test_equilibrium_and_stability(self):
        for delta in [-.95,-.4,0,.4,.95]:
            p,d=equilibrium(delta,2,1,1);a,b=rhs(p,d,delta,2,1,1)
            self.assertAlmostEqual(a,0,places=12);self.assertAlmostEqual(b,0,places=12)
            j=jacobian(p,d,delta,2,1,1);trace=j[0][0]+j[1][1];det=j[0][0]*j[1][1]-j[0][1]*j[1][0]
            self.assertLess(trace,0);self.assertGreater(det,0);self.assertAlmostEqual(det,2*math.cos(p)/(1+d),places=12)
    def test_jacobian_independent_finite_differences(self):
        p,d=.6,.4;j=jacobian(p,d,.8,2,1,1);h=1e-6
        for col in [0,1]:
            xp=[p,d];xm=[p,d];xp[col]+=h;xm[col]-=h
            fp=rhs(*xp,.8,2,1,1);fm=rhs(*xm,.8,2,1,1)
            for row in [0,1]:self.assertAlmostEqual(j[row][col],(fp[row]-fm[row])/(2*h),places=8)
    def test_constant_gain_limit_matches_exact_phase(self):
        r,_=simulate(.8,2,0,1,.01,60,10)
        self.assertAlmostEqual(r['final_phase'],math.asin(.4),places=10);self.assertEqual(r['final_degradation'],0.)
    def test_strict_threshold_excludes_marginal_point(self):
        self.assertIsNone(equilibrium(1,2,1,1));self.assertIsNotNone(equilibrium(.999,2,1,1));self.assertIsNone(equilibrium(1.001,2,1,1))
if __name__=='__main__':unittest.main(verbosity=2)
