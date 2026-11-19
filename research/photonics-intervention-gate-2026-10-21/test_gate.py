import unittest
import numpy as np
from pilot import expected_fidelity,independent_quadrature,score

class GateChecks(unittest.TestCase):
    def test_exact_score_matches_independent_integral(self):
        cov=np.array([[.04,.01],[.01,.02]])
        for m in [[0,0],[.3,-.2],[1,.5]]:
            self.assertAlmostEqual(float(expected_fidelity(m,cov)[0]),independent_quadrature(m,cov),places=12)
    def test_zero_covariance_recovers_deterministic_score(self):
        x=np.random.default_rng(87).normal(size=(100,2))
        np.testing.assert_allclose(expected_fidelity(x,np.zeros((2,2))),score(x),atol=2e-16)
    def test_choice_is_model_optimal_within_two_actions(self):
        m=np.random.default_rng(2).normal(0,.2,(200,2));c=np.diag([.01,.004]);q=.01*np.eye(2)
        hold=expected_fidelity(m,c);act=expected_fidelity([0,0],c+q)[0]
        chosen=np.where(act>hold,act,hold)
        self.assertTrue(np.all(chosen>=hold));self.assertTrue(np.all(chosen>=act))
    def test_noise_can_make_correction_harmful(self):
        c=np.diag([.01,.004]);hold=expected_fidelity([0,0],c)[0];act=expected_fidelity([0,0],c+.04*np.eye(2))[0]
        self.assertGreater(hold,act)
if __name__=='__main__':unittest.main(verbosity=2)
