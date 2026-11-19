import unittest
import numpy as np
from adaptive import R,project,posterior,predict,EW
from benchmark import generate,fidelity

class ScientificChecks(unittest.TestCase):
    def test_covariance_projection_matches_independent_eigendecomposition(self):
        rng=np.random.default_rng(41)
        for _ in range(100):
            x=rng.normal(size=(2,2));x=(x+x.T)/2
            vals,vecs=np.linalg.eigh(x);expected=(vecs*np.maximum(vals,1e-8))@vecs.T
            a,b,c=project(x[0,0],x[0,1],x[1,1])
            np.testing.assert_allclose([[a,b],[b,c]],expected,atol=2e-15)
    def test_posterior_matches_matrix_formula(self):
        rng=np.random.default_rng(19)
        for _ in range(100):
            x=rng.normal(size=(2,2));s=x@x.T;mu=rng.normal(size=2);y=rng.normal(size=2)
            out,ll=posterior(y,mu,(s[0,0],s[0,1],s[1,1]))
            expected=mu+s@np.linalg.solve(s+R*np.eye(2),y-mu)
            np.testing.assert_allclose(out,expected,atol=1e-13)
            self.assertTrue(np.isfinite(ll))
    def test_causal_future_invariance(self):
        cal,y,_=generate(77,'joint_drift');changed=y.copy();changed[60:]+=100
        for kind in ['fixed','ew:0.03','rolling:64','mixture']:
            a=predict(kind,cal,y[:100]);b=predict(kind,cal,changed[:100]);short=predict(kind,cal,y[:60])
            np.testing.assert_allclose(a[:60],b[:60],atol=0)
            np.testing.assert_allclose(short,a[:60],atol=0)
    def test_score_matches_legacy_physical_overlap(self):
        import physical_reference as m
        r=np.random.default_rng(12).normal(size=(200,2));matrix=m.code_overlap(r[:,0]+1j*r[:,1])
        expected=np.abs(np.trace(matrix,axis1=1,axis2=2))**2/4
        np.testing.assert_allclose(fidelity(r,np.zeros_like(r)),expected,atol=3e-16)
        self.assertEqual(float(fidelity(np.zeros((1,2)),np.zeros((1,2)))[0]),1.)
    def test_ew_predicts_before_update(self):
        cal,y,_=generate(55,'stationary');e=EW(cal,.1)
        expected=posterior(y[0],e.m.mean,e.m.covariance())[0]
        np.testing.assert_allclose(e.step(y[0]),expected)
    def test_no_truth_interface(self):
        import inspect
        self.assertEqual(list(inspect.signature(predict).parameters),['kind','calibration','observations'])
    def test_development_and_evaluation_are_disjoint(self):
        from benchmark import P
        self.assertFalse(set(P['development_seeds'])&set(P['evaluation_seeds']))
if __name__=='__main__':unittest.main(verbosity=2)
