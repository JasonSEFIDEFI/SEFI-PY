"""Physical, statistical, and information-boundary tests (stdlib unittest)."""
import unittest
import numpy as np
from bosonic import (annihilation, coherent_state, displacement, code_overlap,
                     erasure_kraus, choi, metrics)
from decoders import (GeometryDecoder, ConventionalPosterior, MaximumLikelihood,
                      AffineCoordinates, NonlinearCoordinates, covariance_from_calibration)
from statistics import interval


class PhysicalTests(unittest.TestCase):
    def test_coherent_poisson_moments_and_tail(self):
        psi, retained=coherent_state(1+.4j,40)
        n=np.arange(40)
        self.assertAlmostEqual(retained,1.,13)
        self.assertAlmostEqual(float(np.sum(n*np.abs(psi)**2)),1.16,12)
        self.assertAlmostEqual(float(np.vdot(psi,psi).real),1.,13)
    def test_ladder_and_displacement_vacuum(self):
        a=annihilation(24)
        self.assertAlmostEqual(a[0,1],1.)
        beta=.3+.2j
        expected,_=coherent_state(beta,24)
        np.testing.assert_allclose(displacement(beta,24)[:,0],expected,atol=2e-14)
    def test_unitary_and_inverse(self):
        d=displacement(.7-.4j,24)
        np.testing.assert_allclose(d.conj().T@d,np.eye(24),atol=1e-13)
        np.testing.assert_allclose(displacement(-.7+.4j,24),d.conj().T,atol=1e-13)
    def test_exact_oscillator_overlap_independent_matrix_reference(self):
        for beta in (0j,.2+.1j,-1.2+.7j,2.2-.4j):
            np.testing.assert_allclose(displacement(beta,48)[:2,:2],code_overlap(beta),atol=1e-12)
    def test_cptp_erasure_completeness_choi_partial_trace(self):
        for beta in (0j,.1+.4j,1.5j,2+.7j):
            ks=erasure_kraus(beta)
            np.testing.assert_allclose(sum(k.conj().T@k for k in ks),np.eye(2),atol=1e-12)
            c=choi(ks)
            self.assertGreaterEqual(np.linalg.eigvalsh(c).min(),-1e-12)
            self.assertAlmostEqual(float(np.trace(c).real),2.,12)
            # Matrix indexed [input,output,input,output] under column vectorization.
            partial=np.einsum('iaja->ij',c.reshape(2,3,2,3))
            np.testing.assert_allclose(partial,np.eye(2),atol=1e-12)
    def test_entangled_reference_fidelity(self):
        beta=.4-.25j; ks=erasure_kraus(beta)
        target=np.array([1,0,0,0,1,0],complex)/np.sqrt(2)
        bell=np.array([1,0,0,1],complex)/np.sqrt(2)
        result=sum(abs(np.vdot(target,np.kron(np.eye(2),k)@bell))**2 for k in ks)
        analytic=metrics(np.array([[beta.real,beta.imag]]))['entanglement_fidelity'][0]
        self.assertAlmostEqual(result,analytic,13)
    def test_six_state_fidelity(self):
        beta=.6+.3j; m=code_overlap(beta)
        states=[np.array([1,0]),np.array([0,1])]
        states += [np.array([1,s])/np.sqrt(2) for s in (1,-1,1j,-1j)]
        direct=np.mean([abs(np.vdot(s,m@s))**2 for s in states])
        self.assertAlmostEqual(direct,metrics(np.array([[.6,.3]]))['state_fidelity'][0],13)
    def test_no_postselection_and_metric_limits(self):
        scores=metrics(np.array([[0,0],[.3,.7],[3,2]]))
        self.assertEqual(scores['entanglement_fidelity'][0],1.)
        for key in ('entanglement_fidelity','state_fidelity','survival','accuracy'):
            self.assertTrue(np.all((scores[key]>=0)&(scores[key]<=1)))
        self.assertLess(scores['survival'][1],1.)
    def test_heterodyne_variance(self):
        rng=np.random.default_rng(910)
        pilots=rng.normal(0,np.sqrt(.5),(20000,16,2))
        np.testing.assert_allclose(np.var(pilots.mean(axis=1),axis=0),[1/32]*2,rtol=.03)
    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):annihilation(1)
        with self.assertRaises(ValueError):metrics(np.array([[np.nan,0]]))
        with self.assertRaises(ValueError):GeometryDecoder(-np.eye(2),np.eye(2))


class InferenceTests(unittest.TestCase):
    def setUp(self):
        self.s=np.array([[.05,.015],[.015,.02]])
        self.r=np.eye(2)/32
        self.y=np.random.default_rng(99).normal(size=(200,2))
    def test_ml_likelihood_stationarity(self):
        prediction=MaximumLikelihood().predict(self.y)
        np.testing.assert_equal(prediction,self.y)
    def test_quadratic_geometry_stationarity(self):
        z=GeometryDecoder(self.s,self.r).predict(self.y)
        residual=(z-self.y)@np.linalg.inv(self.r)+z@np.linalg.inv(self.s)
        np.testing.assert_allclose(residual,0,atol=1e-12)
    def test_conventional_posterior_equivalence(self):
        np.testing.assert_allclose(GeometryDecoder(self.s,self.r).predict(self.y),
                                   ConventionalPosterior(self.s,self.r).predict(self.y),atol=1e-14)
    def test_affine_metric_covariance(self):
        t=np.array([[1.7,.4],[-.3,.65]])
        np.testing.assert_allclose(AffineCoordinates(self.s,self.r,t).predict(self.y),
                                   GeometryDecoder(self.s,self.r).predict(self.y),atol=1e-14)
    def test_nonlinear_information_roundtrip(self):
        g=GeometryDecoder(self.s,self.r)
        np.testing.assert_allclose(NonlinearCoordinates(g).predict(self.y),g.predict(self.y),atol=1e-14)
    def test_zero_prior_and_large_prior_limits(self):
        np.testing.assert_equal(GeometryDecoder(np.zeros((2,2)),self.r).predict(self.y),0)
        np.testing.assert_allclose(GeometryDecoder(np.eye(2)*1e12,self.r).predict(self.y),self.y,atol=1e-12)
    def test_no_truth_argument_or_mutation(self):
        y=self.y.copy(); g=GeometryDecoder(self.s,self.r)
        first=g.predict(y); second=g.predict(y)
        np.testing.assert_equal(first,second);np.testing.assert_equal(y,self.y)
        self.assertEqual(set(g.__dict__),{'weight'})
    def test_calibration_reproducibility(self):
        rng=np.random.default_rng(23);x=rng.multivariate_normal([0,0],self.s,2000)
        np.testing.assert_allclose(covariance_from_calibration(x),np.cov(x,rowvar=False))
    def test_bootstrap_paired_identity_and_reproducibility(self):
        self.assertEqual(interval(np.zeros(20)),(0.,0.,0.))
        self.assertEqual(interval(np.arange(20)),interval(np.arange(20)))


if __name__=='__main__':unittest.main(verbosity=2)
