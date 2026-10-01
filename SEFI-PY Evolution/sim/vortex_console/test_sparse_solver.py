"""Equivalence checks for optional sparse acceleration."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import unittest
import importlib.util
from unittest.mock import patch
from pathlib import Path
import numpy as np
import physics
import ring_operator as ring

@unittest.skipUnless(importlib.util.find_spec('scipy'), 'Install requirements-fast.txt for sparse equivalence tests')
class SparseSolverTests(unittest.TestCase):
    def test_stencil_and_boundary_source(self):
        for n,c in [(16,0),(16,.58),(24,1.2)]:
            dense,b,iu,iw=ring.setup(n,16/n,c)
            sparse,bs,ius,iws=ring.setup(n,16/n,c,sparse=True)
            np.testing.assert_allclose(sparse.toarray(),dense,rtol=0,atol=0)
            for a,z in [(b,bs),(iu,ius),(iw,iws)]:np.testing.assert_array_equal(a,z)

    def test_profiles_and_winding_constraints(self):
        reference,_=physics.load_profile(Path(__file__).parent/'data/charged_ring_refined.npz')
        for n,N in [(16,0),(32,0),(16,1)]:
            p=physics.validate(dict(n=n,N=N,c=.58))
            dense,dr=physics.solve(p,reference,backend='dense')
            sparse,sr=physics.solve(p,reference,backend='sparse')
            self.assertEqual(dr['converged'],sr['converged'])
            self.assertEqual(dr['has_ring'],sr['has_ring'])
            self.assertEqual(dr['reason'],sr['reason'])
            for key in ('u','w','s'):
                np.testing.assert_allclose(sparse[key],dense[key],rtol=1e-7,atol=1e-8)
            if N:self.assertTrue(np.all(sparse['s'][0,:]==0))
            if n==32 and not N:self.assertTrue(sr['converged'])
            self.assertLess(sr['operator_bytes'],dr['operator_bytes']/20)

    def test_missing_optional_dependency_and_cancellation(self):
        reference,_=physics.load_profile(Path(__file__).parent/'data/charged_ring_refined.npz')
        original=__import__
        def no_scipy(name,*args,**kwargs):
            if name.startswith('scipy'):raise ImportError('simulated missing scipy')
            return original(name,*args,**kwargs)
        with patch('builtins.__import__',side_effect=no_scipy):
            _,report=physics.solve(physics.validate(dict(n=16)),reference)
            self.assertEqual(report['solver_backend'],'dense')
            with self.assertRaises(ImportError):physics.solve(physics.validate({}),reference,backend='sparse')
        with self.assertRaises(InterruptedError):
            physics.solve(physics.validate({}),reference,cancelled=lambda:True,backend='sparse')

if __name__=='__main__':unittest.main()
