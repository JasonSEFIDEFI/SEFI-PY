"""Regression checks for the independent fixed-charge research mode."""
import json,tempfile,unittest
from pathlib import Path
import numpy as np
import rest_search,server

class FixedChargeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.old=server.RUNS;server.RUNS=Path(self.tmp.name)
    def tearDown(self):server.RUNS=self.old;self.tmp.cleanup()
    def test_energy_gradient_and_descent(self):
        s=rest_search.Search(L=16,h=1,Q=100,lam=12,N=1)
        self.assertLess(s.gradient_check(),1e-6);before=s.terms(s.y)['energy']
        for _ in range(30):s.descent_step()
        d=s.diagnostics();self.assertLess(d['energy'],before);self.assertAlmostEqual(d['nu']*d['I'],100)
    def test_guardrails(self):
        for p in [dict(Q=0),dict(Q=float('nan')),dict(N=1.5),dict(L=17,h=.6),dict(initial='../file'),dict(steps=30000,every=1)]:
            with self.assertRaises(ValueError):rest_search.validate(p)
    def test_checkpoint_and_deterministic_resume(self):
        p=rest_search.validate(dict(L=16,h=1,Q=100,steps=20,every=10))
        j=server.RestJob(p,'relax');j.run();self.assertEqual(j.state,'completed');self.assertGreater(j.residual,1e-5)
        self.assertIn('unconverged',j.message);self.assertTrue(j.checkpoint)
        frame=json.loads((j.dir/'frame_0000.json').read_text());self.assertEqual(frame['kind'],'relaxation');self.assertNotIn('time',frame['diagnostics']);self.assertEqual(frame['model_version'],rest_search.MODEL_VERSION)
        k=server.RestJob({**p,'steps':40},'restore_rest',restore=j.dir/'rest_checkpoint.npz');k.run()
        direct=rest_search.Search(L=16,h=1,Q=100,lam=12,N=1)
        for _ in range(40):direct.descent_step()
        self.assertEqual(k.solver_iteration,40);np.testing.assert_array_equal(k.search.y,direct.y)
        with np.load(k.dir/'rest_checkpoint.npz') as d:self.assertEqual(str(d['model_version']),rest_search.MODEL_VERSION)
    def test_saved_reference(self):
        with np.load(server.ROOT/'data/fixed_charge_Q1000_refined.npz') as d:
            p=json.loads(str(d['parameters']));s=rest_search.Search(**p);s.y=d['y']
        d=s.diagnostics();self.assertLess(d['residual'],1e-5);self.assertLess(d['nu'],1);self.assertEqual(d['core_radii'],[]);self.assertGreater(d['original_minimum'],.1);self.assertAlmostEqual(d['energy'],907.9235788928818,places=8)

if __name__=='__main__':unittest.main(verbosity=2)
