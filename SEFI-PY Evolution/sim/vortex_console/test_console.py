"""Reproducible numerical and API tests; uses a temporary run directory."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
os.environ['OMP_NUM_THREADS']='2'
import json,tempfile,time,threading,unittest,urllib.request,urllib.error
from unittest.mock import patch
from pathlib import Path
import numpy as np
import physics
import server

ROOT=Path(__file__).resolve().parent

class NumericalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.f,cls.p=physics.load_profile(ROOT/'data/charged_ring_refined.npz')
    def test_validation(self):
        for change in [dict(c=float('nan')),dict(n=20.5),dict(N=-1),dict(dt=.1),dict(duration=100,save_every=.025)]:
            with self.assertRaises(ValueError):physics.validate(change)
        for change in [dict(lam=12),dict(N=1)]:
            with self.assertRaises(ValueError):physics.validate(change,True)
    def test_stationary_control_and_perturbation(self):
        control=physics.Evolution(self.f,{**self.p,'epsilon':0})
        self.assertLess(control.initial_rhs,1e-9)
        for _ in range(40):control.step()
        d=control.diagnostics();self.assertLess(d['phi_distance'],1e-8);self.assertLess(d['charge_drift'],1e-10)
        a=physics.Evolution(self.f,self.p);b=physics.Evolution(self.f,self.p)
        for _ in range(40):a.step(.025)
        for _ in range(80):b.step(.0125)
        delta=float(abs(a.y-b.y).max());self.assertLess(delta,1e-7);self.assertGreater(a.diagnostics()['phi_distance'],1e-4)
        print('NUMERICS',json.dumps(dict(initial_rhs=control.initial_rhs,control=d,dt_max_difference=delta,perturbed=a.diagnostics())),flush=True)
    def test_changed_parameter_solve(self):
        p={**self.p,'c':.58};f,result=physics.solve(p,self.f)
        self.assertTrue(result['converged']);self.assertTrue(result['has_ring']);self.assertGreater(float(abs(f['u']-self.f['u']).max()),1e-4)
        print('CHANGED_PROFILE',json.dumps(result),flush=True)

class APITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();server.RUNS=Path(cls.tmp.name);server.JOBS.clear()
        cls.http=server.ThreadingHTTPServer(('127.0.0.1',0),server.Handler);server.PORT=cls.http.server_port
        cls.url=f'http://127.0.0.1:{server.PORT}'
        threading.Thread(target=cls.http.serve_forever,daemon=True).start()
    @classmethod
    def tearDownClass(cls):cls.http.shutdown();cls.http.server_close();cls.tmp.cleanup()
    def request(self,path,body=None,header=True):
        headers={'Content-Type':'application/json'}
        if header:headers['X-Vortex-Console']='1'
        req=urllib.request.Request(self.url+path,data=json.dumps(body).encode() if body is not None else None,headers=headers)
        with urllib.request.urlopen(req,timeout=10) as r:return json.load(r)
    def wait(self,id,states=('completed','failed','stopped'),limit=20):
        end=time.monotonic()+limit
        while time.monotonic()<end:
            s=self.request('/api/status?id='+id)
            if s['state'] in states:return s
            time.sleep(.03)
        self.fail('Job timed out')
    def test_api_and_checkpoint(self):
        b=self.request('/api/baseline');self.assertAlmostEqual(b['diagnostics']['charge'],7.1333804,places=6)
        for body in [dict(parameters={'lam':12}),dict(action='sweep',values=[.1],field='bad')]:
            with self.assertRaises(urllib.error.HTTPError):self.request('/api/start',body)
        with self.assertRaises(urllib.error.HTTPError):self.request('/api/start',{},False)
        with self.assertRaises(urllib.error.HTTPError):self.request('/api/checkpoint?id=../server.py')
        p={**physics.DEFAULTS,'duration':.5,'save_every':.25}
        j=self.request('/api/start',dict(action='evolve',parameters=p));s=self.wait(j['id']);self.assertEqual(s['state'],'completed');self.assertEqual(s['frames'],3)
        f=self.request('/api/frame?id='+j['id']+'&index=2');self.assertAlmostEqual(f['diagnostics']['time'],.5)
        k=self.request('/api/start',dict(action='restore',checkpoint=j['id'],parameters={**p,'duration':1}));s=self.wait(k['id']);self.assertEqual(s['state'],'completed')
        direct=physics.Evolution(server.REFERENCE,p)
        for _ in range(40):direct.step(.025)
        restored=server.JOBS[k['id']].ev
        self.assertLess(float(abs(restored.y-direct.y).max()),1e-13);self.assertEqual(restored.q0,direct.q0)
        print('CHECKPOINT_MAX_DIFFERENCE',float(abs(restored.y-direct.y).max()),flush=True)
        sweep=self.request('/api/start',dict(action='sweep',parameters=p,field='epsilon',values=[0,.002]));s=self.wait(sweep['id']);self.assertEqual(s['state'],'completed');self.assertEqual(len(s['members']),2)
        self.assertLess(s['members'][0]['final']['phi_distance'],1e-8);self.assertGreater(s['members'][1]['final']['phi_distance'],1e-4)
        # Long run lets lifecycle commands reach step boundaries before completion.
        j=self.request('/api/start',dict(action='evolve',parameters={**p,'duration':100,'save_every':1}))
        self.wait(j['id'],('running',));self.request('/api/control',dict(id=j['id'],command='pause'));self.wait(j['id'],('paused',))
        self.request('/api/control',dict(id=j['id'],command='checkpoint'))
        self.request('/api/control',dict(id=j['id'],command='resume'));self.wait(j['id'],('running',))
        self.request('/api/control',dict(id=j['id'],command='stop'));s=self.wait(j['id']);self.assertEqual(s['state'],'stopped');self.assertTrue(s['checkpoint'])

    def test_failed_profile_blocks_evolution(self):
        p=physics.validate({'c':.59,'duration':.5})
        result=dict(converged=False,has_ring=True,carrier_peak=.3,residual=.1,reason='test-injected nonconvergence')
        with patch.object(physics,'solve',return_value=(server.REFERENCE,result)):
            j=server.Job(p,'evolve');j.run()
        self.assertEqual(j.state,'failed');self.assertIsNone(j.ev);self.assertEqual(len(j.frames),1)

    def test_checkpoint_version_guard(self):
        ident='f'*32;folder=server.RUNS/ident;folder.mkdir(exist_ok=True)
        np.savez(folder/'checkpoint.npz',model_version='incompatible')
        with self.assertRaises(urllib.error.HTTPError):self.request('/api/start',dict(action='restore',checkpoint=ident,parameters={}))

if __name__=='__main__':unittest.main(verbosity=2)
