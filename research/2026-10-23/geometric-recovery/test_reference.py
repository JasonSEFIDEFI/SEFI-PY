import importlib.util
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('geo_reference',Path(__file__).with_name('run.py'))
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def test_projection_matches_nearest_circle():
    y=np.random.default_rng(99).normal(size=(1000,3))
    assert np.allclose(m.project(y),m.conventional(y),atol=2e-15)
    assert np.allclose(m.project(m.project(y)),m.project(y),atol=2e-15)

def test_projection_preserves_unknown_tangent_displacement():
    y=np.array([[1.,.2,.1]])
    assert np.linalg.norm(m.project(y)-[1,0,0])>.1

def test_convergence_orders_and_phase():
    for method,order in [('euler',1),('projected_euler',2),('midpoint',2),('rk4',4)]:
        a=m.oscillator(.01,method,duration=1)['state_error']
        b=m.oscillator(.005,method,duration=1)['state_error']
        assert abs(np.log2(a/b)-order)<.1
    r=m.oscillator(.1,'projected_euler',duration=10)
    assert r['energy_error']<1e-14 and r['phase_error']>.01
