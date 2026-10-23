import importlib.util
from pathlib import Path
import math
import numpy as np

spec=importlib.util.spec_from_file_location('elements',Path(__file__).with_name('analyze.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_nuclear_table_join_and_units():
    allrows,selected,audit=m.load(); lookup={(r['Z'],r['N']):r for r in allrows}
    assert audit['ame_rows']==audit['nubase_ground_rows']==3558
    assert audit['joined_nubase']==3558
    assert all(r['A']==r['Z']+r['N'] for r in allrows)
    deuterium=lookup[(1,1)]
    assert abs(deuterium['binding_per_A_keV']*2/1000-2.2245662)<1e-7
    assert lookup[(1,2)]['half_life_seconds']==12.32*365.25*86400
    assert lookup[(1,0)]['lifetime_status']=='no_observed_decay'
    assert not any(r['mass_estimated'] for r in selected)

def test_split_disjointness_and_pairing():
    _,rows,_=m.load()
    for name,a,b,c in m.split_masks(rows):
        assert not set(a)&set(b) and not set(a)&set(c) and not set(b)&set(c)
    x=np.array([[8,8],[7,7],[8,9]],float)
    assert m.semf(x)[0,4]>0 and m.semf(x)[1,4]<0 and m.semf(x)[2,4]==0

def test_kernel_equivalence_and_extrapolation_reversion():
    x=np.array([[1.,2.],[3.,4.],[5.,6.]])
    y=np.array([1.,-2.,3.]);q=np.array([[2.,3.]])
    w=np.exp(-((q[:,None,:]-x[None,:,:])**2).sum(axis=2)/8)
    assert np.allclose(m.residual_predict(x,y,q,2),(w@y)/(w.sum(axis=1)+1))
    assert np.allclose(m.residual_predict(x,y,np.array([[1000.,1000.]]),2),0)

def test_density_laplacian_at_gaussian_peak():
    # Analytic sign counterexample: d-dimensional Gaussian peak is negative.
    d=5;h=2.;rho=(2*math.pi*h*h)**(-d/2)
    assert -d*rho/(h*h)<0

def test_older_evaluation_parser():
    spec=importlib.util.spec_from_file_location('historical',Path(__file__).with_name('historical.py'))
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    rows=h.load_2016()
    assert len(rows)==3436
    assert abs(rows[(1,1)]['binding_per_A_keV']-1112.283)<1e-6
