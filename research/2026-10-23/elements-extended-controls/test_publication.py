import importlib.util
from pathlib import Path
import numpy as np
spec=importlib.util.spec_from_file_location('publication',Path(__file__).with_name('publication_benchmark.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_gp_mean_against_direct_solve_and_positive_variance():
    x=np.array([[10.,10.],[12.,12.],[14.,14.]])
    y=np.array([1.,-2.,3.]);q=np.array([[11.,11.]])
    g=m.GP(x,y,2.,.1)
    mean,sd=g.predict(q,True)
    k=g.signal**2*m.kernel(q,x,2.,'RBF')
    K=g.signal**2*m.kernel(x,x,2.,'RBF')+np.eye(3)*(g.nugget**2+1e-10)
    assert np.allclose(mean,k@np.linalg.solve(K,y),atol=1e-12)
    expected=g.signal**2-float((k@np.linalg.solve(K,k.T))[0,0])+g.nugget**2
    assert np.allclose(sd**2,expected,atol=1e-12) and np.all(sd>0)

def test_kernel_covariance_and_quantile_rank():
    x=np.random.default_rng(3).normal(size=(20,2))
    for family in ['RBF','Matern32']:
        K=m.kernel(x,x,1.5,family)
        assert np.allclose(K,K.T)
        assert np.linalg.eigvalsh(K).min()>-1e-12
    assert m.interval_quantile(np.arange(9),.1)==8
    assert np.isinf(m.interval_quantile([1,2,3],.1))

def test_partition_no_leakage():
    for _,x,_,tu,_,ca,_,te,_ in m.partitions():
        ids=lambda q:{tuple(z) for z in q}
        groups=[ids(z) for z in [x,tu,ca,te]]
        assert all(not groups[i]&groups[j] for i in range(4) for j in range(i))
