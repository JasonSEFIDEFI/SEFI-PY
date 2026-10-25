"""Leading critical-curve radial problem, inherited from the frequency study."""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp,simpson
from scipy.interpolate import CubicHermiteSpline
HERE=Path(__file__).resolve().parent

def solve_limit(beta,seed=None,X=48.,tol=1e-10):
    x=np.linspace(0,X,1201)
    if seed is None:
        data=np.load(HERE/'reduced_seed.npz')
        profile=data['binding'] if beta>.2 else data['fold']
        curve=CubicHermiteSpline(data['x'],profile[0],profile[1])
        clipped=np.minimum(x,data['x'][-1])
        t=curve(clipped);v=curve(clipped,1)
        t[x>data['x'][-1]]=0;v[x>data['x'][-1]]=0
        y=np.array([t,v])
    else:
        clipped=np.minimum(x,seed.x[-1]);y=seed.sol(clipped)
        y[:,x>seed.x[-1]]=0
    y[1,0]=0
    def fun(x,y):
        t,v=y;sn=np.sin(t);cs=np.cos(t)
        a=3/16*(1+3*cs*cs);ap=-9/8*sn*cs
        up=3/4*sn*cs*(beta-2*sn*sn)
        return np.array([v,up/(a*beta)-ap*v*v/(2*a)])
    def bc(a,b):return np.array([a[1],b[0]])
    S=np.zeros((2,2));S[1,1]=-2
    sol=solve_bvp(fun,bc,x,y,S=S,tol=tol,bc_tol=1e-12,max_nodes=25000)
    assert sol.status==0,(sol.status,sol.message)
    assert .01<sol.y[0,0]<np.pi/2,sol.y[0,0]
    return sol

