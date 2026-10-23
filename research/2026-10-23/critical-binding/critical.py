"""Exact identities and exploratory critical-limit radial BVP; dimensionless toy model."""
from pathlib import Path
import json,time
import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp, simpson

ROOT=Path(__file__).resolve().parent

def identities():
    f,s,lam=sp.symbols('f s lam', real=True)
    U=f**4/4+f**2/2+2*f**2*s**2+lam*(s**2-3/lam)**2/4
    m2=1+12/lam
    square=lam/4*(s**2+(4*f**2-3)/lam)**2+(lam-16)*f**4/(4*lam)
    assert sp.simplify(U-m2*f**2/2-square)==0
    d,a,t=sp.symbols('delta alpha theta', real=True)
    fv=sp.sqrt(3)/2*sp.sin(t);sv=sp.sqrt(3)/4*sp.cos(t)
    potential=(U-(m2-d*a)*f**2/2).subs({f:fv,s:sv,lam:16-d})
    leading=sp.simplify(sp.diff(potential,d).subs(d,0))
    expected=3*a*sp.sin(t)**2/8-9*sp.sin(t)**4/1024
    assert sp.trigsimp(leading-expected)==0
    metric=sp.simplify(sp.diff(fv,t)**2+sp.diff(sv,t)**2)
    assert sp.trigsimp(metric-3*(1+3*sp.cos(t)**2)/16)==0
    return dict(square_identity=True,critical_potential=str(leading),metric=str(metric))

def solve(center,L=240,n=1200,tol=1e-6,guess=None):
    r=np.linspace(0,L,n)
    if guess is None:
        width=18+8/(max(.05,1.58-center))
        theta=center*np.exp(-(r/width)**2)
        y=np.array([theta,-2*r/width**2*theta]); p=[.015]
    else:
        y=guess.sol(r);p=guess.p
    def fun(r,y,p):
        t,v=y;sn=np.sin(t);cs=np.cos(t)
        G=3*(1+3*cs**2)/16;Gp=-9*sn*cs/8
        Vp=3*p[0]*sn*cs/4-9*sn**3*cs/256
        return np.array([v,(Vp-.5*Gp*v*v)/G])
    def bc(ya,yb,p):return np.array([ya[0]-center,ya[1],yb[0]])
    start=time.perf_counter()
    sol=solve_bvp(fun,bc,r,y,p=p,S=np.array([[0.,0.],[0.,-2.]]),tol=tol,max_nodes=30000)
    x=np.linspace(0,L,10001);t,v=sol.sol(x);a=float(sol.p[0]);sn=np.sin(t)
    G=3*(1+3*np.cos(t)**2)/16;V=3*a*sn**2/8-9*sn**4/1024
    norm=4*np.pi*simpson(x*x*.75*sn**2,x=x)
    action=4*np.pi*simpson(x*x*(.5*G*v*v+V),x=x)
    q=np.sqrt(7/4)*norm
    physical=bool(sol.success and 0<a<3/128 and t.min()>-1e-7 and t.max()<np.pi/2 and np.max(v)<1e-7)
    row=dict(center=center,L=L,tol=tol,alpha=a,q_scaled=q,action_scaled=action,
             binding_scaled=action-a*q/(2*np.sqrt(7/4)),success=bool(sol.success),
             admissible=physical,max_rms_residual=float(max(sol.rms_residuals)),
             nodes=len(sol.x),seconds=time.perf_counter()-start,message=sol.message)
    return sol,row

if __name__=='__main__':
    print(identities(),flush=True)
    rows=[];previous=None
    for c in [.3,.5,.7,.9,1.1,1.3,1.45,1.52,1.55,1.565]:
        sol,row=solve(c,guess=previous)
        rows.append(row);print(json.dumps(row),flush=True)
        if row['admissible']:previous=sol
        np.savez_compressed(ROOT/f'profile_{c}.npz',r=sol.x,theta=sol.y[0],derivative=sol.y[1],alpha=sol.p)
        (ROOT/'reduced_results.json').write_text(json.dumps(dict(identities=identities(),runs=rows),indent=2))
