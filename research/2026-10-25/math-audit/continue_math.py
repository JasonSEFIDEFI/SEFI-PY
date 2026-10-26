"""Branch audit and first asymptotic correction in the dimensionless toy model."""
from pathlib import Path
import sys,json,time,platform
import numpy as np
from scipy.integrate import solve_bvp,simpson
from scipy.optimize import minimize_scalar,brentq
import scipy
ROOT=Path(__file__).resolve().parent
OLD=next(p for p in [ROOT.parent/'sefi-critical-binding-2026-10-23',
    ROOT.parent/'critical-binding',ROOT.parents[1]/'2026-10-23/critical-binding'] if p.exists())
sys.path.insert(0,str(OLD))
from critical import solve as central_solve,identities
MC=np.sqrt(7/4)
RUNS=[]

def geom(t):
    return 3*(1+3*np.cos(t)**2)/16,-9*np.sin(t)*np.cos(t)/8

def vp(t,a):
    return 3*a*np.sin(t)*np.cos(t)/4-9*np.sin(t)**3*np.cos(t)/256

def integrate_profile(sol,a,L):
    r=np.linspace(0,L,24001);t,v=sol.sol(r);G,Gp=geom(t)
    lap=(vp(t,a)-.5*Gp*v*v)/G
    T=4*np.pi*simpson(r*r*.5*G*v*v,x=r)
    V=4*np.pi*simpson(r*r*(3*a*np.sin(t)**2/8-9*np.sin(t)**4/1024),x=r)
    norm=4*np.pi*simpson(r*r*.75*np.sin(t)**2,x=r)
    laps=-np.sqrt(3)/4*(np.sin(t)*lap+np.cos(t)*v*v)
    first=4*np.pi*simpson(r*r*(3*np.sin(t)**2*v*v/512
                           -9*np.sin(t)**4/16384-laps*laps/(12*np.cos(t)**2)),x=r)
    return dict(alpha=float(a),center=float(t[0]),q_scaled=float(MC*norm),
                S0=float(T+V),S1=float(first),norm=float(norm),
                binding_scaled=float(T+V-a*norm/2),
                virial_relative=float((T+3*V)/max(1,T)),
                tail_slope=float(v[-1]))

def reduced_alpha(a,guess,L=320,tol=1e-9):
    x=np.linspace(0,L,1800)
    def fun(x,y):
        t,v=y;G,Gp=geom(t)
        return np.array([v,(vp(t,a)-.5*Gp*v*v)/G])
    def bc(ya,yb):return np.array([ya[1],yb[0]])
    start=time.perf_counter()
    sol=solve_bvp(fun,bc,x,guess.sol(x),S=np.diag([0.,-2.]),tol=tol,max_nodes=30000)
    row=integrate_profile(sol,a,L)
    row.update(success=bool(sol.success),max_residual=float(max(sol.rms_residuals)),
               seconds=time.perf_counter()-start,L=L,tol=tol)
    assert sol.success and .05<row['center']<np.pi/2,row
    assert abs(row['virial_relative'])<2e-6,row
    RUNS.append(dict(kind='reduced_alpha',**row))
    return sol,row

def full_alpha(delta,a,guess,L=320,tol=2e-8):
    lam=16-delta;vac=np.sqrt(3/lam);m=np.sqrt(1+12/lam)
    w=np.sqrt(m*m-delta*a);x=np.linspace(0,L,1800)
    t,v=guess.sol(x)
    initial=np.array([np.sqrt(3)/2*np.sin(t),np.sqrt(3)/2*np.cos(t)*v,
                      vac*np.cos(t),-vac*np.sin(t)*v])
    def fun(x,y):
        f,fp,s,sp=y
        return np.array([fp,((1-w*w+f*f+4*s*s)*f)/delta,
                         sp,((-3+4*f*f+lam*s*s)*s)/delta])
    def bc(ya,yb):return np.array([ya[1],ya[3],yb[0],yb[2]-vac])
    start=time.perf_counter()
    sol=solve_bvp(fun,bc,x,initial,S=np.diag([0.,-2.,0.,-2.]),tol=tol,max_nodes=40000)
    r=np.linspace(0,L,24001);f,fp,s,sp=sol.sol(r)
    I=4*np.pi*simpson(r*r*f*f,x=r)
    square=lam/4*(s*s+(4*f*f-3)/lam)**2
    W=(square-delta*f**4/(4*lam))/delta+a*f*f/2
    T=4*np.pi*simpson(r*r*.5*(fp*fp+sp*sp),x=r)
    V=4*np.pi*simpson(r*r*W,x=r)
    row=dict(delta=delta,alpha=a,L=L,tol=tol,success=bool(sol.success),
             q_scaled=float(w*I),action_scaled=float(T+V),
             binding_scaled=float(T+V-a*w*I/(m+w)),
             central_f=float(f[0]),central_s=float(s[0]),
             virial_relative=float((T+3*V)/max(1,T)),
             max_residual=float(max(sol.rms_residuals)),nodes=len(sol.x),
             seconds=time.perf_counter()-start)
    assert sol.success and f.min()>-1e-7 and s.min()>-1e-7 and f[0]>.1,row
    assert np.max(fp)<1e-7 and abs(row['virial_relative'])<2e-5,row
    RUNS.append(dict(kind='full_alpha',**row))
    return sol,row

def symbolic_audit():
    import sympy as sp
    out=identities()
    O,A,S,W,kO,kA,kS,kW,c1,c2,C,D=sp.symbols('O A S W kO kA kS kW c1 c2 C D',real=True)
    V=(kO*O**2+kA*A**2+kS*S**2+kW*W**2)/2+c1*A*C*S+c2*A*D*W
    reduced=V.subs({S:-c1*A*C/kS,W:-c2*A*D/kW})
    assert sp.expand(reduced-(kO*O**2+(kA-c1**2*C**2/kS-c2**2*D**2/kW)*A**2)/2)==0
    H=sp.hessian(V,(A,S,W))
    assert sp.simplify(H.det()-kS*kW*(kA-c1**2*C**2/kS-c2**2*D**2/kW))==0
    d,a,t,v,e,ep=sp.symbols('delta alpha theta v eta eta_prime',real=True)
    lam=16-d;f=sp.sqrt(3)/2*sp.sin(t);s=sp.sqrt(3/lam)*sp.cos(t)+d*e
    fp=sp.diff(f,t)*v;ss=-sp.sqrt(3/lam)*sp.sin(t)*v+d*ep
    z=s*s+(4*f*f-3)/lam
    square1=4*sp.diff(z,d).subs(d,0)**2
    other=(fp*fp+ss*ss)/2-f**4/(4*lam)+a*f*f/2
    first=sp.simplify(sp.diff(other,d).subs(d,0)+square1)
    expected=3*sp.sin(t)**2*v*v/512-sp.sqrt(3)/4*sp.sin(t)*v*ep+3*sp.cos(t)**2*e*e-9*sp.sin(t)**4/16384
    assert sp.trigsimp(first-expected)==0
    out.update(foundations_linear_schur=True,first_correction_local_density=True)
    return out

def main():
    start=time.perf_counter();syms=symbolic_audit()
    print(json.dumps(syms),flush=True)
    warm,warmrow=central_solve(.5,240,1200,1e-7)
    assert warmrow['admissible'],warmrow
    class Extended:
        p=warm.p
        def sol(self,r):
            y=warm.sol(np.minimum(r,240))
            y[:,r>240]=0.
            return y
    previous=Extended();seed=None;scan=[]
    for c in np.arange(.5,1.551,.025):
        sol,row=central_solve(float(c),480,2400,1e-8,previous)
        assert row['admissible'],row
        previous=sol
        if abs(c-1.45)<1e-6:seed=sol
        check=integrate_profile(sol,row['alpha'],480)
        assert abs(check['virial_relative'])<2e-6,check
        scan.append(check)
    (ROOT/'branch_scan.json').write_text(json.dumps(scan,indent=2))
    cache={}
    def at(c):
        key=round(float(c),12)
        if key not in cache:
            sol,row=central_solve(float(c),320,1800,1e-9,seed)
            assert row['admissible'],row
            cache[key]=integrate_profile(sol,row['alpha'],320)
        return cache[key]
    opt=minimize_scalar(lambda c:at(c)['q_scaled'],bounds=(.95,1.25),method='bounded',options={'xatol':1e-8})
    crossing=brentq(lambda c:at(c)['binding_scaled'],1.3,1.5,xtol=1e-10)
    fold=at(opt.x);bound=at(crossing)
    print(json.dumps(dict(fold=fold,energetic_crossing=bound)),flush=True)
    a=float(seed.p[0]);red,base=reduced_alpha(a,seed)
    derivative_checks=[]
    for h in [1e-5,5e-6]:
        _,plus=reduced_alpha(a+h,red);_,minus=reduced_alpha(a-h,red)
        ds=(plus['S0']-minus['S0'])/(2*h)
        derivative_checks.append(dict(h=h,S0_derivative=ds,expected=base['norm']/2,
             relative_envelope_error=ds/(base['norm']/2)-1,S1_derivative=(plus['S1']-minus['S1'])/(2*h)))
    ds1=derivative_checks[-1]['S1_derivative'];q0=base['q_scaled']
    q1=2*MC*ds1+(3/64-a)/MC*base['norm']/2
    predictions=dict(alpha=a,S0=base['S0'],S1=base['S1'],q0=q0,q1=q1,
       statement='Derived from reduced solutions before full-model solves; no fitting to full-model outputs.')
    (ROOT/'predictions_before_full.json').write_text(json.dumps(predictions,indent=2))
    full=[]
    for d in [.125,.0625,.03125,.015625,.0078125]:
        sol,row=full_alpha(d,a,red)
        row.update(S_leading_error=row['action_scaled']-base['S0'],
          S_first_order_error=row['action_scaled']-(base['S0']+d*base['S1']),
          q_leading_error=row['q_scaled']-q0,q_first_order_error=row['q_scaled']-(q0+d*q1))
        full.append(row);print(json.dumps(row),flush=True)
        np.savez_compressed(ROOT/f'full_fixed_alpha_delta_{d}.npz',rho=sol.x,fields=sol.y)
        (ROOT/'partial_full.json').write_text(json.dumps(full,indent=2))
    _,refine=full_alpha(.03125,a,red,L=400,tol=2e-9)
    findings=dict(symbolic=syms,scan=scan,fold=fold,energetic_crossing=bound,
       derivative_checks=derivative_checks,predictions=predictions,full=full,
       full_refinement=refine,runs=RUNS,
       environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
       total_seconds=time.perf_counter()-start)
    (ROOT/'results.json').write_text(json.dumps(findings,indent=2))
    print('COMPLETED',flush=True)

if __name__=='__main__':main()
