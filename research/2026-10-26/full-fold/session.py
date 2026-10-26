"""Full-model charge folds and independently predicted critical corrections.

Dimensionless classical toy model; no particle/nucleus identification.
"""
from pathlib import Path
import sys,json,time,platform
import numpy as np
from scipy.integrate import solve_bvp,simpson
from scipy.optimize import minimize_scalar
import scipy

ROOT=Path(__file__).resolve().parent
AUDIT=next(p for p in [ROOT.parent/'sefi-math-audit-2026-10-25',
    ROOT.parents[1]/'2026-10-25/math-audit'] if p.exists())
sys.path.insert(0,str(AUDIT))
from continue_math import MC,central_solve,reduced_alpha
SOLVES=[]
SOLVE_FILE='all_full_solves.json'

def save(name,data):
    (ROOT/name).write_text(json.dumps(data,indent=2),encoding='utf-8')

def reduced_seed(center=1.105324367332945):
    previous=None
    for c in [.5,.7,.9,center]:
        previous,row=central_solve(c,320,1800,1e-10,previous)
        assert row['admissible'],row
    return previous

def predict(seed):
    reference=json.loads((AUDIT/'landmark_refinement.json').read_text())
    a0=reference['fold']['alpha']
    red,base=reduced_alpha(a0,seed,tol=1e-10)
    rows=[]
    for h in [2e-5,1e-5,5e-6]:
        _,plus=reduced_alpha(a0+h,red,tol=1e-10)
        _,minus=reduced_alpha(a0-h,red,tol=1e-10)
        S1p=(plus['S1']-minus['S1'])/(2*h)
        q1=2*MC*S1p+(3/64-a0)*base['norm']/(2*MC)
        q0pp=(plus['q_scaled']-2*base['q_scaled']+minus['q_scaled'])/h**2
        S1pp=(plus['S1']-2*base['S1']+minus['S1'])/h**2
        # q0'(a0)=0 at the leading fold, so S0''(a0)=0.
        q1p=2*MC*S1pp-base['norm']/(2*MC)
        rows.append(dict(h=h,q1=q1,q0_second_derivative=q0pp,
                         q1_derivative=q1p,alpha1=-q1p/q0pp))
    data=dict(alpha0=a0,q0=base['q_scaled'],q1=rows[-1]['q1'],
         alpha1=rows[-1]['alpha1'],derivative_refinements=rows,
         prediction='q_min(delta)=q0+delta*q1+O(delta^2); alpha_fold=alpha0+delta*alpha1+O(delta^2)',
         protocol='Derived exclusively from reduced solutions before this session full-model evaluation; not external preregistration')
    save('predictions_before_full.json',data)
    return data,seed

def full_center(delta,center,seed,L=320,tol=2e-9):
    lam=16-delta;vac=np.sqrt(3/lam);m=np.sqrt(1+12/lam)
    x=np.linspace(0,L,2000);t,v=seed.sol(x)
    y=np.array([np.sqrt(3)/2*np.sin(t),np.sqrt(3)/2*np.cos(t)*v,
                vac*np.cos(t),-vac*np.sin(t)*v])
    def rhs(x,y,p):
        f,fp,s,sp=y;a=p[0]
        return np.array([fp,(1-m*m+delta*a+f*f+4*s*s)*f/delta,
                         sp,(-3+4*f*f+lam*s*s)*s/delta])
    def bc(ya,yb,p):
        return np.array([ya[0]-np.sqrt(3)/2*np.sin(center),ya[1],ya[3],yb[0],yb[2]-vac])
    start=time.perf_counter()
    sol=solve_bvp(rhs,bc,x,y,p=seed.p,S=np.diag([0.,-2.,0.,-2.]),
                  tol=tol,max_nodes=40000)
    r=np.linspace(0,L,24001);f,fp,s,sp=sol.sol(r);a=float(sol.p[0])
    w=np.sqrt(m*m-delta*a);I=4*np.pi*simpson(r*r*f*f,x=r)
    square=lam/4*(s*s+(4*f*f-3)/lam)**2
    W=(square-delta*f**4/(4*lam))/delta+a*f*f/2
    T=4*np.pi*simpson(r*r*.5*(fp*fp+sp*sp),x=r)
    V=4*np.pi*simpson(r*r*W,x=r)
    admissible=bool(sol.success and 0<a<(m*m-(1+3/np.sqrt(lam)))/delta
        and f.min()>-1e-7 and s.min()>-1e-7 and np.max(fp)<1e-7 and np.min(sp)>-1e-7)
    row=dict(delta=delta,lam=lam,center=center,alpha=a,L=L,tol=tol,
        q_scaled=float(w*I),Q=float(w*I*delta**(-1.5)),omega=float(w),
        I_scaled=float(I),action_scaled=float(T+V),
        binding_scaled=float(T+V-a*w*I/(m+w)),
        virial_relative=float((T+3*V)/max(1,T)),
        success=bool(sol.success),admissible=admissible,
        residual=float(max(sol.rms_residuals)),nodes=len(sol.x),
        seconds=time.perf_counter()-start,message=sol.message)
    rr=np.linspace(.01,L,12001)
    yy=sol.sol(rr);dy=sol.sol(rr,1);expected=rhs(rr,yy,sol.p)
    expected[1]-=2*yy[1]/rr;expected[3]-=2*yy[3]/rr
    row['independent_ode_max_abs']=float(np.max(abs(dy-expected)))
    SOLVES.append(row);save(SOLVE_FILE,SOLVES)
    assert admissible and abs(row['virial_relative'])<2e-6,row
    return sol,row

def locate_fold(delta,seed,L=320,tol=2e-9):
    cache={}
    def evaluate(c):
        key=round(float(c),13)
        if key not in cache:cache[key]=full_center(delta,float(c),seed,L,tol)
        return cache[key]
    opt=minimize_scalar(lambda c:evaluate(c)[1]['q_scaled'],
        bounds=(.95,1.27),method='bounded',options={'xatol':2e-7})
    assert opt.success,opt
    sol,row=evaluate(opt.x)
    _,left=evaluate(opt.x-.001);_,right=evaluate(opt.x+.001)
    row=dict(row,charge_derivative_c=(right['q_scaled']-left['q_scaled'])/.002,
        charge_second_derivative_c=(right['q_scaled']-2*row['q_scaled']+left['q_scaled'])/1e-6,
        optimizer_evaluations=opt.nfev)
    assert left['q_scaled']>row['q_scaled'] and right['q_scaled']>row['q_scaled'],row
    return sol,row

def main():
    started=time.perf_counter();seed=reduced_seed();pred,seed=predict(seed)
    save('protocol.json',dict(deltas=[.5,.25,.125,.0625,.03125,.015625],
        full_domain=320,full_tolerance=2e-9,fold_bounds=[.95,1.27],
        center_parameter='f(0)=sqrt(3)/2*sin(center)',
        endpoint='Minimum scaled charge on the positive monotone spherical solution branch',
        limits='Not a global all-solution minimum, nonlinear stability or physical matter claim'))
    folds=[]
    for delta in [.5,.25,.125,.0625,.03125,.015625]:
        sol,row=locate_fold(delta,seed)
        row.update(q_leading_error=row['q_scaled']-pred['q0'],
            q_first_error=row['q_scaled']-(pred['q0']+delta*pred['q1']),
            alpha_first_error=row['alpha']-(pred['alpha0']+delta*pred['alpha1']))
        folds.append(row);save('folds.json',folds)
        np.savez_compressed(ROOT/f'fold_delta_{delta}.npz',rho=sol.x,fields=sol.y)
        print(json.dumps(row),flush=True)
    _,refine=locate_fold(.0625,seed,L=400,tol=2e-10)
    data=dict(predictions=pred,folds=folds,refinement=refine,
        effective_exponents=[float(np.log(b['Q']/a['Q'])/np.log(a['delta']/b['delta']))
                            for a,b in zip(folds[:-1],folds[1:])],
        full_solve_count=len(SOLVES),seconds=time.perf_counter()-started,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    save('results.json',data);print('COMPLETED',flush=True)

if __name__=='__main__':main()
