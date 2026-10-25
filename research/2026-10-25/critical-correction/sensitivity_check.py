"""Linear profile sensitivity check of the first-order charge coefficients."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp,simpson
from reduced_solver import solve_limit
from derive_correction import BETA_TURN,BETA_BIND,MC,evaluate
HERE=Path(__file__).resolve().parent

def expressions():
    t,p,beta=sp.symbols('t p beta',real=True)
    sn,cs=sp.sin(t),sp.cos(t)
    a=sp.Rational(3,16)*(1+3*cs**2);b=sp.diff(a,t)/(2*a)
    h=sp.Rational(3,4)*sn*cs/a
    g=h*(1-2*sn**2/beta)
    lap_t=h*(beta-2*sn**2)-b*beta*p*p
    lap_s=-sp.sqrt(3)/4*(sn*lap_t+cs*beta*p*p)
    rho=sp.Rational(3,4)*sn**2
    l1=sn**2*beta*p*p/4+sp.Rational(4,9)*rho*rho-lap_s**2/(12*cs**2)
    coefficients=sp.lambdify((t,p,beta),(b,sp.diff(b,t),sp.diff(g,t),sp.diff(g,beta)),modules='numpy',cse=True)
    gradients=sp.lambdify((t,p,beta),(l1,sp.diff(l1,t),sp.diff(l1,p),sp.diff(l1,beta)),modules='numpy',cse=True)
    return coefficients,gradients

def run(beta):
    theta=solve_limit(beta);coef,grad=expressions()
    def fun(x,y):
        t,p=theta.sol(x);b,bp,gt,gb=coef(t,p,beta)
        u,v=y
        return np.array([v,-2*b*p*v-(bp*p*p-gt)*u+gb])
    def bc(a,b):return np.array([a[1],b[0]])
    S=np.zeros((2,2));S[1,1]=-2
    x=np.linspace(0,48,1201)
    sensitivity=solve_bvp(fun,bc,x,np.zeros((2,len(x))),S=S,tol=1e-10,bc_tol=1e-12,max_nodes=25000)
    assert sensitivity.status==0,sensitivity.message
    x=np.linspace(0,48,24001);t,p=theta.sol(x);u,v=sensitivity.sol(x)
    l1,lt,lp,lb=grad(t,p,beta)
    integrate=lambda density:float(simpson(4*np.pi*x*x*density,x=x)/beta**1.5)
    A1=integrate(l1)
    A1p=-1.5*A1/beta+integrate(lt*u+lp*v+lb)
    I0=integrate(.75*np.sin(t)**2)
    I0p=-1.5*I0/beta+integrate(1.5*np.sin(t)*np.cos(t)*u)
    Q1=2*MC*A1p+(2-beta)*I0/(2*MC)
    checks=[]
    for step in [.0001,.00005]:
        lo=solve_limit(beta-step);hi=solve_limit(beta+step)
        fd=(hi.sol(x)[0]-lo.sol(x)[0])/(2*step)
        checks.append(dict(step=step,max_absolute_profile_sensitivity_error=float(max(abs(fd-u)))))
    result=dict(beta=beta,A1=A1,A1_beta_derivative=A1p,I0_beta_derivative=I0p,
        Q0_beta_derivative=MC*I0p,Q1=Q1,
        max_sensitivity_residual=float(max(sensitivity.rms_residuals)),finite_difference_profile_checks=checks)
    return result

if __name__=='__main__':
    out=HERE/'sensitivity_results.json'
    if out.exists():raise FileExistsError(out)
    result=dict(scope=__doc__,runs=[run(BETA_TURN),run(BETA_BIND)])
    predictions=json.loads((HERE/'predictions.json').read_text())
    result['turn_C1_difference_from_frozen_prediction']=result['runs'][0]['Q1']-predictions['coefficients']['turn_C1']
    assert abs(result['turn_C1_difference_from_frozen_prediction'])<1e-6,result
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

