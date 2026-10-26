import numpy as np
import sympy as sp
from continue_math import symbolic_audit

def test_exact_symbolic_identities():
    out=symbolic_audit()
    assert out['square_identity'] and out['first_correction_local_density']
    assert out['foundations_linear_schur']

def test_fixed_charge_second_variation():
    # Includes nonstationary states: checks the charge constraint, not a fitted eigenvalue.
    rng=np.random.default_rng(109)
    f=rng.uniform(.1,1,37);eta=rng.normal(size=37);weights=rng.uniform(.2,2,37)
    Q=5.
    def grad(x):
        I=np.sum(weights*x*x)
        return -Q*Q*weights*x/I**2
    I=np.sum(weights*f*f)
    Heta=-Q*Q*weights*eta/I**2+4*Q*Q*weights*f*np.sum(weights*f*eta)/I**3
    finite=(grad(f+1e-6*eta)-grad(f-1e-6*eta))/(2e-6)
    assert np.allclose(finite,Heta,rtol=1e-8,atol=1e-10)

def test_global_frequency_lower_bound():
    # Independent numerical minimization checks the constrained s^2>=0 branch formula.
    from scipy.optimize import minimize_scalar
    for lam in [2.,8.,12.,15.9,16.,20.]:
        m2=1+12/lam
        def ratio(r):
            t=max(0.,(3-4*r)/lam)
            U=r*r/4+r/2+2*r*t+lam*(t-3/lam)**2/4
            return 2*U/r
        res=minimize_scalar(ratio,bounds=(1e-9,20),method='bounded')
        expected=1+3/np.sqrt(lam) if lam<16 else m2
        assert abs(res.fun-expected)<1e-5

def test_foundations_general_schur_complement():
    A,S,W,kA,kS,kW,C,D=sp.symbols('A S W kA kS kW C D')
    mu=sp.Function('mu')(A);nu=sp.Function('nu')(A)
    V=(kA*A*A+kS*S*S+kW*W*W)/2+mu*C*S+nu*D*W
    K=sp.hessian(V,(A,S,W)).subs({S:-mu*C/kS,W:-nu*D/kW})
    schur=K[0,0]-K[0,1]**2/kS-K[0,2]**2/kW
    expected=kA-C*C*(sp.diff(mu,A)**2+mu*sp.diff(mu,A,2))/kS-D*D*(sp.diff(nu,A)**2+nu*sp.diff(nu,A,2))/kW
    assert sp.simplify(schur-expected)==0

def test_general_quartic_square():
    a,b,c,lam,v,f,s=sp.symbols('a b c lam v f s')
    U=a*f**4/4+b*f*f/2+c*f*f*s*s/2+lam*(s*s-v*v)**2/4
    square=lam/4*(s*s-v*v+c*f*f/lam)**2+(a*lam-c*c)*f**4/(4*lam)
    assert sp.expand(U-(b+c*v*v)*f*f/2-square)==0

def test_radial_two_component_ground_state_identity():
    r=sp.symbols('r',positive=True)
    p=sp.Function('p')(r);q=sp.Function('q')(r)
    z=sp.Function('z')(r);w=sp.Function('w')(r);C=sp.Function('C')(r)
    lap=lambda x:sp.diff(x,r,2)+2*sp.diff(x,r)/r
    V1=lap(p)/p+C*q/p;V2=lap(q)/q+C*p/q
    lhs=r*r*(sp.diff(p*z,r)**2+sp.diff(q*w,r)**2+
                  V1*p*p*z*z+V2*q*q*w*w-2*C*p*q*z*w)
    positive=r*r*(p*p*sp.diff(z,r)**2+q*q*sp.diff(w,r)**2+C*p*q*(z-w)**2)
    boundary=sp.diff(r*r*(p*sp.diff(p,r)*z*z+q*sp.diff(q,r)*w*w),r)
    assert sp.simplify(lhs-positive-boundary)==0
