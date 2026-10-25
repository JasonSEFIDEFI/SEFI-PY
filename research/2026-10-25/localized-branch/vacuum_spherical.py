"""Fixed-Phi-charge spherical relaxation over the true Sigma vacuum.

Separate lambda=12 research, not production equations or a vortex solution.
"""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent
try:from scipy.linalg import solve_banded
except ImportError:
    sys.path.insert(0,str(ROOT.parents[2].parent/'github-sync-test-deps'))
    from scipy.linalg import solve_banded

class Radial:
    def __init__(self,Q=1e5,lam=12.,L=60.,h=.4):
        self.Q=Q;self.lam=lam;self.L=L;self.h=h;n=round(L/h)
        self.r=(np.arange(n)+.5)*h;faces=np.arange(n+1)*h
        self.V=4*np.pi/3*(faces[1:]**3-faces[:-1]**3)
        self.links=4*np.pi*faces[1:-1]**2/h;self.outer=4*np.pi*L*L/(h/2)
        self.diag=np.r_[0,self.links]+np.r_[self.links,self.outer]
        self.vac=np.sqrt(3/lam);self.b=np.zeros(n);self.b[-1]=self.outer*self.vac

    def K(self,a):
        out=self.diag*a;out[:-1]-=self.links*a[1:];out[1:]-=self.links*a[:-1];return out

    def energy(self,f,s):
        I=np.sum(self.V*f*f);rho=f*f;eta=s*s
        K=.5*np.sum(self.links*(np.diff(f)**2+np.diff(s)**2))+.5*self.outer*(f[-1]**2+(s[-1]-self.vac)**2)
        U=np.sum(self.V*(rho*rho/4+rho/2+2*rho*eta+self.lam*(eta-3/self.lam)**2/4))
        return float(K+U+self.Q*self.Q/(2*I))

    def gradient(self,f,s):
        omega=self.Q/np.sum(self.V*f*f)
        gf=-self.K(f)/self.V+(omega*omega-1-f*f-4*s*s)*f
        gs=-(self.K(s)-self.b)/self.V+(3-4*f*f-self.lam*s*s)*s
        return gf,gs,float(omega)

    def precondition(self,g,mass):
        ab=np.zeros((3,len(g)));ab[1]=self.diag+self.V*mass
        ab[0,1:]=-self.links;ab[2,:-1]=-self.links
        return solve_banded((1,1),ab,self.V*g,check_finite=False)

    def initial(self):
        R=27.821717576746355;width=2.;a=.93
        x=np.clip((self.r-R+width)/(2*width),0,1);t=1-10*x**3+15*x**4-6*x**5
        return a*t,self.vac*(1-t)

def run(h=.4,wall=90,maximum=20000,L=60.,tag=None):
    tag=tag or f'vacuum_spherical_h{h}'.replace('.','p')
    if (ROOT/(tag+'.json')).exists():raise FileExistsError(tag)
    status=ROOT/(tag+'.status.json')
    status.write_text(json.dumps(dict(reason='running',h=h,L=L)))
    model=Radial(h=h,L=L);f,s=model.initial();E=model.energy(f,s);initialE=E
    t0=time.monotonic();reason='iteration_limit';history=[]
    for iteration in range(maximum):
        if (ROOT/'STOP_REQUEST').exists():reason='user_interrupted';break
        gf,gs,omega=model.gradient(f,s);res=float(max(abs(gf).max(),abs(gs).max()))
        if iteration%100==0:history.append(dict(iteration=iteration,energy=E,omega=omega,residual=res))
        if res<1e-7:reason='residual_tolerance';break
        if time.monotonic()-t0>wall:reason='wall_limit';break
        df=model.precondition(gf,1+3*f*f+4*s*s+omega*omega)
        ds=model.precondition(gs,1+abs(-3+4*f*f+3*model.lam*s*s))
        descent=np.sum(model.V*(gf*df+gs*ds));step=1.
        while step>2**-20:
            fn=f+step*df;sn=s+step*ds;En=model.energy(fn,sn)
            if En<=E-1e-4*step*descent:break
            step/=2
        if step<=2**-20:reason='line_search_stalled';break
        f,s,E=fn,sn,En
    gf,gs,omega=model.gradient(f,s);res=float(max(abs(gf).max(),abs(gs).max()))
    result=dict(scope='Discrete spherical fixed-charge critical-point search; not a vortex, existence theorem, or stability/observed-matter identification',
        Q=model.Q,lambda_sigma=model.lam,L=model.L,h=h,reason=reason,iterations=iteration,
        residual=res,initial_energy=initialE,energy=E,energy_per_charge=E/model.Q,omega=omega,
        free_phi_mass=float(np.sqrt(1+12/model.lam)),wall_seconds=time.monotonic()-t0,
        converged=res<1e-7,history=history)
    (ROOT/(tag+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(ROOT/(tag+'.npz'),r=model.r,f=f,s=s)
    status.write_text(json.dumps(dict(reason=reason,converged=res<1e-7,iterations=iteration)))
    print(json.dumps({k:v for k,v in result.items() if k!='history'}),flush=True)

if __name__=='__main__':
    # Check the discrete variational gradient before minimizing it.
    model=Radial();f,s=model.initial();rng=np.random.default_rng(4)
    df=.01*rng.normal(size=f.shape);ds=.01*rng.normal(size=s.shape)
    gf,gs,_=model.gradient(f,s);eps=1e-5
    finite=(model.energy(f+eps*df,s+eps*ds)-model.energy(f-eps*df,s-eps*ds))/(2*eps)
    analytic=-np.sum(model.V*(gf*df+gs*ds))
    assert abs(finite-analytic)/max(1,abs(analytic))<1e-6,(finite,analytic)
    print('Discrete gradient checked',flush=True)
    run(.4);run(.2)

