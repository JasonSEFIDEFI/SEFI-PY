"""Fixed-carrier-charge gradient search; fictitious time, NOT physical evolution.

Keeps the Phi background frequency fixed; does not constrain excess Phi charge.
Axisymmetric Phi=u+iw, Sigma=s exp(i nu t+iN azimuth), nu=Q/integral(s^2).
"""
import argparse,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
MODEL_VERSION='fixed-charge-rest-v0.1'
DEFAULTS=dict(Q=1000.,lam=12.,N=1,L=24.,h=.5,steps=18000,every=150,wall_seconds=120.,storage_mb=50.,initial='ring')
SAVED=['lambda12_Q100','lambda25_Q100','Q1000_wide','Q1000_refined']

def validate(raw):
    p=DEFAULTS.copy()
    bounds=dict(Q=(1,3000),lam=(9.01,40),N=(0,3),L=(16,32),h=(.4,1),steps=(1,30000),every=(1,30000),wall_seconds=(5,300),storage_mb=(2,100))
    for key,(lo,hi) in bounds.items():
        value=raw.get(key,p[key])
        if isinstance(value,bool):raise ValueError(key+' must be numeric')
        try:value=float(value)
        except (ValueError,TypeError):raise ValueError(key+' must be numeric')
        if not np.isfinite(value) or not lo<=value<=hi:raise ValueError(f'{key} must be finite and between {lo} and {hi}')
        if key in ['N','steps','every']:
            if value!=int(value):raise ValueError(key+' must be an integer')
            value=int(value)
        p[key]=value
    if abs(p['L']/p['h']-round(p['L']/p['h']))>1e-8:raise ValueError('L/h must be an integer')
    if np.ceil(p['steps']/p['every'])>200:raise ValueError('Choose a save interval allowing at most 200 intervals')
    p['initial']=raw.get('initial','ring')
    if p['initial'] not in ['ring']+SAVED:raise ValueError('Unknown initial profile')
    return p

def load_initial(search,path):
    with np.load(path,allow_pickle=False) as old:
        for k in range(3):
            temp=np.array([np.interp(search.z,old['z'],row) for row in old['y'][k]])
            search.y[k]=np.array([np.interp(search.r,old['r'],temp[:,j]) for j in range(len(search.z))]).T
    search.boundary(search.y)

class Search:
    def __init__(self,L=20,h=.5,Q=100,lam=12,N=1):
        h=float(h)
        self.h=h;self.Q=Q;self.lam=lam;self.N=N;self.r=np.arange(round(L/h)+1,dtype=float)*h;self.z=np.arange(-round(L/h),round(L/h)+1,dtype=float)*h
        self.wt=self.r.copy();self.wt[0]=h/8;self.wt[-1]=self.r[-1]/2;self.weight=2*np.pi*h*h*self.wt[:,None]
        r,z=np.meshgrid(self.r,self.z,indexing='ij');R=5.36
        q=r*r+z*z-R*R+2j*R*z
        phi=q/np.sqrt(abs(q)**2+4*R*R)
        s=.5*(r/R)**N*np.exp(-((r-R)**2+z*z)/8)
        self.y=np.array([phi.real,phi.imag,s]);self.boundary(self.y)
    def boundary(self,y):
        y[0,-1,:]=y[0,:,0]=y[0,:,-1]=1
        y[1:,-1,:]=y[1:,:,0]=y[1:,:,-1]=0
        if self.N:y[2,0,:]=0
    def integral(self,a):return float(np.sum(self.weight*a))
    def lap(self,a):
        h=self.h;out=np.zeros_like(a);out[1:-1,1:-1]=(a[2:,1:-1]+a[:-2,1:-1]-2*a[1:-1,1:-1])/h**2+(a[2:,1:-1]-a[:-2,1:-1])/(2*h*self.r[1:-1,None]);out[0,1:-1]=4*(a[1,1:-1]-a[0,1:-1])/h**2
        out[:-1,1:-1]+=(a[:-1,2:]+a[:-1,:-2]-2*a[:-1,1:-1])/h**2;return out
    def terms(self,y):
        u,w,s=y;x=u*u+w*w;ss=s*s;I=self.integral(ss);nu=self.Q/I
        K=0.
        for a in y:
            K+=np.pi*np.sum((self.r[:-1]+self.h/2)[:,None]*np.diff(a,axis=0)**2)
            K+=np.pi*np.sum(self.wt[:,None]*np.diff(a,axis=1)**2)
        angular=np.zeros_like(s);angular[1:]=self.N**2*ss[1:]/self.r[1:,None]**2
        K+=self.integral(angular)/2
        U=self.integral((x-1)**2/4+self.lam*ss*ss/4+(4*x-3)*ss/2);C=self.Q**2/(2*I)
        return dict(K=float(K),U=U,C=C,energy=K+U+C,I=I,nu=nu,virial=K+3*U-3*C)
    def force(self,y):
        u,w,s=y;rho=u*u+w*w;ss=s*s;nu=self.Q/self.integral(ss)
        out=np.array([self.lap(u)+(1-rho-4*ss)*u,self.lap(w)+(1-rho-4*ss)*w,self.lap(s)+(3+nu*nu-4*rho-self.lam*ss)*s])
        out[2,1:]-=self.N**2*s[1:]/self.r[1:,None]**2
        out[:,-1,:]=out[:,:,0]=out[:,:,-1]=0
        if self.N:out[2,0,:]=0
        return out
    def diagnostics(self):
        d=self.terms(self.y);u,w,s=self.y;j=len(self.z)//2;rs=[]
        for i in range(len(self.r)-1):
            if u[i,j]*u[i+1,j]<0:rs.append(float(self.r[i]-u[i,j]*self.h/(u[i+1,j]-u[i,j])))
        r,z=np.meshgrid(self.r,self.z,indexing='ij');ss=s*s
        d.update(residual=float(abs(self.force(self.y)).max()),core_radii=rs,carrier_peak=float(abs(s).max()),carrier_rms_radius=float(np.sqrt(self.integral(ss*(r*r+z*z))/d['I'])),outer_charge_fraction=self.integral(ss*((r>.8*self.r[-1])|(abs(z)>.8*self.r[-1])))/d['I'],frequency_localized=d['nu']<1)
        deficit=self.integral(u*u+w*w-1)
        d.update(original_minimum=float(np.sqrt(u*u+w*w).min()),charge=self.Q,background_charge_deficit=float(np.sqrt(2)*deficit),physical_relative_energy=d['energy']+2*deficit,energy_per_charge=d['energy']/self.Q)
        return d
    def descent_step(self):
        f=self.force(self.y);energy=self.terms(self.y)['energy'];dt=min(.015,self.h*self.h/12)
        while dt>1e-10:
            trial=self.y+dt*f;new=self.terms(trial)['energy']
            if np.isfinite(new) and new<energy:
                self.y=trial;return dt
            dt/=2
        raise ArithmeticError('Relaxation line search stalled; no solution claimed')
    def gradient_check(self):
        rng=np.random.default_rng(103);a=rng.normal(size=self.y.shape);a[:,-1,:]=a[:,:,0]=a[:,:,-1]=0
        if self.N:a[2,0,:]=0
        eps=1e-6;fd=(self.terms(self.y+eps*a)['energy']-self.terms(self.y-eps*a)['energy'])/(2*eps);exact=-float(np.sum(self.weight*self.force(self.y)*a))
        rel=abs(fd-exact)/max(1,abs(exact));assert rel<1e-6,(fd,exact,rel);return rel
    def run(self,steps=12000,seconds=55):
        start=time.monotonic();hist=[];reason='iteration limit';step=min(.015,self.h*self.h/12);check=self.gradient_check()
        for it in range(steps+1):
            f=self.force(self.y);res=float(abs(f).max())
            if it%500==0:hist.append(dict(iteration=it,**self.diagnostics()))
            if res<1e-5:reason='residual tolerance';break
            if time.monotonic()-start>seconds:reason='runtime limit';break
            if it==steps:break
            energy=self.terms(self.y)['energy'];dt=step
            while dt>1e-10:
                trial=self.y+dt*f;new=self.terms(trial)['energy']
                if np.isfinite(new) and new<energy:break
                dt/=2
            if dt<=1e-10:reason='line search stalled';break
            self.y=trial
        final=self.diagnostics();hist.append(dict(iteration=it,**final))
        return dict(parameters=dict(L=float(self.r[-1]),h=self.h,Q=self.Q,lam=self.lam,N=self.N),gradient_check_relative_error=check,iterations=it,seconds=time.monotonic()-start,reason=reason,final=final,history=hist,scope='Axisymmetric fixed-carrier-charge descent at fixed background frequency. Fictitious relaxation, not physical evolution; no full stability claim.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Fixed-charge relaxation, not physical evolution')
    for key in ['Q','lam','L','h','wall_seconds']:parser.add_argument('--'+key,type=float,default=DEFAULTS[key])
    for key in ['N','steps','every']:parser.add_argument('--'+key,type=int,default=DEFAULTS[key])
    parser.add_argument('--initial',choices=['ring']+SAVED,default='ring');args=parser.parse_args();p=validate(vars(args))
    search=Search(**{k:p[k] for k in ['L','h','Q','lam','N']});base=Path(__file__).resolve().parent
    if p['initial']!='ring':load_initial(search,base/f'data/fixed_charge_{p["initial"]}.npz')
    result=search.run(p['steps'],p['wall_seconds']);result['model_version']=MODEL_VERSION;result['initial']=p['initial']
    from uuid import uuid4
    out=base/'runs'/('cli_'+uuid4().hex);out.mkdir(parents=True)
    (out/'result.json').write_text(json.dumps(result,indent=2));np.savez_compressed(out/'profile.npz',r=search.r,z=search.z,y=search.y,parameters=json.dumps(result['parameters']),model_version=MODEL_VERSION)
    print(json.dumps(dict(output=str(out),reason=result['reason'],final=result['final']),indent=2))


