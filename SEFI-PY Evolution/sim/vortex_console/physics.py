"""Research equations, independent of SEFI-PY engine. NumPy only.

Coordinates T=t, Z=gamma(z-vt). Phi=exp(i Omega t)p(T,r,Z),
Sigma=exp(i nu gamma(t-vz))s(T,r,Z). Speed of the underlying
Minkowski metric is 1. The evolution implemented here is axisymmetric
and restricted to lambda=25 and carrier azimuthal winding N=0.
"""
import math
import numpy as np
import ring_operator as ring

MODEL_VERSION='coupled-ring-axisymmetric-v0.1'

DEFAULTS=dict(c=.6,nu=.3,lam=25.,N=0,L=16.,n=32,seed=5.36,
              epsilon=.001,duration=10.,dt=.025,save_every=.25,
              wall_seconds=120.,storage_mb=50.)

def validate(raw, evolution=False):
    p=DEFAULTS.copy()
    for key in p:
        if key in raw:
            if isinstance(raw[key],bool):raise ValueError(f'{key} must be numeric')
            try:p[key]=float(raw[key])
            except (ValueError,TypeError):raise ValueError(f'{key} must be numeric')
        if not math.isfinite(p[key]):raise ValueError(f'{key} must be finite')
    bounds={'c':(0,1.2),'nu':(0,.95),'lam':(9.01,40),'N':(0,6),
            'L':(12,32),'n':(16,40),'seed':(2,10),'epsilon':(0,.02),
            'duration':(.05,100),'dt':(.002,.025),'save_every':(.025,10),
            'wall_seconds':(5,300),'storage_mb':(2,100)}
    for k,(lo,hi) in bounds.items():
        if not lo<=p[k]<=hi:raise ValueError(f'{k} must be between {lo} and {hi}')
    for k in ['n','N']:
        if int(p[k])!=p[k]:raise ValueError(f'{k} must be an integer')
        p[k]=int(p[k])
    if p['seed']>p['L']-3:raise ValueError('Initial ring must be at least 3 units inside the boundary')
    if p['dt']>p['L']/p['n']/10:raise ValueError('Time step too large for selected spatial grid')
    if p['save_every']<p['dt']:raise ValueError('Frame interval must be at least one time step')
    if math.ceil(p['duration']/p['save_every'])>200:raise ValueError('Choose a frame interval giving at most 200 frames')
    if evolution and (p['lam']!=25 or p['N']!=0):
        raise ValueError('Evolution is validated here only for lambda=25 and N=0. Other settings are profile searches only.')
    return p

def load_profile(path, p=None):
    d=np.load(path,allow_pickle=False)
    fields={k:d[k].copy() for k in ['r','z','u','w','s']}
    settings=DEFAULTS.copy();settings.update(c=float(d['c']),nu=float(d['nu']),
        L=float(d['r'][-1]),n=len(d['r'])-1,
        lam=float(d.get('lambda_sigma',25)),N=int(d.get('carrier_azimuthal_winding',0)))
    if p:settings.update(p)
    return fields,settings

def roots(f):
    r=f['r'];u=f['u'][:,0];out=[]
    for i in range(len(r)-1):
        if u[i]*u[i+1]<0:out.append(float(r[i]-u[i]*(r[i+1]-r[i])/(u[i+1]-u[i])))
    return out

def solve(p, reference, callback=lambda *x:None, cancelled=lambda:False):
    n=p['n'];h=p['L']/n;nu0=n*n;nphi=n*n+n*(n-1);size=nphi+nu0
    A,b,iu,iw=ring.setup(n,h,p['c']);D,_,_,_=ring.setup(n,h,0)
    op=np.zeros((size,size));op[:nphi,:nphi]=A;op[nphi:,nphi:]=D[:nu0,:nu0];bb=np.r_[b,np.zeros(nu0)]
    r=np.arange(n+1)*h;rr,zz=np.meshgrid(r,r,indexing='ij')
    # User-selected seed radius scales only the initial guess, never the result.
    source_radius=roots(reference)[0];scale=p['seed']/source_radius
    def interp(a):
        temp=np.array([np.interp(r,reference['z']*scale,row) for row in a])
        return np.array([np.interp(r,reference['r']*scale,temp[:,j]) for j in range(n+1)]).T
    u=interp(reference['u']);w=interp(reference['w']);s=interp(reference['s'])
    u[-1,:]=u[:,-1]=1;w[-1,:]=w[:,-1]=w[:,0]=0;s[-1,:]=s[:,-1]=0
    if p['N']:
        s[0,:]=0
        for i in range(1,n):
            idx=nphi+i*n+np.arange(n);op[idx,idx]-=p['N']**2/(i*h)**2
    x=np.r_[ring.pack(u,w,n),s[:n,:n].ravel()]
    def unpack(x):
        u,w=ring.unpack(x[:nphi],n);s=np.zeros_like(u);s[:n,:n]=x[nphi:].reshape(n,n);return u,w,s
    def residual(x):
        u,w,s=unpack(x);delta=1-u*u-w*w-4*s*s
        F=op@x+bb+np.r_[ring.pack(delta*u,delta*w,n),((3+p['nu']**2-4*(u*u+w*w)-p['lam']*s*s)*s)[:n,:n].ravel()]
        if p['N']:F[nphi:nphi+n]=x[nphi:nphi+n]
        return F
    hist=[];reason='iteration limit'
    for iteration in range(30):
        if cancelled():raise InterruptedError('Stopped during profile solve')
        F=residual(x);err=float(abs(F).max());hist.append(err);callback(iteration,err)
        if err<1e-9:reason='converged';break
        u,w,s=unpack(x);J=op.copy();ii=np.diag_indices(size)
        J[ii]+=np.r_[(1-3*u*u-w*w-4*s*s)[:n,:n].ravel(),(1-u*u-3*w*w-4*s*s)[:n,1:n].ravel(),(3+p['nu']**2-4*(u*u+w*w)-3*p['lam']*s*s)[:n,:n].ravel()]
        cross=(-2*u*w)[:n,1:n].ravel();J[iu,iw]+=cross;J[iw,iu]+=cross
        ids=nphi+np.arange(nu0);cross=(-8*u*s)[:n,:n].ravel();J[np.arange(nu0),ids]+=cross;J[ids,np.arange(nu0)]+=cross
        ids=nphi+iu;cross=(-8*w*s)[:n,1:n].ravel();J[iw,ids]+=cross;J[ids,iw]+=cross
        if p['N']:J[nphi:nphi+n,:]=0;ids=nphi+np.arange(n);J[ids,ids]=1
        dx=np.linalg.solve(J,-F);step=1
        while step>1/4096:
            test=x+step*dx
            if np.linalg.norm(residual(test))<np.linalg.norm(F)*(1-1e-4*step):break
            step/=2
        if step<=1/4096:reason='line search stalled';break
        x=test
    u,w,s=unpack(x);f=dict(r=r,z=r.copy(),u=u,w=w,s=s)
    err=float(abs(residual(x)).max());rs=roots(f)
    return f,dict(residual=err,converged=err<1e-8,has_ring=bool(rs),radius=rs[0] if rs else None,
                  carrier_peak=float(abs(s).max()),history=hist,reason=reason)

def full_fields(f):
    p=np.concatenate([f['u'][:,1:][:,::-1]-1j*f['w'][:,1:][:,::-1],f['u']+1j*f['w']],axis=1)
    s=np.concatenate([f['s'][:,1:][:,::-1],f['s']],axis=1).astype(complex)
    return p,s,np.r_[-f['z'][1:][::-1],f['z']]

class Evolution:
    def __init__(self,f,p):
        self.p=validate(p,True);self.f=f;self.r=f['r'];self.h=self.r[1];self.base,self.carrier,self.z=full_fields(f)
        self.radius=roots(f)[0];self.v=p['c']/np.sqrt(8+p['c']**2);self.gamma=1/np.sqrt(1-self.v**2)
        rr,zz=np.meshgrid(self.r,self.z,indexing='ij');pert=p['epsilon']*np.exp(-((rr-self.radius)**2+zz**2)/2)
        pert[-1,:]=pert[:,0]=pert[:,-1]=0
        self.y=np.array([self.base+pert,self.carrier+pert,np.zeros_like(self.base),np.zeros_like(self.carrier)])
        self.time=0.;self.steps=0;self.wt=self.r.copy();self.wt[0]=self.h/8;self.wt[-1]=0
        self.q0=self.charge();self.initial_rhs=float(abs(self.rhs(np.array([self.base,self.carrier,np.zeros_like(self.base),np.zeros_like(self.carrier)]))[2:]).max())
    def dz(self,q):
        out=np.zeros_like(q);out[:,1:-1]=(q[:,2:]-q[:,:-2])/(2*self.h);return out
    def lap(self,q):
        h=self.h;out=np.zeros_like(q)
        out[1:-1,1:-1]=(q[2:,1:-1]+q[:-2,1:-1]-2*q[1:-1,1:-1])/h**2+(q[2:,1:-1]-q[:-2,1:-1])/(2*h*self.r[1:-1,None])
        out[0,1:-1]=4*(q[1,1:-1]-q[0,1:-1])/h**2
        out[:-1,1:-1]+=(q[:-1,2:]+q[:-1,:-2]-2*q[:-1,1:-1])/h**2;return out
    def rhs(self,y):
        p,s,pv,sv=y;c=self.p['c'];nu=self.p['nu'];rho=abs(p)**2;eta=abs(s)**2;adv=2*self.gamma*self.v
        pa=adv*self.dz(pv)-2j*np.sqrt(2)*pv+self.lap(p)+1j*c*self.dz(p)+(1-rho-4*eta)*p
        sa=adv*self.dz(sv)-2j*nu*self.gamma*sv+self.lap(s)+(3+nu*nu-4*rho-25*eta)*s
        for a in [pa,sa]:a[-1,:]=a[:,0]=a[:,-1]=0
        return np.array([pv,sv,pa,sa])
    def step(self,dt=None):
        dt=self.p['dt'] if dt is None else dt;y=self.y
        k1=self.rhs(y);k2=self.rhs(y+dt*k1/2);k3=self.rhs(y+dt*k2/2);k4=self.rhs(y+dt*k3)
        self.y=y+dt*(k1+2*k2+2*k3+k4)/6;self.time+=dt;self.steps+=1
        if not np.isfinite(self.y).all():raise FloatingPointError('Nonfinite field: integration stopped')
    def integ(self,q):return float(2*np.pi*self.h**2*np.sum(self.wt[:,None]*q))
    def norm(self,q):return np.sqrt(self.integ(abs(q)**2))
    def charge(self):
        s,sv=self.y[1],self.y[3]
        return self.integ(self.p['nu']*abs(s)**2+np.imag(s.conj()*sv)/self.gamma-self.v*np.imag(s.conj()*self.dz(s)))
    def diagnostics(self):
        p,s=self.y[:2];phase=np.angle(np.sum(self.wt[:,None]*self.carrier.conj()*s))
        th=np.linspace(0,2*np.pi,513);rq=self.radius+1.5*np.cos(th);zq=1.5*np.sin(th)
        if rq.min()<0 or rq.max()>self.r[-1]:wind=None;contour=None
        else:
            i=np.floor(rq/self.h).astype(int);j=np.floor((zq-self.z[0])/self.h).astype(int);aa=rq/self.h-i;bb=(zq-self.z[0])/self.h-j
            q=(1-aa)*(1-bb)*p[i,j]+aa*(1-bb)*p[i+1,j]+(1-aa)*bb*p[i,j+1]+aa*bb*p[i+1,j+1]
            contour=float(abs(q).min());wind=float((np.unwrap(np.angle(q))[-1]-np.unwrap(np.angle(q))[0])/(2*np.pi)) if contour>1e-6 else None
        Q=self.charge()
        return dict(time=self.time,charge=Q,charge_drift=abs(Q-self.q0)/max(abs(self.q0),1e-12),
                    phi_distance=float(self.norm(p-self.base)),sigma_distance=float(self.norm(s-self.carrier)),
                    aligned_distance=float(self.norm(s*np.exp(-1j*phase)-self.carrier)),phase_shift=float(phase),
                    winding=wind,contour_min=contour,peak=float(abs(s).max()))

def frame(r,z,p,s,diag=None):
    return dict(r=np.round(r,6).tolist(),z=np.round(z,6).tolist(),
                phi=np.round(abs(p)**2,6).tolist(),sigma=np.round(abs(s)**2,6).tolist(),
                phase=np.round(np.angle(p),6).tolist(),diagnostics=diag or {})
