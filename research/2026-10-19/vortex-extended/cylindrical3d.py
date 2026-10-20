"""Experimental nonlinear angular extension of the preserved N=0 equations.

Independent research prototype. No changes to the production console.
Finite-difference r,z and Fourier theta on a periodic angular grid.
"""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent
PRIOR=ROOT.parent/'vortex-survival'
sys.path.insert(0,str(PRIOR))
from study import physics,CONSOLE,cores

class Cylinder:
    def __init__(self,reference='refined',ntheta=16,dt=.025):
        if ntheta<1 or (ntheta!=1 and ntheta%2):raise ValueError('Use one or an even number of angles')
        self.f,self.p=physics.load_profile(CONSOLE/f'data/charged_ring_{reference}.npz')
        self.p.update(epsilon=0,dt=dt)
        self.axis=physics.Evolution(self.f,self.p)
        e=self.axis;self.r=e.r;self.z=e.z;self.h=e.h;self.ntheta=ntheta
        self.theta=np.arange(ntheta)*2*np.pi/ntheta
        self.m=np.fft.fftfreq(ntheta,1/ntheta)
        self.y=np.repeat(e.y[:,:,None,:],ntheta,axis=2)
        self.base=self.y.copy();self.time=0.;self.dt=dt

    def dz(self,q):
        out=np.zeros_like(q);out[...,1:-1]=(q[...,2:]-q[...,:-2])/(2*self.h);return out

    def lap(self,q):
        out=np.zeros_like(q);h=self.h
        out[1:-1,:,1:-1]=(q[2:,:,1:-1]+q[:-2,:,1:-1]-2*q[1:-1,:,1:-1])/h**2+(q[2:,:,1:-1]-q[:-2,:,1:-1])/(2*h*self.r[1:-1,None,None])
        out[0,:,1:-1]=4*(q[1,:,1:-1].mean(axis=0)-q[0,:,1:-1].mean(axis=0))/h**2
        out[:-1,:,1:-1]+=(q[:-1,:,2:]+q[:-1,:,:-2]-2*q[:-1,:,1:-1])/h**2
        if self.ntheta>1:
            angular=np.fft.ifft(-self.m[None,:,None]**2*np.fft.fft(q,axis=1),axis=1)
            out[1:-1,:,1:-1]+=angular[1:-1,:,1:-1]/self.r[1:-1,None,None]**2
        return out

    def rhs(self,y):
        p,s,pv,sv=y;e=self.axis
        rho=abs(p)**2;eta=abs(s)**2
        pa=2*e.gamma*e.v*self.dz(pv)-2j*np.sqrt(2)*pv+self.lap(p)+1j*self.p['c']*self.dz(p)+(1-rho-4*eta)*p
        sa=2*e.gamma*e.v*self.dz(sv)-2j*self.p['nu']*e.gamma*sv+self.lap(s)+(3+self.p['nu']**2-4*rho-25*eta)*s
        for a in [pa,sa]:
            a[-1,:,:]=a[:,:,0]=a[:,:,-1]=0
            a[0,:,:]=a[0,:,:].mean(axis=0,keepdims=True)
        return np.array([pv,sv,pa,sa])

    def step(self,dt=None):
        dt=self.dt if dt is None else dt;y=self.y
        a=self.rhs(y);b=self.rhs(y+dt*a/2);c=self.rhs(y+dt*b/2);d=self.rhs(y+dt*c)
        self.y=y+dt*(a+2*b+2*c+d)/6;self.time+=dt
        if not np.isfinite(self.y).all():raise FloatingPointError('Nonfinite 3D field')

    def integ(self,q):return self.axis.integ(np.mean(q,axis=1))

    def norm(self,q):return float(np.sqrt(self.integ(abs(q)**2)))

    def energy(self):
        p,s,pv,sv=self.y;d=p-1;rho=abs(p)**2;eta=abs(s)**2
        return (.5*self.integ(abs(pv)**2+abs(sv)**2)-.5*self.integ(np.real(d.conj()*self.lap(d)+s.conj()*self.lap(s)))
            +self.integ((rho-1)**2/4+25*eta**2/4+(4*rho-3-self.p['nu']**2)*eta/2)
            +.5*self.p['c']*self.integ(np.imag(d.conj()*self.dz(d))))

    def charge(self,field=1):
        q,v=self.y[field],self.y[field+2];e=self.axis
        omega=np.sqrt(2) if field==0 else self.p['nu']*e.gamma
        density=(omega*(abs(q)**2-(1 if field==0 else 0))+np.imag(q.conj()*v))/e.gamma-e.v*np.imag(q.conj()*self.dz(q))
        return self.integ(density)

    def perturb(self,epsilon=.01,m=2,kind='amplitude'):
        rr,tt,zz=np.meshgrid(self.r,self.theta,self.z,indexing='ij')
        g=(rr/self.axis.radius)**m*np.exp(-((rr-self.axis.radius)**2+zz**2)/2)*np.cos(m*tt)
        g[0,:,:]=g[-1,:,:]=g[:,:,0]=g[:,:,-1]=0
        if kind=='amplitude':self.y[0]+=epsilon*g;self.y[1]+=epsilon*(.4+.3j)*g
        elif kind=='phase':self.y[0]*=np.exp(1j*epsilon*g);self.y[1]*=np.exp(-1j*epsilon*g)
        else:raise ValueError('Unknown perturbation')

    def diagnostics(self,core_scan=True):
        p,s=self.y[:2];meanp=p.mean(axis=1,keepdims=True);means=s.mean(axis=1,keepdims=True)
        rho=abs(s)**2;I=self.integ(rho)
        rr=self.r[:,None,None];zz=self.z[None,None,:]
        result=dict(time=self.time,energy=self.energy(),carrier_charge=self.charge(),phi_excess_charge=self.charge(0),
            phi_distance=self.norm(p-self.base[0]),sigma_distance=self.norm(s-self.base[1]),
            angular_asymmetry=float(np.sqrt(self.norm(p-meanp)**2+self.norm(s-means)**2)),
            carrier_rms_radius=float(np.sqrt(self.integ((rr**2+zz**2)*rho)/I)),
            carrier_edge_fraction=self.integ(rho*((rr>.8*self.r[-1])|(abs(zz)>.8*self.z[-1])))/I,
            max_amplitude=float(max(abs(p).max(),abs(s).max())))
        if core_scan:
            counts=[];ambiguous=0
            for j in range(self.ntheta):
                found,bad=cores(p[:,j,:],self.r,self.z);counts.append(len(found));ambiguous+=bad
            result.update(angles_with_detected_core=sum(n>0 for n in counts),minimum_core_count=min(counts),ambiguous_cells=ambiguous)
        # Angular power on the non-axisymmetric original-field deviation.
        power=np.mean(abs(np.fft.fft(p-meanp,axis=1)/self.ntheta)**2,axis=(0,2))
        total=float(power.sum());high=abs(self.m)>=self.ntheta/3
        result['high_angular_power_fraction']=float(power[high].sum()/total) if total>1e-30 else 0.
        return result

def run_case(name,T,ntheta=16,dt=.025,epsilon=.01,m=2,reference='refined',wall_limit=180,initial=None):
    status_path=ROOT/(name+'.status.json')
    def status(value):
        temporary=status_path.with_suffix('.tmp')
        temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(status_path)
    if (ROOT/(name+'.json')).exists():raise FileExistsError('Use a new run ID; preserve the prior result')
    status(dict(name=name,status='running',T_requested=T,ntheta=ntheta,dt=dt,epsilon=epsilon,m=m,reference=reference))
    c=Cylinder(reference,ntheta,dt);c.perturb(epsilon,m)
    if initial is not None:
        old=np.load(initial);assert old['final'].shape==c.y.shape;c.y=old['final'].copy()
    t0=time.monotonic();samples=[c.diagnostics()];initial_y=c.y.copy();reason='completed'
    for i in range(round(T/dt)):
        if (ROOT/'STOP_REQUEST').exists():reason='user_interrupted';break
        if time.monotonic()-t0>wall_limit:reason='wall_limit';break
        try:c.step()
        except KeyboardInterrupt:reason='user_interrupted';break
        except FloatingPointError:reason='nonfinite';break
        if (i+1)%round(1/dt)==0:
            d=c.diagnostics();samples.append(d)
            if d['max_amplitude']>10:reason='amplitude_guard';break
    if samples[-1]['time']!=c.time:samples.append(c.diagnostics())
    q0=samples[0]['carrier_charge'];E0=samples[0]['energy']
    result=dict(name=name,T_requested=T,ntheta=ntheta,dt=dt,epsilon=epsilon,m=m,reference=reference,
        reason=reason,wall_seconds=time.monotonic()-t0,initial_source=str(initial) if initial else 'saved_reference_plus_perturbation',
        min_angles_with_core=min(x['angles_with_detected_core'] for x in samples),
        maximum_asymmetry=max(x['angular_asymmetry'] for x in samples),
        asymmetry_gain=max(x['angular_asymmetry'] for x in samples)/max(samples[0]['angular_asymmetry'],1e-30),
        max_relative_energy_drift=max(abs(x['energy']-E0) for x in samples)/max(abs(E0),1e-30),
        max_relative_carrier_charge_drift=max(abs(x['carrier_charge']-q0) for x in samples)/max(abs(q0),1e-30),
        samples=samples)
    (ROOT/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(ROOT/(name+'.npz'),initial=initial_y,final=c.y,r=c.r,z=c.z,theta=c.theta)
    status(dict(name=name,status=reason,T_requested=T,achieved_time=c.time,checkpoint_saved=True,eligible_for_completed_comparison=reason=='completed' and abs(c.time-T)<1e-8))
    print(json.dumps({k:v for k,v in result.items() if k!='samples'}),flush=True)
    return result
