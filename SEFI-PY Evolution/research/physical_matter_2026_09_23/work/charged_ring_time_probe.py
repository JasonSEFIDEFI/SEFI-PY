"""Full two-field axisymmetric evolution in an exactly transformed moving frame.

T=t, Z=gamma(z-vt); Phi=exp(i Omega t) p(T,r,Z),
Sigma=exp(i nu gamma(t-vz)) s(T,r,Z). Both p and s evolve complex.
This is a limited perturbation test, not a stability theorem.
"""
import argparse,json
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--dt',type=float,default=.025);ap.add_argument('--amplitude',type=float,default=.001);ap.add_argument('--duration',type=float,default=100);ap.add_argument('--profile',default='charged_ring_refined');ap.add_argument('--name',default='')
a=ap.parse_args();root=Path(__file__).resolve().parents[1];d=np.load(root/'work'/f'{a.profile}.npz');r=d['r'];h=float(r[1]);z=np.r_[-d['z'][1:][::-1],d['z']]
if float(d.get('lambda_sigma',25))!=25 or int(d.get('carrier_azimuthal_winding',0))!=0:
    raise ValueError('This evolution probe supports only lambda=25 and zero carrier azimuthal winding.')
base=np.concatenate([d['u'][:,1:][:,::-1]-1j*d['w'][:,1:][:,::-1],d['u']+1j*d['w']],axis=1)
carrier=np.concatenate([d['s'][:,1:][:,::-1],d['s']],axis=1).astype(complex)
c=float(d['c']);nu=float(d['nu']);Omega=np.sqrt(2);v=c/np.sqrt(8+c*c);gamma=1/np.sqrt(1-v*v);adv=2*gamma*v
radius=json.loads((root/'work'/f'{a.profile}.json').read_text())['core_radii'][0]
rr,zz=np.meshgrid(r,z,indexing='ij');pert=a.amplitude*np.exp(-((rr-radius)**2+zz**2)/2);pert[-1,:]=pert[:,0]=pert[:,-1]=0
p=base+pert;s=carrier+pert;pv=np.zeros_like(p);sv=np.zeros_like(s)
wt=r.copy();wt[0]=h/8;wt[-1]=0
def dz(q):
    out=np.zeros_like(q);out[:,1:-1]=(q[:,2:]-q[:,:-2])/(2*h);return out
def lap(q):
    out=np.zeros_like(q)
    out[1:-1,1:-1]=(q[2:,1:-1]+q[:-2,1:-1]-2*q[1:-1,1:-1])/h**2+(q[2:,1:-1]-q[:-2,1:-1])/(2*h*r[1:-1,None])
    out[0,1:-1]=4*(q[1,1:-1]-q[0,1:-1])/h**2
    out[:-1,1:-1]+=(q[:-1,2:]+q[:-1,:-2]-2*q[:-1,1:-1])/h**2
    return out
def rhs(y):
    p,s,pv,sv=y;rho=abs(p)**2;eta=abs(s)**2
    pa=adv*dz(pv)-2j*Omega*pv+lap(p)+1j*c*dz(p)+(1-rho-4*eta)*p
    sa=adv*dz(sv)-2j*nu*gamma*sv+lap(s)+(3+nu*nu-4*rho-25*eta)*s
    for arr in [pa,sa]:arr[-1,:]=arr[:,0]=arr[:,-1]=0
    return np.array([pv,sv,pa,sa])
def integ(q):return float(2*np.pi*h*h*np.sum(wt[:,None]*q))
def norm(q):return float(np.sqrt(integ(abs(q)**2)))
def charge(s,sv):return integ(nu*abs(s)**2+np.imag(s.conj()*sv)/gamma-v*np.imag(s.conj()*dz(s)))
def winding(p):
    th=np.linspace(0,2*np.pi,513);rq=radius+1.5*np.cos(th);zq=1.5*np.sin(th);i=np.floor(rq/h).astype(int);j=np.floor((zq-z[0])/h).astype(int);aa=rq/h-i;bb=(zq-z[0])/h-j
    q=(1-aa)*(1-bb)*p[i,j]+aa*(1-bb)*p[i+1,j]+(1-aa)*bb*p[i,j+1]+aa*bb*p[i+1,j+1]
    return float((np.unwrap(np.angle(q))[-1]-np.unwrap(np.angle(q))[0])/(2*np.pi)),float(abs(q).min())
y=np.array([p,s,pv,sv]);hist=[];steps=round(a.duration/a.dt)
initial_rhs=rhs(np.array([base,carrier,pv,sv]));initial_res=float(abs(initial_rhs[2:]).max())
for it in range(steps+1):
    if it%max(1,steps//200)==0:
        p,s,pv,sv=y;wind,contour_min=winding(p)
        hist.append({'time':it*a.dt,'phi_distance':norm(p-base),'sigma_distance':norm(s-carrier),'carrier_charge':charge(s,sv),'winding':wind,'contour_min_amplitude':contour_min,'carrier_peak':float(abs(s).max())})
    if it==steps:break
    k1=rhs(y);k2=rhs(y+a.dt*k1/2);k3=rhs(y+a.dt*k2/2);k4=rhs(y+a.dt*k3);y+=a.dt*(k1+2*k2+2*k3+k4)/6
    if not np.isfinite(y).all():raise RuntimeError('Nonfinite solution')
report={'profile':a.profile,'dt':a.dt,'duration':steps*a.dt,'amplitude':a.amplitude,'background_acceleration_residual':initial_res,'initial':hist[0],'final':hist[-1],'max_relative_charge_drift':max(abs(t['carrier_charge']/hist[0]['carrier_charge']-1) for t in hist),'history':hist,'scope':'One axisymmetric perturbation family; finite reflecting cylinder, fixed spatial grid and finite duration. No full stability proof.'}
name=a.name or f'charged_time_dt{a.dt}_amp{a.amplitude}';np.savez(root/'work'/f'{name}.npz',phi=y[0],sigma=y[1],phi_velocity=y[2],sigma_velocity=y[3],r=r,z=z)
(root/'work'/f'{name}.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='history'},indent=2))
