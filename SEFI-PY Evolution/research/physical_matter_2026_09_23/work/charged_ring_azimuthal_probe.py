"""Linear initial-value probe of one azimuthal sector of the charged ring.

Four real field perturbations times cos(m phi); no background carrier winding.
This samples an initial condition, not the entire mode spectrum.
"""
import argparse,json
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--m',type=int,default=2);ap.add_argument('--dt',type=float,default=.025);ap.add_argument('--duration',type=float,default=100);ap.add_argument('--profile',default='charged_ring_refined')
a=ap.parse_args();root=Path(__file__).resolve().parents[1];d=np.load(root/'work'/f'{a.profile}.npz');r=d['r'];h=float(r[1]);z=np.r_[-d['z'][1:][::-1],d['z']];m=a.m
if float(d.get('lambda_sigma',25))!=25 or int(d.get('carrier_azimuthal_winding',0))!=0:
    raise ValueError('This linear probe supports only lambda=25 and zero carrier azimuthal winding.')
u=np.concatenate([d['u'][:,1:][:,::-1],d['u']],axis=1);w=np.concatenate([-d['w'][:,1:][:,::-1],d['w']],axis=1);s=np.concatenate([d['s'][:,1:][:,::-1],d['s']],axis=1)
c=float(d['c']);nu=float(d['nu']);Omega=np.sqrt(2);v=c/np.sqrt(8+c*c);gamma=1/np.sqrt(1-v*v);adv=2*gamma*v
rr,zz=np.meshgrid(r,z,indexing='ij');radius=json.loads((root/'work'/f'{a.profile}.json').read_text())['core_radii'][0]
env=np.exp(-((rr-radius)**2+zz**2)/2);env[-1,:]=env[:,0]=env[:,-1]=0
if m:env[0,:]=0
q=np.array([env, .3*env*np.cos(zz), .7*env, .2*env*np.sin(zz)]);vel=np.zeros_like(q);y=np.array([q,vel])
wt=r.copy();wt[0]=h/8;wt[-1]=0
def dz(q):
    out=np.zeros_like(q);out[:,:,1:-1]=(q[:,:,2:]-q[:,:,:-2])/(2*h);return out
def lap(q):
    out=np.zeros_like(q);out[:,1:-1,1:-1]=(q[:,2:,1:-1]+q[:,:-2,1:-1]-2*q[:,1:-1,1:-1])/h**2+(q[:,2:,1:-1]-q[:,:-2,1:-1])/(2*h*r[None,1:-1,None])
    if not m:out[:,0,1:-1]=4*(q[:,1,1:-1]-q[:,0,1:-1])/h**2
    out[:,:-1,1:-1]+=(q[:,:-1,2:]+q[:,:-1,:-2]-2*q[:,:-1,1:-1])/h**2
    if m:out[:,1:-1,:]-=m*m*q[:,1:-1,:]/r[None,1:-1,None]**2
    return out
def rhs(y):
    q,vv=y;aa,bb,dd,ee=q;rho=u*u+w*w
    acc=adv*dz(vv)+lap(q);der=dz(q)
    acc[0]+=2*Omega*vv[1]-c*der[1]+(1-3*u*u-w*w-4*s*s)*aa-2*u*w*bb-8*u*s*dd
    acc[1]+=-2*Omega*vv[0]+c*der[0]+(1-u*u-3*w*w-4*s*s)*bb-2*u*w*aa-8*w*s*dd
    acc[2]+=2*nu*gamma*vv[3]+(3+nu*nu-4*rho-75*s*s)*dd-8*u*s*aa-8*w*s*bb
    acc[3]+=-2*nu*gamma*vv[2]+(3+nu*nu-4*rho-25*s*s)*ee
    acc[:,-1,:]=acc[:,:,0]=acc[:,:,-1]=0
    if m:acc[:,0,:]=0
    return np.array([vv,acc])
def norm(q):return float(np.sqrt(2*np.pi*h*h*np.sum(wt[None,:,None]*q*q)))
hist=[];steps=round(a.duration/a.dt);initial=norm(q)
for it in range(steps+1):
    if it%max(1,steps//200)==0:hist.append({'time':it*a.dt,'field_norm_ratio':norm(y[0])/initial,'velocity_norm_ratio':norm(y[1])/initial})
    if it==steps:break
    k1=rhs(y);k2=rhs(y+a.dt*k1/2);k3=rhs(y+a.dt*k2/2);k4=rhs(y+a.dt*k3);y+=a.dt*(k1+2*k2+2*k3+k4)/6
report={'m':m,'dt':a.dt,'duration':steps*a.dt,'profile':a.profile,'final_field_norm_ratio':hist[-1]['field_norm_ratio'],'max_field_norm_ratio':max(x['field_norm_ratio'] for x in hist),'history':hist,'scope':'One initial perturbation in one azimuthal Fourier sector of the linearized equations, finite grid and reflecting domain. Not an eigenvalue search or nonlinear 3D evolution.'}
name=f'charged_azimuthal_m{m}_dt{a.dt}';np.savez(root/'work'/f'{name}.npz',perturbations=y[0],velocities=y[1]);(root/'work'/f'{name}.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='history'},indent=2))
