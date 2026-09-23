"""Exact traveling two-field ansatz; finite-cylinder Newton test, not stability."""
import argparse,json
from pathlib import Path
import numpy as np
import solve_traveling_ring as ring
ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=24);ap.add_argument('--h',type=float,default=2/3);ap.add_argument('--c',type=float,default=.6);ap.add_argument('--nu',type=float,default=.3);ap.add_argument('--initial',required=True);ap.add_argument('--name',default='charged_ring_trial');ap.add_argument('--steps',type=int,default=35);ap.add_argument('--lam',type=float,default=25);ap.add_argument('--winding',type=int,default=0);ap.add_argument('--carrier-scale',type=float,default=1)
a=ap.parse_args();n=a.n;root=Path(__file__).resolve().parents[1];nu0=n*n;ns=nu0;nphi=n*n+n*(n-1);nt=nphi+ns
A,b,iu,iw=ring.setup(n,a.h,a.c)
D,_,_,_=ring.setup(n,a.h,0)
L=np.zeros((nt,nt));L[:nphi,:nphi]=A;L[nphi:,nphi:]=D[:nu0,:nu0]
bb=np.r_[b,np.zeros(ns)];rr,zz=np.meshgrid(np.arange(n+1)*a.h,np.arange(n+1)*a.h,indexing='ij')
if a.winding:
    for i in range(1,n):
        ids=nphi+i*n+np.arange(n);L[ids,ids]-=a.winding**2/(i*a.h)**2
xphi=ring.initial(n,a.h,5.5,a.initial)
old=np.load(a.initial)
if 's' in old:
    ro=old['r'];zo=old['z'];az=np.array([np.interp(zz[0],zo,row) for row in old['s']]);s=np.array([np.interp(rr[:,0],ro,az[:,j]) for j in range(n+1)]).T
else:s=.31*np.exp(-((rr-5.5)**2+zz**2)/3)
s*=a.carrier_scale;s[-1,:]=s[:,-1]=0
if a.winding:s[0,:]=0
x=np.r_[xphi,s[:n,:n].ravel()]
def unpack(x):
    u,w=ring.unpack(x[:nphi],n);s=np.zeros_like(u);s[:n,:n]=x[nphi:].reshape(n,n);return u,w,s
def residual(x):
    u,w,s=unpack(x);d=1-u*u-w*w-4*s*s
    out=L@x+bb+np.r_[ring.pack(d*u,d*w,n),((3+a.nu*a.nu-4*(u*u+w*w)-a.lam*s*s)*s)[:n,:n].ravel()]
    if a.winding:out[nphi:nphi+n]=x[nphi:nphi+n]
    return out
isall=nphi+np.arange(ns);isw=nphi+iu;ii=np.diag_indices(nt);history=[]
for it in range(a.steps):
    F=residual(x);err=float(abs(F).max());history.append(err);print(json.dumps({'iteration':it,'residual':err}),flush=True)
    if err<1e-9:break
    u,w,s=unpack(x);J=L.copy();J[ii]+=np.r_[(1-3*u*u-w*w-4*s*s)[:n,:n].ravel(),(1-u*u-3*w*w-4*s*s)[:n,1:n].ravel(),(3+a.nu*a.nu-4*(u*u+w*w)-3*a.lam*s*s)[:n,:n].ravel()]
    cross=(-2*u*w)[:n,1:n].ravel();J[iu,iw]+=cross;J[iw,iu]+=cross
    cross=(-8*u*s)[:n,:n].ravel();J[np.arange(nu0),isall]+=cross;J[isall,np.arange(nu0)]+=cross
    cross=(-8*w*s)[:n,1:n].ravel();J[iw,isw]+=cross;J[isw,iw]+=cross
    if a.winding:
        J[nphi:nphi+n,:]=0;axis=nphi+np.arange(n);J[axis,axis]=1
    dx=np.linalg.solve(J,-F);fac=1
    while fac>1/4096:
        test=x+fac*dx
        if np.linalg.norm(residual(test))<np.linalg.norm(F)*(1-1e-4*fac):break
        fac/=2
    if fac<=1/4096:break
    x=test
u,w,s=unpack(x);r=rr[:,0];z=zz[0];K=float(4*np.pi*np.trapezoid(np.trapezoid(rr*s*s,z,axis=1),r))
report={'n':n,'h':a.h,'L':n*a.h,'c':a.c,'nu':a.nu,'lambda_sigma':a.lam,'carrier_azimuthal_winding':a.winding,'v':float(a.c/np.sqrt(8+a.c*a.c)),'residual':float(abs(residual(x)).max()),'carrier_max':float(s.max()),'carrier_min':float(s.min()),'carrier_norm_comoving':K,'carrier_charge':a.nu*K,'core_radii':ring.diagnostics(x[:nphi],n,a.h,a.c)['core_radii'],'history':history,'scope':'Finite-domain coupled classical profile search; convergence and nonzero winding must be checked; no dynamical stability test.'}
np.savez(root/'work'/f'{a.name}.npz',r=r,z=z,u=u,w=w,s=s,c=a.c,nu=a.nu,lambda_sigma=a.lam,carrier_azimuthal_winding=a.winding)
(root/'work'/f'{a.name}.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
