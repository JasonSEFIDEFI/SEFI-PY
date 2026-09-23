"""Coupled straight-vortex/core-condensate profiles for an explicit two-field extension."""
import argparse,json
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=400);ap.add_argument('--L',type=float,default=20);ap.add_argument('--omega',type=float,default=0)
a=ap.parse_args();root=Path(__file__).resolve().parents[1];n=a.n;L=a.L;h=L/n;r=np.arange(n+1)*h
source=np.genfromtxt(root/'outputs'/'vortex_profile.csv',delimiter=',',names=True)
kappa=4.;mu2=3.+a.omega**2;lam=25.;beta=1.;Delta=1.
nf=n-1;ns=n;size=nf+ns
A=np.zeros((size,size));b=np.zeros(size)
for i in range(1,n):
    row=i-1;A[row,row]=-2/h**2-1/r[i]**2
    if i>1:A[row,row-1]=(1-1/(2*i))/h**2
    if i<n-1:A[row,row+1]=(1+1/(2*i))/h**2
    else:b[row]=(1+1/(2*i))/h**2
for i in range(n):
    row=nf+i
    if i==0:A[row,row]=-4/h**2;A[row,row+1]=4/h**2
    else:
        A[row,row]=-2/h**2;A[row,row-1]=(1-1/(2*i))/h**2
        if i<n-1:A[row,row+1]=(1+1/(2*i))/h**2
def unpack(x):
    f=np.r_[0,x[:nf],1.];s=np.r_[x[nf:],0.];return f,s
def residual(x):
    f,s=unpack(x)
    return A@x+b+np.r_[((Delta-beta*f*f-kappa*s*s)*f)[1:n],((mu2-kappa*f*f-lam*s*s)*s)[:n]]
x=np.r_[np.interp(r[1:n],source['radius'],source['amplitude']),.32*np.exp(-r[:n]**2/3)]
ii=np.arange(size);fi=np.arange(nf);si=nf+np.arange(1,n)
for it in range(30):
    f,s=unpack(x);F=residual(x);err=float(abs(F).max())
    print(json.dumps({'iteration':it,'residual':err}),flush=True)
    J=A.copy();J[ii,ii]+=np.r_[(Delta-3*beta*f*f-kappa*s*s)[1:n],(mu2-kappa*f*f-3*lam*s*s)[:n]]
    cross=-2*kappa*f[1:n]*s[1:n];J[fi,si]+=cross;J[si,fi]+=cross
    if err<1e-9:break
    step=np.linalg.solve(J,-F);fac=1
    while fac>1/4096:
        test=x+fac*step
        if np.linalg.norm(residual(test))<np.linalg.norm(F):break
        fac/=2
    if fac<=1/4096:raise RuntimeError('Stalled')
    x=test
if err>=1e-9:raise RuntimeError('Not converged')
weights=np.r_[r[1:n]*h,np.r_[h*h/8,r[1:n]*h]]
H=-J*np.sqrt(weights[:,None]/weights[None,:])
assert np.max(abs(H-H.T))<1e-9
ev=np.linalg.eigvalsh(H)
f,s=unpack(x);norm=float(2*np.pi*np.trapezoid(r*s*s,r))
report={'n':n,'L':L,'h':h,'omega':a.omega,'carrier_charge_per_length':a.omega*norm,'tail_decay_rate':float(np.sqrt(1-a.omega**2)),'residual':err,'carrier_amplitude_at_core':float(s[0]),
        'carrier_transverse_norm':norm,'carrier_at_r10':float(np.interp(10,r,s)),
        'lowest_radial_real_amplitude_hessian_eigenvalues':ev[:5].tolist(),
        'scope':'Straight string, specified carrier phase frequency and zero longitudinal current. Radial real-amplitude Hessian of the frequency-constrained functional only; not full dynamical stability or a charged loop.'}
name=f'core_condensate_L{L:g}_n{n}'
if a.omega:name+=f'_omega{a.omega:g}'
np.savez(root/'work'/(name+'.npz'),r=r,f=f,s=s)
(root/'work'/(name+'.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
