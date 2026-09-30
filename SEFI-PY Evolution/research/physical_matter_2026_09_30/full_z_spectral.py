"""Full-Z quadratic Bogoliubov solver, isolated from the v0.1 implementation.

Finite cylinder, fixed outer perturbations, regular cylindrical axis. The
stationary background has parity; eigenvectors have no imposed axial parity.
No positive-real eigenvalue is accepted before all calibration/refinement gates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse as sp
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import eigs, spsolve


def operators(n, h, ell=0):
    """Node grid r=0..L, Z=-L..L, with outer Dirichlet nodes eliminated.

    For ell=0 the axis Laplacian is 4(f(h)-f(0))/h²; ell!=0
    eliminates the axis (regular angular harmonics vanish there).
    """
    i0 = 0 if ell == 0 else 1
    nz = 2*n-1
    ids = np.arange((n-i0)*nz).reshape(n-i0, nz)
    rows, cols, vals = [], [], []
    dr, dc, dv = [], [], []
    bc = np.zeros(ids.size)
    def add(row, i, j, value):
        if i == n or j in (-n, n):
            bc[row] += value  # inhomogeneous unit background for u only
        elif i >= i0:
            rows.append(row); cols.append(ids[i-i0, j+n-1]); vals.append(value)
    for i in range(i0, n):
        for j in range(-n+1, n):
            k = ids[i-i0, j+n-1]
            if i == 0:
                add(k, 1, j, 4/h**2); add(k, 0, j, -4/h**2)
            else:
                add(k, i-1, j, (1-1/(2*i))/h**2)
                add(k, i+1, j, (1+1/(2*i))/h**2)
                add(k, i, j, -2/h**2-ell**2/(i*h)**2)
            add(k, i, j-1, 1/h**2); add(k, i, j+1, 1/h**2)
            add(k, i, j, -2/h**2)
            for jj, v in ((j-1, -1/(2*h)), (j+1, 1/(2*h))):
                if -n < jj < n:
                    dr.append(k); dc.append(ids[i-i0, jj+n-1]); dv.append(v)
    shape = (ids.size, ids.size)
    return sp.csr_matrix((vals,(rows,cols)),shape=shape), sp.csr_matrix((dv,(dr,dc)),shape=shape), bc


def load_profile(path, n=None, h=None):
    d = np.load(path)
    oldr, oldz = d['r'], d['z']
    if oldz[0] < 0:
        raise ValueError('Input must be a retained half-domain profile')
    if n is None: n = len(oldr)-1
    if h is None: h = float(oldr[1]-oldr[0])
    r, z = np.arange(n+1)*h, np.arange(-n,n+1)*h
    rr, zz = np.meshgrid(r,z,indexing='ij')
    pts = np.stack((rr.ravel(),abs(zz).ravel()),axis=1)
    fields=[]
    for key, bg in (('u',1.),('w',0.),('s',0.)):
        f=RegularGridInterpolator((oldr,oldz),d[key],bounds_error=False,fill_value=bg)(pts).reshape(rr.shape)
        if key=='w': f *= np.sign(zz)
        f[-1,:]=bg; f[:,0]=bg; f[:,-1]=bg
        fields.append(f)
    return n,h,fields,float(d['c']),float(d['nu']),float(d.get('lambda_sigma',25)),int(d.get('carrier_azimuthal_winding',0))


def stationary(fields,n,h,c,nu,lam,steps=20):
    """Independent sparse Newton solve of the N=0 full-domain equations."""
    lap,dz,bc=operators(n,h)
    sl=np.s_[:n,1:-1]
    u,w,s=[f[sl].ravel().copy() for f in fields]
    size=u.size
    derivative_bc=np.zeros((n,2*n-1))
    derivative_bc[:,0]=-1/(2*h); derivative_bc[:,-1]=1/(2*h)
    derivative_bc=derivative_bc.ravel()
    def residual(x):
        u,w,s=np.split(x,3); density=u*u+w*w
        a=1-density-4*s*s
        return np.r_[lap@u+bc-c*(dz@w)+a*u,
                     lap@w+c*(dz@u+derivative_bc)+a*w,
                     lap@s+(3+nu*nu-4*density-lam*s*s)*s]
    x=np.r_[u,w,s]; hist=[]
    for it in range(steps+1):
        f=residual(x); err=float(abs(f).max()); hist.append(err)
        if err<1e-13 or it==steps: break
        u,w,s=np.split(x,3)
        diag=sp.diags
        uw=diag(-2*u*w); us=diag(-8*u*s); ws=diag(-8*w*s)
        j=sp.bmat([[lap+diag(1-3*u*u-w*w-4*s*s),uw-c*dz,us],
                   [uw+c*dz,lap+diag(1-u*u-3*w*w-4*s*s),ws],
                   [us,ws,lap+diag(3+nu*nu-4*(u*u+w*w)-3*lam*s*s)]],format='csc')
        dx=spsolve(j,-f); factor=1.
        while factor>=1/4096:
            test=x+factor*dx
            if np.linalg.norm(residual(test))<np.linalg.norm(f)*(1-1e-4*factor): break
            factor/=2
        if factor<1/4096: break
        x=test
    for f, v in zip(fields,np.split(x,3)): f[sl]=v.reshape(n,2*n-1)
    return fields,hist


def spectral_operator(fields,n,h,c,nu,lam,m):
    """Y_TT = B Y + G Y_T, Y=(xi,xi_bar,eta,eta_bar).

    N=0: all four components have angular index m. Complexification retains
    both real perturbation directions without restricting Z parity.
    """
    lap,dz,_=operators(n,h,m)
    sl=np.s_[0 if m==0 else 1:n,1:-1]
    u,w,s=[f[sl].ravel() for f in fields]
    p=u+1j*w; density=abs(p)**2
    diag=sp.diags
    a=lap+diag(1-2*density-4*s*s)
    b=lap+diag(3+nu*nu-4*density-2*lam*s*s)
    ps=diag(-4*p*s); pcs=diag(-4*p.conj()*s)
    bmat=sp.bmat([[a+1j*c*dz,diag(-p*p),ps,ps],
                 [diag(-p.conj()**2),a-1j*c*dz,pcs,pcs],
                 [pcs,ps,b,diag(-lam*s*s)],
                 [pcs,ps,diag(-lam*s*s),b]],format='csr')
    omega=np.sqrt(2); gv=c/(2*omega); gamma=np.sqrt(1+gv*gv)
    ident=sp.eye(lap.shape[0],format='csr')
    g=sp.block_diag([2*gv*dz-2j*omega*ident,2*gv*dz+2j*omega*ident,
                     2*gv*dz-2j*nu*gamma*ident,2*gv*dz+2j*nu*gamma*ident],format='csr')
    full=sp.bmat([[None,sp.eye(bmat.shape[0])],[bmat,g]],format='csc')
    # Physical longitudinal translation includes the carrier phase derivative.
    if m==0:
        pz=np.gradient(fields[0]+1j*fields[1],h,axis=1,edge_order=2)[sl].ravel()
        sz=np.gradient(fields[2],h,axis=1,edge_order=2)[sl].ravel()
        phase=np.r_[0*p,0*p,1j*s,-1j*s]
        translation=np.r_[gamma*pz,gamma*pz.conj(),gamma*sz-1j*nu*gv*s,gamma*sz+1j*nu*gv*s]
        modes={'carrier_phase':phase,'longitudinal_translation':translation}
    elif m==1:
        pr=np.gradient(fields[0]+1j*fields[1],h,axis=0,edge_order=2)[sl].ravel()
        sr=np.gradient(fields[2],h,axis=0,edge_order=2)[sl].ravel()
        modes={'transverse_translation':np.r_[pr,pr.conj(),sr,sr]}
    else: modes={}
    return full,bmat,g,modes,sl


def measure_null(b,y):
    return {'relative_l2':float(np.linalg.norm(b@y)/np.linalg.norm(y)),
            'relative_max':float(abs(b@y).max()/abs(y).max())}


def spectrum(full,b,g,modes,k,shift,tol):
    start=time.perf_counter()
    ev,vec=eigs(full,k=k,sigma=shift,which='LM',tol=tol,maxiter=6000,
                v0=np.random.default_rng(20260930).standard_normal(full.shape[0]))
    order=np.argsort(abs(ev)); ev=ev[order]; vec=vec[:,order]
    size=b.shape[0]; records=[]
    for value,v in zip(ev,vec.T):
        y=v[:size]; vel=v[size:]
        qres=b@y+value*(g@y)-value*value*y
        denom=(np.linalg.norm(b@y)+abs(value)*np.linalg.norm(g@y)+abs(value)**2*np.linalg.norm(y))
        records.append({'real':float(value.real),'imag':float(value.imag),'abs':float(abs(value)),
                        'first_order_residual':float(np.linalg.norm(full@v-value*v)/np.linalg.norm(v)),
                        'quadratic_absolute_residual':float(np.linalg.norm(qres)/np.linalg.norm(y)),
                        'quadratic_scaled_residual':float(np.linalg.norm(qres)/max(denom,1e-30)),
                        'velocity_consistency':float(np.linalg.norm(vel-value*y)/np.linalg.norm(y)),
                        'overlap':{name:float(abs(np.vdot(mode,y))/(np.linalg.norm(mode)*np.linalg.norm(y))) for name,mode in modes.items()}})
    return records,ev,vec,time.perf_counter()-start


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('profile',type=Path); ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--n',type=int); ap.add_argument('--h',type=float)
    ap.add_argument('--m',type=int,nargs='+',default=[0,1])
    ap.add_argument('--k',type=int,default=16);ap.add_argument('--shift',type=complex,default=0.0001)
    ap.add_argument('--tol',type=float,default=1e-10);ap.add_argument('--newton-steps',type=int,default=20)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    n,h,fields,c,nu,lam,winding=load_profile(a.profile,a.n,a.h)
    if n<4 or h<=0 or any(m<0 for m in a.m):
        raise ValueError('Require n>=4, h>0, m>=0')
    if winding!=0: raise ValueError('Only audited N=0 is implemented; shifted N±m sectors require separate axis spaces')
    fields,hist=stationary(fields,n,h,c,nu,lam,a.newton_steps)
    if not np.isfinite(hist[-1]) or hist[-1]>1e-9 or fields[2].max()<1e-4:
        raise RuntimeError(f'Nontrivial stationary profile not converged: residual={hist[-1]}')
    np.savez(a.output/'profile.npz',r=np.arange(n+1)*h,z=np.arange(-n,n+1)*h,u=fields[0],w=fields[1],s=fields[2],c=c,nu=nu,lambda_sigma=lam)
    report={'source':str(a.profile),'source_sha256':hashlib.sha256(a.profile.read_bytes()).hexdigest(),
            'n':n,'h':h,'L':n*h,'c':c,'nu':nu,'lambda_sigma':lam,'N':winding,
            'newton_history':hist,'stationary_residual':hist[-1],
            'carrier_max':float(fields[2].max()),'shift':[a.shift.real,a.shift.imag],
            'acceptance':'UNVALIDATED: translation grid/domain convergence and eigenfunction checks required; no instability claim',
            'sectors':{}}
    for m in a.m:
        full,b,g,modes,sl=spectral_operator(fields,n,h,c,nu,lam,m)
        nulls={name:measure_null(b,y) for name,y in modes.items()}
        print(json.dumps({'m':m,'matrix_shape':full.shape,'neutral_residuals':nulls}),flush=True)
        records,ev,vec,elapsed=spectrum(full,b,g,modes,a.k,a.shift,a.tol)
        np.savez_compressed(a.output/f'm{m}_eigenpairs.npz',eigenvalues=ev,eigenvectors=vec)
        report['sectors'][str(m)]={'dimension':full.shape[0],'neutral_residuals':nulls,'seconds':elapsed,'eigenvalues':records}
        (a.output/'spectrum.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(json.dumps({'m':m,'seconds':elapsed,'eigenvalues':[[x['real'],x['imag']] for x in records]}),flush=True)


if __name__=='__main__': main()
