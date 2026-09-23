"""Exploratory axisymmetric GP profile solve. Not a stability solver."""
import argparse, json, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def setup(n,h,c):
    nu=n*n; nw=n*(n-1); size=nu+nw
    A=np.zeros((size,size)); b=np.zeros(size)
    def idx(i,j,imag):
        if i==n or j==n or (imag and j==0):return None
        return nu+i*(n-1)+j-1 if imag else i*n+j
    def put(row,i,j,imag,coef):
        k=idx(i,j,imag)
        if k is None:
            if not imag:b[row]+=coef
        else:A[row,k]+=coef
    for imag in [False,True]:
        for i in range(n):
            for j in range(1 if imag else 0,n):
                row=idx(i,j,imag)
                if i==0:
                    put(row,1,j,imag,4/h**2);put(row,0,j,imag,-4/h**2)
                else:
                    put(row,i-1,j,imag,(1-1/(2*i))/h**2)
                    put(row,i+1,j,imag,(1+1/(2*i))/h**2)
                    put(row,i,j,imag,-2/h**2)
                if j==0:
                    put(row,i,1,imag,2/h**2);put(row,i,0,imag,-2/h**2)
                else:
                    put(row,i,j-1,imag,1/h**2);put(row,i,j+1,imag,1/h**2)
                    put(row,i,j,imag,-2/h**2)
                if imag:
                    put(row,i,j+1,False,c/(2*h));put(row,i,j-1,False,-c/(2*h))
                elif j==0:put(row,i,1,True,-c/h)
                else:
                    put(row,i,j+1,True,-c/(2*h));put(row,i,j-1,True,c/(2*h))
    iu=np.array([i*n+j for i in range(n) for j in range(1,n)])
    iw=np.arange(nu,size)
    return A,b,iu,iw

def unpack(x,n):
    u=np.ones((n+1,n+1));w=np.zeros_like(u)
    u[:n,:n]=x[:n*n].reshape(n,n)
    w[:n,1:n]=x[n*n:].reshape(n,n-1)
    return u,w

def pack(u,w,n):return np.r_[u[:n,:n].ravel(),w[:n,1:n].ravel()]

def initial(n,h,radius,old=None):
    r=np.arange(n+1)*h;z=r.copy();rr,zz=np.meshgrid(r,z,indexing='ij')
    if old:
        d=np.load(old);ro=d['r'];zo=d['z']
        def interp(a):
            az=np.array([np.interp(z,zo,row) for row in a])
            return np.array([np.interp(r,ro,az[:,j]) for j in range(n+1)]).T
        u=interp(d['u']);w=interp(d['w'])
        # Outside a previous smaller box, extend to background rather than edge values.
        u[(rr>ro[-1])|(zz>zo[-1])]=1
        w[(rr>ro[-1])|(zz>zo[-1])]=0
    else:
        ph=np.arctan2(2*radius*zz,rr*rr+zz*zz-radius*radius)
        amp=np.tanh(np.sqrt((rr-radius)**2+zz**2)/np.sqrt(2))
        u=amp*np.cos(ph);w=amp*np.sin(ph)
    u[-1,:]=u[:,-1]=1;w[-1,:]=w[:,-1]=w[:,0]=0
    return pack(u,w,n)

def diagnostics(x,n,h,c):
    u,w=unpack(x,n);r=np.arange(n+1)*h;z=r.copy()
    ur,uz=np.gradient(u,h,h,edge_order=2);wr,wz=np.gradient(w,h,h,edge_order=2)
    ur[0,:]=wr[0,:]=uz[:,0]=0
    sq=u*u+w*w
    def integral(a):return float(4*np.pi*np.trapezoid(np.trapezoid(a*r[:,None],z,axis=1),r))
    Kperp=integral(ur*ur+wr*wr); Az=integral(uz*uz+wz*wz)
    V=integral((1-sq)**2)/4
    E=(Kperp+Az)/2+V
    p=-integral((u-1)*wz-w*uz)
    roots=[]
    for i in range(n):
        if u[i,0]*u[i+1,0]<0:roots.append(float(r[i]-u[i,0]*h/(u[i+1,0]-u[i,0])))
    return dict(c=c,n=n,h=h,domain_radius=float(r[-1]),domain_half_length=float(z[-1]),
        GP_energy=E,GP_renormalized_impulse=p,axial_gradient_integral=Az,
        transverse_gradient_integral=Kperp,potential_energy=V,
        core_radii=roots,center_real=float(u[0,0]),min_density=float(sq.min()),
        boundary_note='Finite cylinder, fixed psi=1 on outer radial and axial boundaries; real-even, imaginary-odd axial parity. These energies are finite-box estimates, not infinite-domain values.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=32);ap.add_argument('--h',type=float,default=.5)
    ap.add_argument('--c',type=float,default=.6);ap.add_argument('--radius',type=float,default=5.)
    ap.add_argument('--initial');ap.add_argument('--name',default='ring_trial');ap.add_argument('--steps',type=int,default=25)
    a=ap.parse_args();n=a.n
    A,b,iu,iw=setup(n,a.h,a.c);x=initial(n,a.h,a.radius,a.initial)
    diag=np.diag_indices_from(A)
    def residual(x):
        u,w=unpack(x,n);d=1-u*u-w*w
        return A@x+b+pack(d*u,d*w,n)
    hist=[]
    for it in range(a.steps):
        f=residual(x);err=float(np.max(np.abs(f)));norm=float(np.linalg.norm(f))
        print(json.dumps({'iteration':it,'max_residual':err,'norm':norm}),flush=True)
        hist.append({'iteration':it,'max_residual':err,'norm':norm})
        if err<1e-9:break
        u,w=unpack(x,n)
        J=A.copy();J[diag]+=np.r_[(1-3*u[:n,:n]**2-w[:n,:n]**2).ravel(),(1-u[:n,1:n]**2-3*w[:n,1:n]**2).ravel()]
        cross=(-2*u[:n,1:n]*w[:n,1:n]).ravel();J[iu,iw]+=cross;J[iw,iu]+=cross
        dx=np.linalg.solve(J,-f)
        step=1.
        while step>1/2048:
            new=x+step*dx
            if np.linalg.norm(residual(new)) < norm*(1-1e-4*step):break
            step/=2
        if step<=1/2048:
            print('Line search stalled',flush=True);break
        x=new
    d=diagnostics(x,n,a.h,a.c);d['max_residual']=float(np.max(np.abs(residual(x))));d['history']=hist
    d['converged']=d['max_residual']<1e-8
    u,w=unpack(x,n);p=ROOT/'work'/a.name
    np.savez(str(p)+'.npz',u=u,w=w,r=np.arange(n+1)*a.h,z=np.arange(n+1)*a.h,c=a.c)
    Path(str(p)+'.json').write_text(json.dumps(d,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k!='history'},indent=2),flush=True)

if __name__=='__main__':main()
