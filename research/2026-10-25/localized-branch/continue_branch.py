"""Fixed-charge continuation and restricted fragmentation comparisons.
All values dimensionless. No continuum/global/nonlinear stability claim.
"""
from pathlib import Path
import argparse, hashlib, json, platform, time
import numpy as np
import scipy
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import splu
from vacuum_spherical import Radial

HERE=Path(__file__).resolve().parent

def polish(model,f,s,tol):
    root=np.sqrt(model.V)
    n=len(f)
    for iteration in range(30):
        gf,gs,w=model.gradient(f,s)
        residual=max(np.max(abs(gf)),np.max(abs(gs)))
        if residual<tol: return f,s,iteration
        kin=diags([-model.links/(root[:-1]*root[1:]),model.diag/model.V,
                   -model.links/(root[:-1]*root[1:])],[-1,0,1],format='csc')
        A=bmat([[kin+diags(1+3*f*f+4*s*s-w*w),diags(8*f*s)],
                [diags(8*f*s),kin+diags(-3+4*f*f+3*model.lam*s*s)]],format='csc')
        u=np.r_[root*f,np.zeros(n)]
        c=4*w*w/np.sum(model.V*f*f)
        lu=splu(A)
        rhs=np.r_[root*gf,root*gs]
        z=lu.solve(rhs);v=lu.solve(u)
        denom=1+c*np.dot(u,v)
        if abs(denom)<1e-12: break
        delta=z-v*(c*np.dot(u,z)/denom)
        df,ds=np.split(delta,2);df=df/root;ds=ds/root
        step=1.
        for _ in range(22):
            fn=f+step*df;sn=s+step*ds
            gn,hn,_=model.gradient(fn,sn)
            if max(np.max(abs(gn)),np.max(abs(hn)))<residual:
                f,s=fn,sn;break
            step*=.5
        else: break
    return f,s,iteration

def solve(Q,lam,L,h,seed=None,tol=1e-9,maximum=4000):
    t0=time.monotonic()
    m=Radial(Q=Q,lam=lam,L=L,h=h)
    if seed is None:
        f,s=m.initial()
        scale=(Q/1e5)**(1/3)
        f=np.interp(m.r/scale,m.r,f,left=f[0],right=0)
        s=np.interp(m.r/scale,m.r,s,left=s[0],right=m.vac)
    else:
        old,fo,so=seed
        scale=(Q/old.Q)**(1/3)
        f=np.interp(m.r/scale,old.r,fo,left=fo[0],right=0)
        s=np.interp(m.r/scale,old.r,so,left=so[0],right=old.vac)*m.vac/old.vac
    E=m.energy(f,s)
    history=[];reason='iteration_limit'
    for it in range(maximum):
        gf,gs,w=m.gradient(f,s)
        res=float(max(abs(gf).max(),abs(gs).max()))
        if it%100==0: history.append(dict(iteration=it,energy=E,residual=res))
        if res<1e-6:reason='polish';break
        df=m.precondition(gf,1+3*f*f+4*s*s+w*w)
        ds=m.precondition(gs,1+abs(-3+4*f*f+3*lam*s*s))
        descent=np.sum(m.V*(gf*df+gs*ds));step=1.
        for _ in range(24):
            fn=f+step*df;sn=s+step*ds;En=m.energy(fn,sn)
            if En<=E-1e-4*step*descent:break
            step*=.5
        else:reason='line_search_stalled';break
        f,s,E=fn,sn,En
    f,s,nit=polish(m,f,s,tol)
    gf,gs,w=m.gradient(f,s);res=float(max(abs(gf).max(),abs(gs).max()))
    E=m.energy(f,s);I=np.sum(m.V*f*f);C=Q*Q/(2*I)
    rho=f*f;eta=s*s
    U=np.sum(m.V*(rho*rho/4+rho/2+2*rho*eta+lam*(eta-3/lam)**2/4))
    K=E-U-C
    rms=np.sqrt(np.sum(m.V*m.r*m.r*f*f)/I)
    tail=float(np.sum(m.V[m.r>.8*L]*f[m.r>.8*L]**2)/I)
    converged=res<tol
    result=dict(Q=Q,lambda_sigma=lam,L=L,h=h,energy=float(E),
        energy_per_charge=float(E/Q),omega=w,
        free_phi_mass=float(np.sqrt(1+12/lam)),
        rms_radius=float(rms),max_residual=res,converged=converged,
        classification=('below_free_threshold_candidate' if E<Q*np.sqrt(1+12/lam)
                        else 'not_below_free_threshold') if converged else 'unresolved',
        virial_relative=float((K+3*U-3*C)/E),tail_charge_fraction=tail,
        relaxation_iterations=it,newton_iterations=nit,
        relaxation_stop=reason,wall_seconds=time.monotonic()-t0,history=history)
    return m,f,s,result

def checks():
    rng=np.random.default_rng(42);results=[]
    for Q,lam in [(5e4,12.),(1e5,15.)]:
        m=Radial(Q=Q,lam=lam,L=60,h=.4);f,s=m.initial()
        df=.01*rng.normal(size=f.shape);ds=.01*rng.normal(size=s.shape)
        gf,gs,_=m.gradient(f,s);eps=1e-5
        fd=(m.energy(f+eps*df,s+eps*ds)-m.energy(f-eps*df,s-eps*ds))/(2*eps)
        an=-np.sum(m.V*(gf*df+gs*ds))
        err=abs(fd-an)/max(1,abs(an))
        assert err<2e-6,(err,fd,an)
        results.append(dict(Q=Q,lambda_sigma=lam,relative_gradient_error=float(err)))
    return results

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--h',type=float,default=.2)
    p.add_argument('--L',type=float,default=80.)
    p.add_argument('--lam',type=float,default=12.)
    p.add_argument('--charges',default='100000,80000,60000,50000,40000,30000,25000,20000,15000,10000,5000,2500,1250')
    p.add_argument('--out',required=True)
    args=p.parse_args()
    if args.h<=0 or args.L<=0 or abs(args.L/args.h-round(args.L/args.h))>1e-8:
        raise ValueError('Positive h,L with integer L/h required')
    charges=[float(x) for x in args.charges.split(',')]
    if args.lam<=0 or any(q<=0 for q in charges):raise ValueError('Positive lambda and charges required')
    out=HERE/args.out
    out.mkdir(parents=True,exist_ok=False)
    report=dict(scope=__doc__,source_commit='847a071351ed2b1c54e240fd9a520a1cce238178',
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        configuration=vars(args),gradient_checks=checks(),runs=[])
    seed=None
    for Q in charges:
        m,f,s,r=solve(Q,args.lam,args.L,args.h,seed)
        tag='Q'+format(Q,'.10g')
        np.savez_compressed(out/(tag+'.npz'),r=m.r,f=f,s=s)
        report['runs'].append(r)
        if r['converged']: seed=(m,f,s)
        (out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in r.items() if k!='history'}),flush=True)
    energies={r['Q']:r['energy'] for r in report['runs'] if r['converged']}
    report['fragmentation']=[]
    for Q,E in sorted(energies.items()):
        for q,e in sorted(energies.items()):
            if q<=Q/2 and Q-q in energies:
                margin=e+energies[Q-q]-E
                report['fragmentation'].append(dict(Q=Q,Q1=q,Q2=Q-q,
                    daughter_minus_parent_energy=margin,relative_margin=margin/E,
                    scope='Compared spherical stationary candidates only; not a global fission bound'))
    for name in ['continue_branch.py','vacuum_spherical.py']:
        report.setdefault('source_sha256',{})[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    (out/'results.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()

