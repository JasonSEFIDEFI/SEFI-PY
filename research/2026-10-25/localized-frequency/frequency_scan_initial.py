"""Frequency continuation of spherical two-complex-scalar stationary states.
Exploratory / Not Promoted. All quantities dimensionless.
"""
from pathlib import Path
import argparse, hashlib, json, platform, time
import numpy as np
import scipy
from scipy.integrate import solve_bvp,simpson
from scipy.optimize import minimize_scalar
HERE=Path(__file__).resolve().parent

def scales(lam,beta):
    if not 0<lam<16 or not 0<beta<1:
        raise ValueError('This subthreshold frequency chart requires 0<lambda<16 and 0<beta<1')
    m2=1+12/lam
    delta=12/lam-3/np.sqrt(lam)
    kappa=np.sqrt(beta*delta)
    return m2,delta,kappa,np.sqrt(m2-kappa*kappa)

def solve(lam,beta,seed=None,X=32.,tol=1e-7):
    t0=time.monotonic()
    m2,delta,kappa,w=scales(lam,beta);vac=np.sqrt(3/lam)
    x=np.linspace(0,X,801)
    if seed is None:
        data=np.load(HERE/'seed_Q1250.npz')
        old_kappa=np.sqrt(2-1.401223084350689**2)
        f=np.interp(x/old_kappa,data['r'],data['f'],right=0)
        s=np.interp(x/old_kappa,data['r'],data['s'],right=.5)
    else:
        old,sol=seed
        z=sol.sol(np.minimum(x,old['X']))
        f=z[0].copy();s=z[2].copy()*vac/np.sqrt(3/old['lambda_sigma'])
        # At the high-frequency weak-amplitude end, this scaling stays near the branch.
        if beta<old['beta']:
            ratio=np.sqrt(beta/old['beta'])
            f*=ratio
            s=vac+(s-vac)*ratio*ratio
        outside=x>old['X'];f[outside]=0;s[outside]=vac
    y=np.array([f,np.gradient(f,x),s,np.gradient(s,x)])
    y[1,0]=0;y[3,0]=0
    def fun(x,y):
        f,fx,s,sx=y
        return np.array([fx,(1+f*f+4*s*s-w*w)*f/kappa**2,
                         sx,(-3+4*f*f+lam*s*s)*s/kappa**2])
    def bc(a,b):return np.array([a[1],a[3],b[0],b[2]-vac])
    S=np.zeros((4,4));S[1,1]=-2;S[3,3]=-2
    sol=solve_bvp(fun,bc,x,y,S=S,tol=tol,bc_tol=1e-10,max_nodes=25000)
    xx=np.linspace(0,X,16001)
    f,fx,s,sx=sol.sol(xx)
    r=xx/kappa;fp=kappa*fx;sp=kappa*sx
    I=simpson(4*np.pi*r*r*f*f,x=r);Q=w*I
    U=simpson(4*np.pi*r*r*(f**4/4+f*f/2+2*f*f*s*s+lam*(s*s-3/lam)**2/4),x=r)
    K=simpson(4*np.pi*r*r*(fp*fp+sp*sp)/2,x=r)
    C=w*w*I/2;E=K+U+C
    maxf=float(np.max(abs(f)))
    nodes_in_f=int(np.count_nonzero(np.diff(np.sign(f[abs(f)>max(1e-8,maxf*1e-5)]))))
    tail=simpson(4*np.pi*r[-3201:]**2*f[-3201:]**2,x=r[-3201:])/max(I,1e-300)
    nontrivial=maxf>1e-3 and Q>1.
    accepted=sol.status==0 and nontrivial and nodes_in_f==0 and f[0]>0 and tail<1e-6
    result=dict(lambda_sigma=lam,beta=beta,omega=w,m_phi=float(np.sqrt(m2)),
        delta=delta,kappa=kappa,X=X,L=X/kappa,tolerance=tol,nodes=len(sol.x),
        status=int(sol.status),message=sol.message,accepted=bool(accepted),
        nontrivial=bool(nontrivial),radial_nodes=nodes_in_f,central_f=float(f[0]),central_s=float(s[0]),
        Q=float(Q),energy=float(E),energy_per_charge=float(E/Q) if Q>0 else None,
        binding_difference=float(E/Q-np.sqrt(m2)) if Q>0 else None,
        rms_radius=float(np.sqrt(simpson(4*np.pi*r**4*f*f,x=r)/I)) if I>0 else None,
        scaled_charge=float(Q*delta**1.5),scaled_radius=float(np.sqrt(delta)*np.sqrt(simpson(4*np.pi*r**4*f*f,x=r)/I)) if I>0 else None,
        max_rms_residual=float(max(sol.rms_residuals)),boundary_residual=float(max(abs(bc(sol.y[:,0],sol.y[:,-1])))),
        virial_relative=float((K+3*U-3*C)/E) if E else None,
        tail_charge_fraction=float(tail),wall_seconds=time.monotonic()-t0)
    return result,sol

def save_profile(path,result,sol):
    x=np.linspace(0,result['X'],1601)
    np.savez_compressed(path,x=x,y=sol.sol(x),lambda_sigma=result['lambda_sigma'],
                        beta=result['beta'],kappa=result['kappa'],omega=result['omega'])

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',required=True)
    p.add_argument('--lambdas',default='12,13,14,15,15.5,15.75,15.9')
    p.add_argument('--X',type=float,default=32)
    p.add_argument('--tol',type=float,default=1e-7)
    args=p.parse_args();out=HERE/args.out;out.mkdir(exist_ok=False)
    report=dict(scope=__doc__,python=platform.python_version(),numpy=np.__version__,
        scipy=scipy.__version__,configuration=vars(args),runs=[],folds=[])
    def persist():(out/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    anchor=None
    betas_low=[.275,.25,.225,.2,.175,.15,.125,.1,.075,.05,.025,.0125,.00625]
    betas_high=[.3,.35,.4,.45,.5,.6,.7,.8,.9]
    for lam in [float(x) for x in args.lambdas.split(',')]:
        cache=[]
        seed=anchor
        for beta in betas_low:
            result,sol=solve(lam,beta,seed,X=args.X,tol=args.tol)
            report['runs'].append(result)
            tag='lambda'+format(lam,'.6g')+'_beta'+format(beta,'.6g')
            save_profile(out/(tag+'.npz'),result,sol)
            persist();print(json.dumps(result),flush=True)
            if result['accepted']:
                seed=(result,sol);cache.append(seed)
                if beta==.275:anchor=seed
            else:break
        seed=anchor
        for beta in betas_high:
            result,sol=solve(lam,beta,seed,X=args.X,tol=args.tol)
            report['runs'].append(result)
            tag='lambda'+format(lam,'.6g')+'_beta'+format(beta,'.6g')
            save_profile(out/(tag+'.npz'),result,sol)
            persist();print(json.dumps(result),flush=True)
            if result['accepted']:seed=(result,sol);cache.append(seed)
            else:break
        if len(cache)<4:continue
        def objective(beta):
            near=min(cache,key=lambda z:abs(z[0]['beta']-beta))
            result,sol=solve(lam,beta,near,X=args.X,tol=args.tol)
            report.setdefault('fold_search_evaluations',[]).append(result);persist()
            if not result['accepted']:raise RuntimeError('Fold solve failed: '+str(result))
            cache.append((result,sol))
            return result['Q']
        opt=minimize_scalar(objective,bounds=(.025,.5),method='bounded',
                            options=dict(xatol=2e-6,maxiter=40))
        best=min(cache,key=lambda z:z[0]['Q'])
        fold=dict(lambda_sigma=lam,success=bool(opt.success),beta=float(opt.x),
            Q=float(opt.fun),profile=best[0],
            scope='Numerical minimum of Q(beta) in declared interval; not global existence threshold')
        report['folds'].append(fold);persist();print(json.dumps(fold),flush=True)
    report['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    persist()

if __name__=='__main__':main()

