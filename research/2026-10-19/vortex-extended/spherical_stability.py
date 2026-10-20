"""Fixed-charge second variation, separated into spherical harmonics.

Real-amplitude block includes the fixed-charge rank-one term only for l=0.
Imaginary Phi/Sigma blocks test phase directions. These are energy curvatures,
not physical mode frequencies or a nonlinear stability theorem.
"""
from vacuum_spherical import *
from scipy.sparse import diags,bmat
from scipy.sparse.linalg import LinearOperator,eigsh

def operators(model,f,s,l):
    V=model.V;root=np.sqrt(V);omega=model.Q/np.sum(V*f*f)
    kin=diags([-model.links/(root[:-1]*root[1:]),model.diag/V,
               -model.links/(root[:-1]*root[1:])],[-1,0,1],format='csr')
    kin=kin+diags(l*(l+1)/model.r**2)
    a=kin+diags(1+3*f*f+4*s*s-omega**2)
    d=kin+diags(-3+4*f*f+3*model.lam*s*s)
    cross=diags(8*f*s)
    base=bmat([[a,cross],[cross,d]],format='csr')
    u=np.r_[root*f,np.zeros(len(f))]
    coeff=4*omega**2/np.sum(V*f*f) if l==0 else 0.
    amplitude=LinearOperator(base.shape,matvec=lambda x:base@x+coeff*u*np.dot(u,x),dtype=float)
    return amplitude,kin+diags(1+f*f+4*s*s-omega**2),kin+diags(-3+4*f*f+model.lam*s*s)

def check_hessian():
    model=Radial(h=.4);f,s=model.initial();op,_,_=operators(model,f,s,0)
    rng=np.random.default_rng(92);q=rng.normal(size=2*len(f));q/=np.linalg.norm(q)
    root=np.sqrt(model.V);df,ds=np.split(q.copy(),2);df/=root;ds/=root
    def grad(t):
        gf,gs,_=model.gradient(f+t*df,s+t*ds)
        return -np.r_[root*gf,root*gs]
    actual=(grad(1e-4)-grad(-1e-4))/(2e-4)
    error=np.linalg.norm(actual-op@q)/np.linalg.norm(op@q)
    assert error<1e-7,error
    return float(error)

def run(tag):
    out=ROOT/(tag+'_curvatures.json')
    if out.exists():raise FileExistsError(out)
    status=ROOT/(tag+'_curvatures.status.json');status.write_text('{"reason":"running"}')
    meta=json.loads((ROOT/(tag+'.json')).read_text());data=np.load(ROOT/(tag+'.npz'))
    assert meta['converged']
    model=Radial(h=meta['h'],L=meta['L']);f=data['f'];s=data['s']
    result=dict(source=tag,scope=__doc__,hessian_directional_relative_error=check_hessian(),sectors=[])
    for l in range(5):
        if (ROOT/'STOP_REQUEST').exists():
            status.write_text('{"reason":"user_interrupted"}');return
        entry=dict(l=l)
        for name,op in zip(['real_amplitudes','phi_phase','sigma_phase'],operators(model,f,s,l)):
            vals,vecs=eigsh(op,k=3,which='SA',tol=1e-9,maxiter=20000,v0=np.random.default_rng(12).normal(size=op.shape[0]))
            entry[name]=dict(eigenvalues=vals.tolist(),max_residual=float(max(np.linalg.norm(op@vecs[:,i]-vals[i]*vecs[:,i]) for i in range(3))))
        result['sectors'].append(entry)
        print(tag,l,{k:v['eigenvalues'][0] for k,v in entry.items() if isinstance(v,dict)},flush=True)
    # Continuum dilation identity; its residual also measures discretization.
    I=np.sum(model.V*f*f);C=model.Q**2/(2*I)
    rho=f*f;eta=s*s
    U=np.sum(model.V*(rho*rho/4+rho/2+2*rho*eta+model.lam*(eta-3/model.lam)**2/4))
    E=model.energy(f,s);K=E-U-C
    result['dilation']=dict(K=float(K),U=float(U),charge_kinetic=float(C),relative_virial_residual=float((K+3*U-3*C)/E))
    result['reason']='completed';out.write_text(json.dumps(result,indent=2)+'\n')
    status.write_text('{"reason":"completed"}')

if __name__=='__main__':
    for tag in ['vacuum_spherical_h0p4','vacuum_spherical_h0p2','vacuum_spherical_h0p1_L60','vacuum_spherical_h0p1_L80']:run(tag)
