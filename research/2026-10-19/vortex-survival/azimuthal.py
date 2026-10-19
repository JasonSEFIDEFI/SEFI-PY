"""Linear non-axisymmetric cos(m theta) perturbations of an N=0 ring.

This is a tangent evolution, not full nonlinear 3D evolution.
"""
from study import *

def tangent(e,y,m):
    if m<1 or int(m)!=m:raise ValueError('Positive integer angular mode required')
    p,s=e.base,e.carrier;dp,ds,vp,vs=y
    def lap(q):
        out=e.lap(q)
        out[1:-1,1:-1]-=m*m*q[1:-1,1:-1]/e.r[1:-1,None]**2
        out[0,:]=0
        return out
    drho=2*np.real(p.conj()*dp);deta=2*np.real(s.conj()*ds)
    pa=2*e.gamma*e.v*e.dz(vp)-2j*np.sqrt(2)*vp+lap(dp)+1j*e.p['c']*e.dz(dp)+(1-abs(p)**2-4*abs(s)**2)*dp-(drho+4*deta)*p
    sa=2*e.gamma*e.v*e.dz(vs)-2j*e.p['nu']*e.gamma*vs+lap(ds)+(3+e.p['nu']**2-4*abs(p)**2-25*abs(s)**2)*ds-(4*drho+25*deta)*s
    for a in [pa,sa]:a[0,:]=a[-1,:]=a[:,0]=a[:,-1]=0
    return np.array([vp,vs,pa,sa])

def step(e,y,m,dt):
    a=tangent(e,y,m);b=tangent(e,y+dt*a/2,m);c=tangent(e,y+dt*b/2,m);d=tangent(e,y+dt*c,m)
    return y+dt*(a+2*b+2*c+d)/6

def run_mode(m,dt=.025,T=24.,reference='refined'):
    f,p=physics.load_profile(CONSOLE/f'data/charged_ring_{reference}.npz');p.update(epsilon=0,dt=dt)
    e=physics.Evolution(f,p)
    rr,zz=np.meshgrid(e.r,e.z,indexing='ij')
    # r^m factor imposes regularity on the axis; support is smooth near the ring.
    g=(rr/e.radius)**m*np.exp(-((rr-e.radius)**2+zz**2)/2)
    g[0,:]=g[-1,:]=g[:,0]=g[:,-1]=0
    y=np.array([g.astype(complex),(.4+.3j)*g,np.zeros_like(g),np.zeros_like(g)])
    def norm(y):return float(np.sqrt(sum(e.integ(abs(a)**2) for a in y)))
    n0=norm(y);y/=n0;samples=[dict(time=0.,norm=1.)]
    for i in range(round(T/dt)):
        y=step(e,y,m,dt)
        if (i+1)%round(.5/dt)==0:samples.append(dict(time=(i+1)*dt,norm=norm(y)))
    result=dict(m=m,reference=reference,dt=dt,T=T,scope='Linear cos(m theta) disturbance; selected initial direction, not full mode spectrum',final_gain=samples[-1]['norm'],maximum_gain=max(x['norm'] for x in samples),samples=samples)
    tag=f'angular_m{m}'+('_half_dt' if dt<.025 else '')+('' if reference=='refined' else '_'+reference)
    (ROOT/(tag+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(ROOT/(tag+'.npz'),final=y)
    print(json.dumps({k:v for k,v in result.items() if k!='samples'}),flush=True)
    return result

if __name__=='__main__':
    results=[run_mode(m) for m in [1,2,3,4]]+[run_mode(2,.0125)]
    a=np.load(ROOT/'angular_m2.npz')['final'];b=np.load(ROOT/'angular_m2_half_dt.npz')['final']
    (ROOT/'angular_summary.json').write_text(json.dumps({'cases':[{k:v for k,v in x.items() if k!='samples'} for x in results],'m2_temporal_refinement_max_difference':float(abs(a-b).max())},indent=2)+'\n')
