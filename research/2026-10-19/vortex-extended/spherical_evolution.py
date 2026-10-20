"""Physical radial time evolution of the separate lambda=12 vacuum candidate.

Uses Phi=exp(i omega t)F; Sigma boundary is fixed to its real vacuum value.
RK4 evolves the equations, not the fictitious relaxation. Radial symmetry is
retained: these finite-time tests are not unrestricted nonlinear stability.
"""
from vacuum_spherical import *

def run(name,kind='amplitude',epsilon=.01,dt=.02,T=100.,wall=300):
    if (ROOT/(name+'.json')).exists():raise FileExistsError(name)
    status=ROOT/(name+'.status.json')
    def mark(reason,**kw):
        temp=status.with_suffix('.tmp');temp.write_text(json.dumps(dict(reason=reason,T_requested=T,dt=dt,kind=kind,epsilon=epsilon,**kw)));temp.replace(status)
    mark('running');start=time.monotonic()
    tag='vacuum_spherical_h0p2_L80';meta=json.loads((ROOT/(tag+'.json')).read_text())
    data=np.load(ROOT/(tag+'.npz'));model=Radial(h=meta['h'],L=meta['L']);omega=meta['omega']
    f=data['f'];s=data['s'];bump=np.exp(-((model.r-27)/5)**2)
    F=f*(1+epsilon*bump*(1j if kind=='phase' else 1));S=s.astype(complex)
    if kind=='phase':S+=1j*epsilon*bump*s
    else:S*=1-epsilon*bump
    # Set the physical Phi velocity to retain the requested initial charge.
    rotation=model.Q/np.sum(model.V*abs(F)**2)
    y=np.array([F,1j*(rotation-omega)*F,S,np.zeros_like(F)])
    def rhs(y):
        F,P,S,H=y;rho=abs(F)**2;eta=abs(S)**2
        return np.array([P,-model.K(F)/model.V-2j*omega*P+(omega**2-1-rho-4*eta)*F,
                         H,-(model.K(S)-model.b)/model.V+(3-4*rho-model.lam*eta)*S])
    def measure(y,t):
        F,P,S,H=y;rho=abs(F)**2;eta=abs(S)**2;physical=P+1j*omega*F
        kinetic=.5*np.sum(model.V*(abs(physical)**2+abs(H)**2))
        gradient=.5*np.sum(model.links*(abs(np.diff(F))**2+abs(np.diff(S))**2))+.5*model.outer*(abs(F[-1])**2+abs(S[-1]-model.vac)**2)
        potential=np.sum(model.V*(rho*rho/4+rho/2+2*rho*eta+model.lam*(eta-3/model.lam)**2/4))
        return dict(time=t,energy=float(kinetic+gradient+potential),charge=float(np.sum(model.V*np.imag(F.conj()*physical))),
            sigma_charge=float(np.sum(model.V*np.imag(S.conj()*H))),rms_phi_radius=float(np.sqrt(np.sum(model.V*model.r**2*rho)/np.sum(model.V*rho))),
            max_phi=float(abs(F).max()),amplitude_profile_relative_change=float(np.sqrt(np.sum(model.V*(abs(F)-f)**2)/np.sum(model.V*f*f))))
    n=round(T/dt);assert abs(n*dt-T)<1e-9
    samples=[measure(y,0)];reason='completed'
    for step in range(1,n+1):
        if (ROOT/'STOP_REQUEST').exists():reason='user_interrupted';break
        if time.monotonic()-start>wall:reason='wall_limit';break
        k1=rhs(y);k2=rhs(y+dt*k1/2);k3=rhs(y+dt*k2/2);k4=rhs(y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        if not np.isfinite(y).all():reason='nonfinite';break
        if step%max(1,n//200)==0:samples.append(measure(y,step*dt))
    achieved=(step if reason=='completed' else step-1)*dt
    if samples[-1]['time']!=achieved:samples.append(measure(y,achieved))
    result=dict(name=name,source=tag,scope=__doc__,T_requested=T,dt=dt,kind=kind,epsilon=epsilon,reason=reason,
        wall_seconds=time.monotonic()-start,samples=samples,
        max_relative_energy_drift=max(abs(v['energy']/samples[0]['energy']-1) for v in samples),
        max_relative_charge_drift=max(abs(v['charge']/samples[0]['charge']-1) for v in samples))
    (ROOT/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(ROOT/(name+'.npz'),state=y,r=model.r)
    mark(reason,achieved_time=achieved,eligible_completed_run=reason=='completed' and achieved==T)
    print(json.dumps({k:v for k,v in result.items() if k not in ['samples','scope']}),flush=True)

if __name__=='__main__':
    run('radial_control_T100',epsilon=0)
    run('radial_amplitude_T100')
    run('radial_amplitude_halfdt_T100',dt=.01)
    run('radial_phase_T100',kind='phase')
