"""Independent discrete energy and original-field charge balance audit."""
from study import *

def energy(e,y):
    p,s,pv,sv=y;d=p-1;rho=abs(p)**2;eta=abs(s)**2
    kinetic=.5*e.integ(abs(pv)**2+abs(sv)**2)
    gradient=-.5*e.integ(np.real(d.conj()*e.lap(d)+s.conj()*e.lap(s)))
    potential=e.integ((rho-1)**2/4+25*eta**2/4+(4*rho-3-e.p['nu']**2)*eta/2)
    current=.5*e.p['c']*e.integ(np.imag(d.conj()*e.dz(d)))
    return kinetic+gradient+potential+current

def charge(e,y):
    p,_,pv,_=y
    return e.integ((np.sqrt(2)*(abs(p)**2-1)+np.imag(p.conj()*pv))/e.gamma-e.v*np.imag(p.conj()*e.dz(p)))

def charge_rate(e,y):
    p,_,pv,_=y
    # Derived after cancellation of real constitutive forces and rotation.
    return (e.integ(np.imag(p.conj()*e.lap(p)))/e.gamma
        +2*np.sqrt(2)*e.v*e.integ(np.real(p.conj()*e.dz(p)))
        +e.v*e.integ(np.imag(p.conj()*e.dz(pv)-pv.conj()*e.dz(p))))

def audit(reference='refined',dt=.025,T=24.):
    f,p=physics.load_profile(CONSOLE/f'data/charged_ring_{reference}.npz');p.update(epsilon=.01,dt=dt)
    e=physics.Evolution(f,p);Q0=charge(e,e.y);E0=energy(e,e.y);integrated=0.;samples=[]
    for i in range(round(T/dt)+1):
        if i%round(.5/dt)==0:
            samples.append(dict(time=e.time,energy=energy(e,e.y),phi_charge=charge(e,e.y),integrated_charge_rate=integrated))
        if i==round(T/dt):break
        y=e.y;a=e.rhs(y);yb=y+dt*a/2;b=e.rhs(yb);yc=y+dt*b/2;c=e.rhs(yc);yd=y+dt*c;d=e.rhs(yd)
        integrated+=dt*(charge_rate(e,y)+2*charge_rate(e,yb)+2*charge_rate(e,yc)+charge_rate(e,yd))/6
        e.y=y+dt*(a+2*b+2*c+d)/6;e.time+=dt
    result=dict(reference=reference,dt=dt,T=T,initial_canonical_energy=E0,
        max_relative_energy_drift=max(abs(x['energy']-E0) for x in samples)/max(abs(E0),1e-12),
        phi_charge_change=charge(e,e.y)-Q0,integrated_discrete_boundary_balance=integrated,
        charge_balance_residual=charge(e,e.y)-Q0-integrated,
        interpretation='Conserved rotating/moving-frame discrete energy is not isolated-particle rest energy. Fixed Phi boundary permits charge exchange; integrated finite-box balance is checked.',samples=samples)
    (ROOT/f'balance_{reference}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='samples'}),flush=True)
    return result

if __name__=='__main__':
    for reference in ['refined','wide']:audit(reference)
