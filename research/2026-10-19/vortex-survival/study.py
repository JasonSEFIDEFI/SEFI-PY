"""Finite-time axisymmetric survival audit of preserved moving-ring equations.

No new field laws; no extension of the production evolution guardrails.
"""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parent
CONSOLE=ROOT.parents[2]/'SEFI-PY Evolution/sim/vortex_console'
sys.path.insert(0,str(CONSOLE))
import physics

def cores(field,r,z):
    """Phase circulation on every cell, with an amplitude-validity check.

    Reports containing cells, not exact continuum root locations.
    """
    q=[field[:-1,:-1],field[1:,:-1],field[1:,1:],field[:-1,1:]]
    circulation=sum(np.angle(q[(i+1)%4]*q[i].conj()) for i in range(4))/(2*np.pi)
    minimum=np.minimum.reduce([abs(x) for x in q])
    # Zeros exactly on vertices need a separate warning; do not silently count.
    valid=(minimum>1e-8)
    ix=np.argwhere(valid & (abs(circulation)>.5))
    found=[dict(r=float((r[i]+r[i+1])/2),z=float((z[j]+z[j+1])/2),winding=int(round(circulation[i,j]))) for i,j in ix]
    return found,int(np.sum(~valid))

def diagnostics(e):
    d=e.diagnostics();p,s,pv,sv=e.y
    cells,ambiguous=cores(p,e.r,e.z)
    rr,zz=np.meshgrid(e.r,e.z,indexing='ij');rho=abs(s)**2;I=e.integ(rho)
    # Phi charge is excess over the fixed rotating boundary background.
    # Boundary conditions are not invariant under global phase rotations:
    # finite-box charge exchange must be measured, not assumed absent.
    qphi=e.integ((np.sqrt(2)*(abs(p)**2-1)+np.imag(p.conj()*pv))/e.gamma-e.v*np.imag(p.conj()*e.dz(p)))
    edge=(rr>.8*e.r[-1]) | (abs(zz)>.8*e.z[-1])
    d.update(core_cells=cells,ambiguous_cells=ambiguous,
             phi_excess_charge=qphi,carrier_norm=I,
             carrier_rms_radius=float(np.sqrt(e.integ((rr*rr+zz*zz)*rho)/I)),
             carrier_edge_fraction=float(e.integ(rho*edge)/I),
             max_velocity=float(max(abs(pv).max(),abs(sv).max())))
    return d

def run(name,reference,epsilon,dt,T=24.,kind='amplitude'):
    f,p=physics.load_profile(CONSOLE/f'data/charged_ring_{reference}.npz')
    p.update(epsilon=0.,dt=dt,duration=T,save_every=.5)
    e=physics.Evolution(f,p)
    rr,zz=np.meshgrid(e.r,e.z,indexing='ij')
    bump=np.exp(-((rr-e.radius)**2+zz**2)/2)
    bump[-1,:]=bump[:,0]=bump[:,-1]=0
    if kind=='amplitude':e.y[0]+=epsilon*bump;e.y[1]+=epsilon*bump
    elif kind=='phase':e.y[0]*=np.exp(1j*epsilon*bump);e.y[1]*=np.exp(-1j*epsilon*bump)
    e.q0=e.charge();initial=e.y.copy();samples=[diagnostics(e)]
    t0=time.monotonic();steps=round(T/dt);every=round(.5/dt);reason='completed'
    try:
        for i in range(steps):
            e.step(dt)
            if (i+1)%every==0:
                samples.append(diagnostics(e))
                if samples[-1]['max_velocity']>100:
                    reason='large_velocity_guard';break
    except (FloatingPointError,OverflowError) as exc:reason=type(exc).__name__
    initialcores=samples[0]['core_cells'];survival=all(len(d['core_cells'])>0 for d in samples)
    # Upper characteristic speed in the moving coordinates; rough warning,
    # not a rigorous causal-isolation proof (profiles have noncompact tails).
    boundary_warning=(p['L']-e.radius)/(e.gamma*(1+e.v))
    result=dict(name=name,reference=reference,epsilon=epsilon,kind=kind,dt=dt,T_requested=T,
        grid_intervals=len(e.r)-1,h=e.h,L=p['L'],reason=reason,runtime_seconds=time.monotonic()-t0,
        boundary_influence_warning_time=boundary_warning,core_detected_all_samples=survival,
        max_carrier_charge_drift=max(d['charge_drift'] for d in samples),
        phi_excess_charge_change=samples[-1]['phi_excess_charge']-samples[0]['phi_excess_charge'],
        samples=samples)
    np.savez_compressed(ROOT/(name+'.npz'),initial=initial,final=e.y,r=e.r,z=e.z)
    (ROOT/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='samples'}),flush=True)
    return result

if __name__=='__main__':
    cases=[('control','refined',0,.025,'amplitude'),
           ('small','refined',.001,.025,'amplitude'),
           ('medium','refined',.01,.025,'amplitude'),
           ('strong','refined',.02,.025,'amplitude'),
           ('phase','refined',.02,.025,'phase'),
           ('medium_half_dt','refined',.01,.0125,'amplitude'),
           ('coarse_medium','trial',.01,.025,'amplitude'),
           ('wide_medium','wide',.01,.025,'amplitude')]
    results=[run(a,b,c,d,kind=e) for a,b,c,d,e in cases]
    summary=[]
    for x in results:
        a=x['samples'][0];b=x['samples'][-1]
        summary.append({k:v for k,v in x.items() if k!='samples'}|dict(final_time=b['time'],initial_core_cells=a['core_cells'],final_core_cells=b['core_cells'],initial_rms=a['carrier_rms_radius'],final_rms=b['carrier_rms_radius'],final_edge_fraction=b['carrier_edge_fraction'],final_phi_distance=b['phi_distance']))
    data={'scope':'Finite-time axisymmetric numerical study; no proof of full stability or physical matter','cases':summary,'source_sha256':{n:hashlib.sha256((CONSOLE/n).read_bytes()).hexdigest() for n in ['physics.py','ring_operator.py']}}
    # Temporal refinement compares full complex fields on the same grid/time.
    a=np.load(ROOT/'medium.npz')['final'];b=np.load(ROOT/'medium_half_dt.npz')['final']
    data['temporal_refinement_max_field_difference']=float(abs(a-b).max())
    (ROOT/'summary.json').write_text(json.dumps(data,indent=2)+'\n')
