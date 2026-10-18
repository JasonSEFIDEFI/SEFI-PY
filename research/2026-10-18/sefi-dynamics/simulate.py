"""1D candidate scalar-field recoil simulator. Not a validated SEFI engine.

Dimensionless c=1; nonrelativistic rigid source carrier; Gaussian sources;
Fourier Galerkin spatial representation, velocity Verlet time integration.
No external momentum sink. A prescribed source modulation supplies work.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import json, math, time
import numpy as np

ROOT=Path(__file__).resolve().parent

@dataclass
class Parameters:
    length: float=160.
    grid: int=512
    dt: float=.02
    duration: float=52.
    pulse_duration: float=24.
    mass: float=50.
    Z: float=1.
    mu: float=.2
    sigma: float=.65
    amplitude: float=.3
    omega: float=1.
    separation: float=math.pi/(2*math.sqrt(1-.2**2))

CASES={'off':('Sources off',0.),
       'symmetric':('Symmetric pulse',0.),
       'phased':('Phased pulse',math.pi/2),
       'reversed':('Reversed phase',-math.pi/2)}

def drive(t,p,phase,off=False):
    if off or not 0<t<p.pulse_duration:
        return np.zeros(2),np.zeros(2)
    envelope=math.sin(math.pi*t/p.pulse_duration)**2
    rate=math.pi/p.pulse_duration*math.sin(2*math.pi*t/p.pulse_duration)
    angles=p.omega*t+np.array([phase/2,-phase/2])
    return (p.amplitude*envelope*np.cos(angles),
            p.amplitude*(rate*np.cos(angles)-envelope*p.omega*np.sin(angles)))

def run(p,case,keep_frames=False):
    _,phase=CASES[case]
    k=2*math.pi*np.arange(p.grid//2)/p.length
    omega2=k*k+p.mu*p.mu
    weights=np.full(k.size,2.); weights[0]=1.
    smooth=np.exp(-.5*(p.sigma*k)**2)
    offsets=np.array([-p.separation/2,p.separation/2])
    w=np.zeros(k.size,dtype=complex); v=w.copy()
    position=0.; momentum=0.; work=0.; peak_speed=0.; peak_interaction=0.
    records=[]; frames=[]
    nsteps=round(p.duration/p.dt)
    assert abs(nsteps*p.dt-p.duration)<1e-9
    assert p.dt*math.sqrt(omega2[-1])<1.8, 'Time step too large'
    stride=max(1,round(.2/p.dt))

    def forces(t,x,w):
        g,dg=drive(t,p,phase,case=='off')
        phase_matrix=np.exp(1j*np.outer(offsets+x,k))
        avg=(phase_matrix@(weights*smooth*w)).real
        grad=(phase_matrix@(weights*smooth*1j*k*w)).real
        j=(g@phase_matrix.conjugate())*smooth/p.length
        return -omega2*w-j/p.Z, float(-g@grad), float(g@avg),float(dg@avg)

    def diagnostics(t,potential,power):
        nonlocal peak_speed,peak_interaction
        energy=.5*p.Z*p.length*float(np.sum(weights*(abs(v)**2+omega2*abs(w)**2)))
        field_p=-p.Z*p.length*float(np.sum(weights*(v.conjugate()*1j*k*w).real))
        kinetic=momentum**2/(2*p.mass)
        peak_speed=max(peak_speed,abs(momentum/p.mass))
        peak_interaction=max(peak_interaction,abs(potential)/p.mass)
        return {'t':t,'x':position,'v':momentum/p.mass,'vehicle_p':momentum,
            'field_p':field_p,'total_p':momentum+field_p,'field_energy':energy,
            'kinetic':kinetic,'interaction':potential,'work':work,'drive_power':power,
            'energy_residual':energy+kinetic+potential-work}

    acceleration,force,potential,power=forces(0,position,w)
    records.append(diagnostics(0,potential,power))
    for step in range(nsteps):
        v+=.5*p.dt*acceleration; momentum+=.5*p.dt*force
        w+=p.dt*v; position+=p.dt*momentum/p.mass
        t=(step+1)*p.dt
        acceleration,force,potential,next_power=forces(t,position,w)
        v+=.5*p.dt*acceleration; momentum+=.5*p.dt*force
        work+=.5*p.dt*(power+next_power); power=next_power
        if (step+1)%stride==0 or step+1==nsteps:
            records.append(diagnostics(t,potential,power))
            if keep_frames:
                # 256 displayed samples reconstructed from the retained modes.
                xs=np.linspace(-p.length/2,p.length/2,256,endpoint=False)
                coeff=w*weights
                field=(np.exp(1j*np.outer(xs,k))@coeff).real
                frames.append({'t':round(t,5),'x':round(position,9),
                               'field':np.round(field,6).tolist()})
    final=records[-1]
    scale=max(final['field_energy']+final['kinetic'],1e-30)
    summary={'name':CASES[case][0],'parameters':asdict(p),'final':final,
       'max_abs_total_momentum':max(abs(r['total_p']) for r in records),
       'max_abs_energy_residual':max(abs(r['energy_residual']) for r in records),
       'final_relative_energy_residual':abs(final['energy_residual'])/scale,
       'radiation_momentum_per_energy':abs(final['field_p'])/max(final['field_energy'],1e-30),
       'peak_sampled_speed_over_c':peak_speed,
       'peak_sampled_interaction_over_rest_energy':peak_interaction,
       'light_crossing_half_domain':p.length/2,
       'before_boundary_return':p.duration+p.separation/2+8*p.sigma<p.length/2}
    return {'summary':summary,'history':records,'frames':frames}

def main():
    started=time.time()
    p=Parameters(); runs={key:run(p,key,True) for key in CASES}
    convergence=[]
    for dt,grid in [(.01,512),(.005,512),(.005,1024)]:
        fine=Parameters(dt=dt,grid=grid)
        convergence.append(run(fine,'phased')['summary'])
    convergence.append(run(Parameters(dt=.005,grid=1024,length=320.),'phased')['summary'])
    a=runs['phased']['summary']; b=runs['reversed']['summary']
    tests={
      'off_is_stationary':abs(runs['off']['summary']['final']['vehicle_p'])<1e-14,
      'symmetric_has_no_net_recoil':abs(runs['symmetric']['summary']['final']['vehicle_p'])<1e-12,
      'phase_reversal_reverses_momentum':abs(a['final']['vehicle_p']+b['final']['vehicle_p'])<1e-12,
      'all_runs_conserve_total_momentum':all(r['summary']['max_abs_total_momentum']<1e-10 for r in runs.values()),
      'positive_final_outgoing_energy':a['final']['field_energy']>0,
      'radiation_momentum_bound':a['radiation_momentum_per_energy']<=1+1e-10,
      'weak_nonrelativistic_regime':a['peak_sampled_speed_over_c']<.01 and a['peak_sampled_interaction_over_rest_energy']<.01,
      'no_boundary_return_in_run':a['before_boundary_return'],
      'energy_error_falls_under_time_refinement':convergence[1]['final_relative_energy_residual']<convergence[0]['final_relative_energy_residual']/3,
      'spatial_refinement_agrees':abs(convergence[2]['final']['vehicle_p']-convergence[1]['final']['vehicle_p'])<1e-9,
      'doubled_domain_agrees':abs(convergence[3]['final']['vehicle_p']-convergence[1]['final']['vehicle_p'])<1e-8,
    }
    assert all(tests.values()),tests
    results={'model_status':'Additional hypothetical scalar coupling; not validated SEFI physics. Dimensionless 1D dynamics, not flight.',
       'tests':tests,'runs':{key:r['summary'] for key,r in runs.items()},
       'convergence':convergence,'elapsed_seconds':time.time()-started}
    (ROOT/'results.json').write_text(json.dumps(results,indent=2))
    # Downsample histories for inline plots. Numerical checks above use all recorded data.
    data={key:{'summary':r['summary'],'history':r['history'][::2],
               'frames':r['frames'][::4]+[r['frames'][-1]]} for key,r in runs.items()}
    (ROOT/'playback-data.json').write_text(json.dumps(data,separators=(',',':')))
    print(json.dumps(results,indent=2))

if __name__=='__main__': main()
