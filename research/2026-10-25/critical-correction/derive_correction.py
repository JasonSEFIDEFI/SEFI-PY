"""Derive charge correction coefficients from the leading variational problem.
No full-model charge data enter these coefficients.
"""
from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
import numpy as np
from scipy.integrate import simpson
from reduced_solver import solve_limit
HERE=Path(__file__).resolve().parent
MC=np.sqrt(1.75)
BETA_TURN=.12437038000563955
BETA_BIND=.26513102811715317

def evaluate(beta,X=48.,tol=1e-10):
    sol=solve_limit(beta,X=X,tol=tol)
    x=np.linspace(0,X,24001);theta,v=sol.sol(x);z=x/np.sqrt(beta)
    sn=np.sin(theta);cs=np.cos(theta);rho=.75*sn*sn
    a=3/16*(1+3*cs*cs);ap=-9/8*sn*cs
    up=3/4*sn*cs*(beta-2*sn*sn)
    grad2=beta*v*v
    lap_theta=up/a-ap*grad2/(2*a)
    lap_s0=-np.sqrt(3)/4*(sn*lap_theta+cs*grad2)
    chi=lap_s0/(6*cs*cs)
    geometric=.25*sn*sn*grad2
    potential=4/9*rho*rho
    transverse=-lap_s0*lap_s0/(12*cs*cs)
    measure=4*np.pi*z*z
    integrate=lambda density:float(simpson(measure*density,x=z))
    I0=integrate(rho);K=integrate(a*grad2/2)
    A0=K+integrate(beta*rho/2-2*rho*rho/3)
    A1=integrate(geometric+potential+transverse)
    Q0=MC*I0
    result=dict(beta=beta,X=X,tolerance=tol,nodes=len(sol.x),
        max_rms_residual=float(max(sol.rms_residuals)),min_cos_theta=float(min(cs)),
        I0=I0,A0=A0,A1=A1,Q0=Q0,
        A1_terms=dict(geometric=integrate(geometric),potential=integrate(potential),transverse=integrate(transverse)),
        max_abs_chi=float(max(abs(chi))),central_chi=float(chi[0]))
    return result,dict(x=x,z=z,theta=theta,v=v,chi=chi,lap_s0=lap_s0)

def main():
    out=HERE/'predictions.json'
    if out.exists():raise FileExistsError(out)
    turn,profile=evaluate(BETA_TURN);binding,_=evaluate(BETA_BIND)
    derivative_checks=[]
    for h in [.001,.0005,.00025,.000125]:
        vals=[]
        for k in [-2,-1,1,2]:vals.append(evaluate(BETA_TURN+k*h)[0])
        a,b,c,d=vals
        derivative=(a['A1']-8*b['A1']+8*c['A1']-d['A1'])/(12*h)
        C1=2*MC*derivative+(2-BETA_TURN)*turn['I0']/(2*MC)
        derivative_checks.append(dict(step=h,A1_beta_derivative=derivative,turn_C1=C1,evaluations=vals))
    # Binding root displacement cancels Q1 in the charge at the displaced crossing.
    C1_binding=2*MC*binding['A1']/BETA_BIND+binding['Q0']*(4-BETA_BIND)/(4*MC*MC)
    result=dict(scope=__doc__,status='Exploratory / Not Promoted',
        created_utc=datetime.now(timezone.utc).isoformat(),
        turn=turn,binding=binding,derivative_checks=derivative_checks,
        coefficients=dict(turn_C0=turn['Q0'],turn_C1=derivative_checks[-1]['turn_C1'],
                          binding_C0=binding['Q0'],binding_C1=C1_binding),
        declared_heldout_couplings=[15.95,15.975],
        heldout_policy='Coefficients use reduced equations only; held-out full-model charges are computed afterward without fitting.',
        limitations='Formal regular expansion on a chart with cos(theta)>0; no uniform singular-limit theorem.')
    result['source_sha256']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['derive_correction.py','reduced_solver.py']}
    out.write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(HERE/'correction_profile.npz',**profile)
    print(json.dumps(dict(coefficients=result['coefficients'],A1_turn=turn['A1'],
        A1_binding=binding['A1'],derivative_checks=[(r['step'],r['turn_C1']) for r in derivative_checks])),flush=True)
if __name__=='__main__':main()

