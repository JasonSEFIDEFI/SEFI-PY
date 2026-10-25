"""Symbolic local expansion and direct finite-delta functional checks."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.integrate import simpson
from derive_correction import evaluate,BETA_TURN,BETA_BIND
HERE=Path(__file__).resolve().parent

def main():
    out=HERE/'expansion_validation.json'
    if out.exists():raise FileExistsError(out)
    d,theta,p,chi,q,beta,lap=sp.symbols('d theta p chi q beta lap',real=True)
    v=(1+sp.sqrt(1+sp.Rational(16,3)*d))/8
    lam=1/v**2
    f=sp.sqrt(3)/2*sp.sin(theta);s0=sp.sqrt(3)*v*sp.cos(theta);s=s0+d*chi
    fp=sp.diff(f,theta)*p;sz=sp.diff(s0,theta)*p+d*q
    rho=f*f
    # Expand the squared stiff term separately to avoid cancellation at d=0.
    square=lam*(s*s-s0*s0)**2/(4*d)
    potential=((sp.Rational(1,4)-4*v*v)/d)*rho**2+beta*rho/2
    density=(fp*fp+sz*sz)/2+square+potential
    coefficient=sp.simplify(sp.diff(sp.series(density,d,0,2).removeO(),d).subs(d,0))
    expected=sp.sin(theta)**2*p*p/4+sp.Rational(4,9)*rho**2-sp.sqrt(3)/4*sp.sin(theta)*p*q+3*sp.cos(theta)**2*chi**2
    difference=sp.trigsimp(coefficient-expected)
    square_check=sp.simplify(3*sp.cos(theta)**2*chi**2-lap*chi+lap**2/(12*sp.cos(theta)**2)
        -3*sp.cos(theta)**2*(chi-lap/(6*sp.cos(theta)**2))**2)
    assert difference==0 and square_check==0,(difference,square_check)
    result=dict(scope=__doc__,symbolic_checks=dict(local_first_order=str(difference),
        transverse_completion_of_square=str(square_check)),functional_checks=[])
    for beta0 in [BETA_TURN,BETA_BIND]:
        base,profile=evaluate(beta0)
        z=profile['z'];t=profile['theta'];theta_z=np.sqrt(beta0)*profile['v'];chi0=profile['chi']
        chi_z=np.gradient(chi0,z,edge_order=2);chi_z[0]=0
        rho=.75*np.sin(t)**2
        fp=np.sqrt(3)/2*np.cos(t)*theta_z
        rows=[]
        for delta in [.01,.005,.0025,.00125]:
            vv=(1+np.sqrt(1+16*delta/3))/8
            lambda0=1/vv**2
            s0=np.sqrt(3)*vv*np.cos(t)
            ss=s0+delta*chi0
            sz=-np.sqrt(3)*vv*np.sin(t)*theta_z+delta*chi_z
            # Use the algebraically stable exact coefficient -(4+sqrt(lambda))/12.
            density=(fp*fp+sz*sz)/2+lambda0*(ss*ss-s0*s0)**2/(4*delta)
            density+=-(4+np.sqrt(lambda0))/12*rho*rho+beta0*rho/2
            actual=float(simpson(4*np.pi*z*z*density,x=z))
            residual=actual-base['A0']-delta*base['A1']
            rows.append(dict(delta=delta,lambda_sigma=lambda0,actual_action=actual,
                first_order_action=base['A0']+delta*base['A1'],
                remainder=residual,remainder_over_delta_squared=residual/delta**2))
        result['functional_checks'].append(dict(beta=beta0,rows=rows))
    predictions=json.loads((HERE/'predictions.json').read_text())
    env=[]
    for r in predictions['derivative_checks']:
        a,b,c,d=r['evaluations'];h=r['step']
        deriv=(a['A0']-8*b['A0']+8*c['A0']-d['A0'])/(12*h)
        err=abs(deriv-predictions['turn']['I0']/2)
        env.append(dict(step=h,derivative=deriv,expected=predictions['turn']['I0']/2,absolute_error=err))
    result['leading_envelope_checks']=env
    assert env[-1]['absolute_error']<1e-6,env
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(symbolic=result['symbolic_checks'],envelope=env,
        remainders=[(r['beta'],[x['remainder_over_delta_squared'] for x in r['rows']]) for r in result['functional_checks']])),flush=True)
if __name__=='__main__':main()

