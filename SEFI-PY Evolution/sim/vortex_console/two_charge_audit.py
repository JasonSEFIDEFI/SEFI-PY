"""Construct charge-matched initial velocities. No dynamics or stability claim."""
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
(ROOT/'runs').mkdir(exist_ok=True)
from rest_search import Search,load_initial

source=ROOT/'data/fixed_charge_Q1000_refined.npz'
records=[]
for support in [20.,28.,36.,48.,64.]:
    L=max(28,support+4);a=Search(L=L,h=.4,Q=1000,lam=12,N=1);load_initial(a,source)
    u,w,s=a.y;F=u+1j*w;rho=abs(F)**2;r,z=np.meshgrid(a.r,a.z,indexing='ij')
    window=np.maximum(1-(r*r+z*z)/support**2,0)**4
    B=a.integral(rho-1);q0=np.sqrt(2)*B;A=a.integral(window*rho);D=a.integral(window**2*rho)
    domega=-q0/A;pv=1j*domega*window*F
    q1=a.integral(np.sqrt(2)*(rho-1)+np.imag(F.conj()*pv));terms=a.terms(a.y)
    susceptibility=A*A/D;extra=q0*q0/(2*susceptibility);cost=a.integral(abs(pv)**2)/2
    assert abs(q1)<1e-8 and abs(extra-cost)<1e-8
    E_initial=terms['energy']+np.sqrt(2)*q0
    E_corrected=E_initial+a.integral(np.sqrt(2)*np.imag(F.conj()*pv)+abs(pv)**2/2)
    assert abs(E_corrected-terms['energy']-extra)<1e-8
    records.append(dict(support_radius=support,L=L,h=.4,background_charge_before=q0,background_charge_after=q1,carrier_charge=terms['nu']*terms['I'],phase_rate_peak=domega,susceptibility=susceptibility,compensation_energy_cost=extra,rotating_functional=terms['energy'],charge_matched_relative_energy=E_corrected,energy_per_carrier_charge=E_corrected/1000,below_dilute_carrier_threshold=bool(E_corrected<1000)))
report=dict(source=str(source.relative_to(ROOT)),ansatz='d_T F = i deltaOmega f F; d_T s=0; f=max(1-(r^2+z^2)/R^2,0)^4',formula='deltaOmega=-DeltaQ_Phi/A; susceptibility=A^2/D; E_matched=R_Q+(DeltaQ_Phi)^2/(2*susceptibility)',results=records,limitations=['Charge-matched initial data, not stationary solutions or physical evolution.','Spatial fields are interpolated from a finite-domain candidate and extended by the homogeneous boundary values.','Both U(1) charges are matched; axial momentum and full constrained stability are not analyzed.','Angular momentum about the symmetry axis is unchanged by the axisymmetric original-field phase-rate correction.','Compensation support depends on the surrounding medium; this does not establish intrinsic particle mass.'])
(ROOT/'runs/two_charge_torus_audit.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))


