"""Exact homogeneous dispersion versus derivative-expanded phase cones."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvalsh


def potential_hessian(phi,sigma,lam=25):
    A=np.dot(phi,phi);B=np.dot(sigma,sigma)
    return np.block([[(1+A+4*B)*np.eye(2)+2*np.outer(phi,phi),8*np.outer(phi,sigma)],
        [8*np.outer(sigma,phi),(lam*B+4*A-3)*np.eye(2)+2*lam*np.outer(sigma,sigma)]])


def coupled_phase_speeds(A,B,lam=25):
    X=1+A+4*B;Y=-3+4*A+lam*B
    if min(A,B,X,Y)<=0 or lam<=16:raise ValueError('Requires positive two-condensate branch and positive amplitude susceptibility')
    omega,nu=np.sqrt(X),np.sqrt(Y)
    spatial=np.diag([A,B])
    temporal=spatial+2/(lam-16)*np.array([[lam*omega**2,-4*omega*nu],[-4*omega*nu,nu**2]])
    return eigvalsh(spatial,temporal),temporal,spatial


def exact_coupled_frequencies(k,A,B,lam=25):
    omega=np.sqrt(1+A+4*B);nu=np.sqrt(-3+4*A+lam*B)
    cross=-4*np.sqrt(A*B)
    b=np.array([[-k*k-A,-A,cross,cross],[-A,-k*k-A,cross,cross],
        [cross,cross,-k*k-lam*B,-lam*B],[cross,cross,-lam*B,-k*k-lam*B]],complex)
    g=np.diag([-2j*omega,2j*omega,-2j*nu,2j*nu])
    a=np.block([[np.zeros((4,4)),np.eye(4)],[b,g]])
    values=np.linalg.eigvals(a)
    if abs(values.real).max()>1e-7:raise RuntimeError('Homogeneous example has nonimaginary frequency; phase-speed comparison cannot assume stability')
    return np.sort(values.imag[values.imag>1e-8])


def main():
    out=Path(__file__).with_name('three_track_runs');out.mkdir(exist_ok=True)
    speeds,T,D=coupled_phase_speeds(1.,.1)
    rows=[]
    for k in (.001,.003,.01,.1,1.,10.,100.,1000.,10000.):
        # Rationalized gapless root avoids catastrophic cancellation at small k.
        root=np.sqrt(25+8*k*k)
        phase2=(k*k*(k*k+2))/(k*k+5+root)
        rows.append({'k':k,'single_phi_gapless_frequency':float(np.sqrt(phase2)),
            'single_phi_gapped_frequency':float(np.sqrt(k*k+5+root)),
            'carrier_lab_frequency':float(np.sqrt(k*k+1)),
            'two_condensate_exact_frequencies':exact_coupled_frequencies(k,1.,.1).tolist()})
    result={'microscopic_principal_symbol':'eta^{mu nu} k_mu k_nu times I_4 in Cartesian real fields',
        'microscopic_characteristic_determinant':'(eta^{mu nu} k_mu k_nu)^4',
        'microscopic_front_speed':1.,'homogeneous_asymptotic_phi_phase_sound_speed_squared':.2,
        'carrier_mass_squared_at_asymptotic_background':1.,
        'two_condensate_example':{'A':1.,'B':.1,'X':2.4,'Y':3.5,'lambda':25,
            'phase_speeds_squared':speeds.tolist(),'temporal_matrix':T.tolist(),'spatial_matrix':D.tolist(),
            'is_the_localized_branch_asymptotic_background':False},
        'dispersion':rows,
        'conclusion':'Canonical microscopic cones are universally the already-assumed Minkowski cone. The infrared Phi phase cone differs; carrier/amplitude dispersion and coupled two-phase cones do not share that infrared metric. No universal emergent geometry or gravitational field equation is established.'}
    (out/'characteristics.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'phase_speeds_squared':speeds.tolist(),'microscopic_front_speed':1.,'asymptotic_single_phase_speed_squared':.2}))


if __name__=='__main__':main()
