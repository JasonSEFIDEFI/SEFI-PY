"""Check physical-potential reconstruction of the existing rotating equations."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def U0(rho,eta,lam=25):
    return (rho-1)**2/4+lam*eta**2/4+(4*rho-3)*eta/2

def Uphysical(rho,eta,lam=25):return U0(rho,eta,lam)+rho

def force(p,s,lam=25):
    rho=abs(p)**2;eta=abs(s)**2
    return np.array([-(1+rho+4*eta)*p,(3-4*rho-lam*eta)*s])

if __name__=='__main__':
    # The physical field is Phi=exp(i sqrt(2)t)F.
    # Phi_tt=exp(i sqrt(2)t)(F_tt+2i sqrt(2)F_t-2F).
    p=.8+.2j;s=.13-.07j;h=1e-6
    values=np.array([p.real,p.imag,s.real,s.imag])
    def energy(x):return Uphysical(x[0]**2+x[1]**2,x[2]**2+x[3]**2)
    grad=[]
    for i in range(4):
        d=np.zeros(4);d[i]=h;grad.append((energy(values+d)-energy(values-d))/(2*h))
    f=force(p,s);np.testing.assert_allclose(-np.array(grad),[f[0].real,f[0].imag,f[1].real,f[1].imag],atol=1e-9)
    rho=np.linspace(0,4,101);eta=np.linspace(0,1,101);x,y=np.meshgrid(rho,eta)
    # Exact identity proves the lower bound for rho,eta >= 0 at lambda=25.
    # Uphysical - .16 = rho^2/4+rho/2+2rho eta+25(eta-.12)^2/4.
    residual=Uphysical(x,y)-.16-(x*x/4+x/2+2*x*y+25*(y-.12)**2/4)
    assert np.max(abs(residual))<1e-13
    # At the true physical vacuum the Phi small-amplitude mass is sqrt(1.48).
    # A stronger exact identity excludes below-threshold Phi-charge binding
    # for this lambda=25 potential, without asserting nonexistence of all
    # possible metastable or medium-supported structures.
    binding_residual=(Uphysical(x,y)-.16-.74*x
        -(25/4*(y-.12+4*x/25)**2+.09*x*x))
    assert abs(binding_residual).max()<1e-13
    out=dict(scope='Potential audit of existing dimensionless scalar model; not a new equation set',
        omega_squared=2,lambda_sigma=25,
        physical_ground_state=dict(phi_squared=0,sigma_squared=.12,potential=.16,phase_degeneracy='Sigma phase arbitrary'),
        rotating_background=dict(phi_squared=1,sigma_squared=0,phi_charge_density=float(np.sqrt(2)),potential=1,kinetic_energy_density=1,total_energy_density=2),
        interpretation='The console ring is an excitation of a charged rotating background, not yet an isolated particle over the zero-charge physical ground state.',
        vacuum_phi_mass_squared=1.48,
        vacuum_binding_audit=dict(lambda25='Uphys-Umin >= 1.48 |Phi|^2/2, hence E_relative >= sqrt(1.48) |Q_Phi| in the zero-charge vacuum sector',
            lambda12_comparison=dict(vacuum_phi_mass_squared=2.,pointwise_charge_binding_threshold_squared=float(1+3/np.sqrt(12)),scope='Energetically permitted by the pointwise bound, not an existence or stability result; not the current console parameter set')),
        checks=dict(force_gradient_max_difference=float(np.max(abs(-np.array(grad)-np.array([f[0].real,f[0].imag,f[1].real,f[1].imag])))),potential_identity_max_residual=float(abs(residual).max()),binding_identity_max_residual=float(abs(binding_residual).max())))
    (ROOT/'physical_potential_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
