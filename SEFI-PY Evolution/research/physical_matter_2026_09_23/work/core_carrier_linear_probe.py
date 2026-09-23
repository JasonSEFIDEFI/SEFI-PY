"""Provisional two-field extension: test the onset of a core condensate only."""
import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
data=np.genfromtxt(root/'outputs'/'vortex_profile.csv',delimiter=',',names=True)
kappa=4.;mu2=3.;lam=25.;beta=1.;Rinf2=1.
results=[]
for L,n in [(20.,400),(20.,800),(30.,600)]:
    h=L/n;r=(np.arange(n)+.5)*h
    f=np.interp(r,data['radius'],data['amplitude'])
    potential=kappa*f*f-mu2
    left=np.arange(n)*h;right=(np.arange(n)+1)*h
    diag=(left+right)/(r*h*h)+potential
    diag[-1]+=right[-1]/(r[-1]*h*h) # outer value zero half a cell away
    off=-right[:-1]/(h*h*np.sqrt(r[:-1]*r[1:]))
    S=np.diag(diag)+np.diag(off,1)+np.diag(off,-1)
    values,vectors=np.linalg.eigh(S);mode=vectors[:,0]/np.sqrt(r*h)
    if mode[0]<0:mode=-mode
    probability=r*h*mode*mode;probability/=probability.sum()
    result={'L':L,'n':n,'h':h,'lowest_eigenvalues':values[:4].tolist(),
      'ground_mode_rms_radius':float(np.sqrt(np.sum(probability*r*r))),
      'probability_outside_radius_8':float(probability[r>8].sum()),
      'symmetric_eigen_residual':float(np.linalg.norm(S@vectors[:,0]-values[0]*vectors[:,0]))}
    results.append(result)
    if n==800:
        np.savetxt(root/'outputs'/'core_carrier_linear_mode.csv',np.c_[r,f,potential,mode],delimiter=',',header='radius,original_vortex_amplitude,carrier_potential,linear_mode',comments='')
report={'status':'New explicit model assumption; linear instability of the zero-carrier vortex toward a core condensate. Not a solved nonlinear condensate or stable loop.',
 'potential':'V=alpha/2 |Phi|^2+beta/4 |Phi|^4+lambda/4 |Sigma|^4+(kappa |Phi|^2-mu2)/2 |Sigma|^2',
 'parameters':{'alpha':1,'beta':beta,'Omega':float(np.sqrt(2)),'R_infinity_squared':Rinf2,'kappa':kappa,'mu2':mu2,'lambda':lam},
 'bulk_carrier_mass_squared':kappa*Rinf2-mu2,
 'core_carrier_mass_squared':-mu2,
 'rotating_potential_hessian_determinant_in_squared_amplitudes_up_to_positive_factor':beta*lam-kappa*kappa,
 'results':results,
 'limitations':['Second complex field and cross coupling are added assumptions.','Carrier field is ungauged; no electric charge or fermions have been derived.','Negative eigenvalue signals instability of Sigma=0, not stability of the extended system.','Original vortex profile is held fixed at linear order. Nonlinear backreaction and loop curvature remain untested.','Original background profile is interpolated from spacing 0.05, limiting eigenvalue accuracy.']}
(root/'outputs'/'core_carrier_linear_test.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
