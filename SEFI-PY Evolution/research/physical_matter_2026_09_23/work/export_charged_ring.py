import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
names=['charged_ring_trial','charged_ring_refined','charged_ring_wide']
records=[]
for name in names:
    rec=json.loads((root/'work'/f'{name}.json').read_text());p=np.load(root/'work'/f'{name}.npz');r=p['r'];z=p['z'];h=r[1]-r[0]
    radius=rec['core_radii'][0];theta=np.linspace(0,2*np.pi,513);rp=radius+1.5*np.cos(theta);zp=1.5*np.sin(theta)
    def sample(arr):
        i=np.minimum((rp/h).astype(int),len(r)-2);j=np.minimum((abs(zp)/h).astype(int),len(z)-2);a=rp/h-i;b=abs(zp)/h-j
        return (1-a)*(1-b)*arr[i,j]+a*(1-b)*arr[i+1,j]+(1-a)*b*arr[i,j+1]+a*b*arr[i+1,j+1]
    psi=sample(p['u'])+1j*np.sign(zp)*sample(p['w']);rec['meridional_winding']=float((np.unwrap(np.angle(psi))[-1]-np.unwrap(np.angle(psi))[0])/(2*np.pi));rec['minimum_amplitude_on_winding_contour']=float(abs(psi).min())
    s=p['s'];density=r[:,None]*s*s;dist=np.sqrt((r[:,None]-radius)**2+z[None,:]**2)
    integ=lambda f:float(4*np.pi*np.trapezoid(np.trapezoid(f,z,axis=1),r))
    rec['carrier_norm_fraction_beyond_core_distance_6']=integ(density*(dist>6))/integ(density)
    records.append(rec)
np.savez_compressed(root/'outputs'/'charged_vortex_ring_profile.npz',**dict(p))
(root/'outputs'/'charged_vortex_ring_results.json').write_text(json.dumps({'parameters':{'alpha':1,'beta':1,'Omega_phi':float(np.sqrt(2)),'kappa':4,'mu_sigma_squared':3,'lambda_sigma':25},'runs':records,'limitations':'Finite-cylinder candidates. No stability, rest-particle, quantum or universal-geometry result.'},indent=2))
print(json.dumps(records,indent=2))
