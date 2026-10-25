"""Generate comparison tables, symbolic checks, figure and file manifest."""
from pathlib import Path
import csv,hashlib,json,os
import numpy as np
import sympy as sp
HERE=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(HERE.parent.parent/'work'/'mpl-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def read(path):return json.loads((HERE/path).read_text())
scan=read('scan_X32/results.json');rec=read('critical_recovery/results.json')
ver=read('verification_refined/results.json');lim=read('critical_limit_refined.json')
folds=[r for r in scan['folds'] if r['lambda_sigma']<15.5]+rec['folds']
bind={r['lambda_sigma']:r for r in ver['binding_crossings']}
rows=[]
for entry in folds:
    p=entry['profile'];lam=entry['lambda_sigma'];delta=p['delta']
    rows.append(dict(lambda_sigma=lam,beta_fold=entry['beta'],Q_fold=entry['Q'],
        omega_fold=p['omega'],Q_binding=bind[lam]['Q'],beta_binding=bind[lam]['beta'],
        delta=delta,scaled_fold=entry['Q']*delta**1.5,
        scaled_binding=bind[lam]['Q']*delta**1.5,
        fold_limit_relative_error=(entry['Q']*delta**1.5/lim['fold_prediction']['scaled_charge']-1),
        binding_limit_relative_error=(bind[lam]['Q']*delta**1.5/lim['binding_prediction']['scaled_charge']-1)))
with (HERE/'coupling_table.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
qrows=[]
for r in scan['runs']+rec['runs']:
    qrows.append({k:r[k] for k in ['lambda_sigma','beta','omega','Q','energy_per_charge','accepted','status','max_rms_residual']})
with (HERE/'frequency_table.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(qrows[0]));writer.writeheader();writer.writerows(qrows)
rho,eta,t,beta=sp.symbols('rho eta t beta',positive=True)
lam=t*t;delta=3*(4-t)/t**2
V=rho**2/4+rho/2+2*rho*eta+lam*(eta-3/lam)**2/4
m2=1+12/lam
square=lam*(eta-3/lam+4*rho/lam)**2/4+(sp.Rational(1,4)-4/lam)*rho**2
identity=sp.simplify(V-m2*rho/2-square)
valley=sp.simplify(((V-(m2-delta*beta)*rho/2)/delta).subs(eta,(3-4*rho)/lam))
valley_check=sp.simplify(valley-(beta*rho/2-(4+t)*rho*rho/12))
theta=sp.symbols('theta',real=True)
metric=sp.simplify(sp.diff(sp.sqrt(3)/2*sp.sin(theta),theta)**2+
                   sp.diff(sp.sqrt(3)/4*sp.cos(theta),theta)**2)
metric_check=sp.trigsimp(metric-sp.Rational(3,16)*(1+3*sp.cos(theta)**2))
assert identity==0 and valley_check==0 and metric_check==0
hess=read('paired_hessian_refined.json')
summary=dict(scope='Exploratory / Not Promoted; model-dependent classical stationary calculations',
    initial_attempts=len(scan['runs']),initial_rejections=sum(not r['accepted'] for r in scan['runs']),
    recovery_attempts=len(rec['runs']),recovery_rejections=sum(not r['accepted'] for r in rec['runs']),
    refined_fold_lambda12=ver['refined_fold'],paired_states=ver['paired_states'],
    reduced_fold_prediction=lim['fold_prediction'],reduced_binding_prediction=lim['binding_prediction'],
    coupling_rows=rows,symbolic_checks=dict(potential_square=str(identity),
        restricted_potential=str(valley_check),kinetic_metric=str(metric_check)),
    hessian_lowest=[{k:r[k] for k in ['branch_label','h','L']}|dict(lowest=r['eigenvalues'][0]) for r in hess['runs']])
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
valid=[r for r in scan['runs'] if r['accepted'] and r['lambda_sigma']==12 and r['beta']<=.4]
valid.sort(key=lambda r:r['omega'])
axs[0,0].plot([r['omega'] for r in valid],[r['Q'] for r in valid],'o-',ms=4)
axs[0,0].scatter([ver['refined_fold']['omega']],[ver['refined_fold']['Q']],color='red',label='charge minimum')
axs[0,0].set(xlabel='Frequency omega',ylabel='Charge Q',title='lambda=12: both sides of charge turn',ylim=(850,3000))
axs[0,0].legend(fontsize=8)
axs[0,1].semilogy([r['lambda_sigma'] for r in rows],[r['Q_fold'] for r in rows],'o-',label='charge minimum')
axs[0,1].semilogy([r['lambda_sigma'] for r in rows],[r['Q_binding'] for r in rows],'s-',label='free-threshold crossing')
axs[0,1].set(xlabel='Coupling lambda',ylabel='Charge',title='Distinct thresholds approaching lambda=16')
axs[0,1].legend(fontsize=8)
axs[1,0].plot([16-r['lambda_sigma'] for r in rows],[r['scaled_fold'] for r in rows],'o-',label='full model minimum')
axs[1,0].plot([16-r['lambda_sigma'] for r in rows],[r['scaled_binding'] for r in rows],'s-',label='full model binding')
axs[1,0].axhline(lim['fold_prediction']['scaled_charge'],ls='--',color='C0',label='reduced minimum prediction')
axs[1,0].axhline(lim['binding_prediction']['scaled_charge'],ls='--',color='C1',label='reduced binding prediction')
axs[1,0].set(xlabel='16 - lambda',ylabel='Q times delta^(3/2)',title='Formal critical-limit predictions')
axs[1,0].legend(fontsize=7)
for label,title in [('high_frequency','high frequency'),('low_frequency','low frequency')]:
    rr=[r for r in hess['runs'] if r['branch_label']==label and r['L']==120]
    axs[1,1].plot([r['h']**2 for r in rr],[r['eigenvalues'][0] for r in rr],'o-',label=title)
axs[1,1].axhline(0,color='black',lw=1)
axs[1,1].set(xlabel='h squared',ylabel='Lowest radial energy curvature',title='Equal Q=1000, different energy character')
axs[1,1].legend(fontsize=8)
fig.suptitle('Frequency continuation and near-critical scaling — exploratory research',fontsize=12)
fig.savefig(HERE/'frequency_evidence.png',dpi=160);plt.close(fig)
manifest={}
for f in sorted(HERE.rglob('*')):
    if f.is_file() and f.suffix in ['.py','.json','.csv','.npz','.png','.md','.txt'] and f.name!='manifest.json' and not f.relative_to(HERE).as_posix().startswith('verification/'):
        manifest[str(f.relative_to(HERE)).replace('\\','/')]=hashlib.sha256(f.read_bytes()).hexdigest()
(HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(coupling_rows=rows,reduced_fold=lim['fold_prediction']['scaled_charge'],
    reduced_binding=lim['binding_prediction']['scaled_charge'],symbolic=summary['symbolic_checks'])),flush=True)

