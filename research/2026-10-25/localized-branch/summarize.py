"""Build tables, convergence diagnostics, analytic checks and the research figure."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
def read(name):return json.loads((HERE/name).read_text())
def indexed(name):return {r['Q']:r for r in read(name+'/results.json')['runs']}
a,b,c=[indexed(name) for name in ['h0p4_L80','h0p2_L80','h0p1_L80']]
rows=[]
for Q,r in sorted(c.items()):
    ratio=(a[Q]['energy']-b[Q]['energy'])/(b[Q]['energy']-r['energy'])
    rows.append(dict(Q=Q,energy_per_charge=r['energy_per_charge'],omega=r['omega'],
        rms_radius=r['rms_radius'],binding_fraction=1-r['energy_per_charge']/np.sqrt(2),
        three_grid_difference_ratio=ratio,
        richardson_energy_per_charge=r['energy_per_charge']+(r['energy_per_charge']-b[Q]['energy_per_charge'])/3))
with (HERE/'branch_table.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
frags=read('h0p1_L80/results.json')['fragmentation']
with (HERE/'fragmentation_table.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(frags[0]));w.writeheader();w.writerows(frags)
rho,eta,lam=sp.symbols('rho eta lam',positive=True)
V=rho*rho/4+rho/2+2*rho*eta+lam*(eta-3/lam)**2/4
square=lam*(eta-3/lam+4*rho/lam)**2/4+(sp.Rational(1,4)-4/lam)*rho**2
err=sp.simplify(V-(1+12/lam)*rho/2-square)
assert err==0
summary=dict(scope='Exploratory model-dependent numerical research, not promoted',
    counts=dict(main_profiles=sum(len(indexed(x)) for x in ['h0p4_L80','h0p2_L80','h0p1_L80']),
                sampled_partitions=len(frags)),
    convergence_ratio_range=[min(r['three_grid_difference_ratio'] for r in rows),
                             max(r['three_grid_difference_ratio'] for r in rows)],
    min_sampled_fragmentation_margin=min(r['daughter_minus_parent_energy'] for r in frags),
    min_sampled_fragmentation_relative_margin=min(r['relative_margin'] for r in frags),
    symbolic_square_identity=str(err),branch_rows=rows)
if (HERE/'additional_checks.json').exists():
    extra=read('additional_checks.json')
    # A shifted target charge uses a neighboring seed, not a matched-charge FV comparison.
    def clean(obj):
        if isinstance(obj,dict):
            if 'source_folder' in obj and 'Q' in obj:
                choices=list(indexed(obj['source_folder']))
                seed=obj.get('seed_Q',min(choices,key=lambda q:abs(q-obj['Q'])))
                obj['seed_Q']=seed
                if obj['Q']!=seed:
                    obj['finite_volume_energy_per_charge']=None
                    obj['relative_energy_difference']=None
            for v in obj.values():clean(v)
        elif isinstance(obj,list):
            for v in obj:clean(v)
    clean(extra)
    (HERE/'additional_checks.json').write_text(json.dumps(extra,indent=2)+'\n')
    summary['binding_crossing']=extra.get('binding_crossing')
    summary['first_law_errors']=[dict(Q=x['Q'],errors=[y['absolute_error'] for y in x['checks']]) for x in extra['first_law']]
    summary['spectra_completed']=len(extra['spectra'])
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
q=np.array([r['Q'] for r in rows])
axs[0,0].semilogx(q,[r['energy_per_charge'] for r in rows],'o-',label='finite volume h=0.1')
axs[0,0].axhline(np.sqrt(2),color='black',ls='--',label='free threshold')
axs[0,0].set(xlabel='Charge Q',ylabel='E/Q',title='Selected spherical branch')
axs[0,0].legend(fontsize=8)
threshold=indexed('threshold_h0p05_L120')
qq=sorted(q for q in threshold if q<1500)
axs[0,1].plot(qq,[threshold[q]['energy_per_charge']-np.sqrt(2) for q in qq],'o-')
axs[0,1].axhline(0,color='black',ls='--')
axs[0,1].set(xlabel='Charge Q',ylabel='E/Q - sqrt(2)',title='Refined binding crossing (h=0.05, L=120)')
axs[1,0].scatter([r['Q'] for r in frags],[r['relative_margin']*100 for r in frags],s=20)
axs[1,0].set_xscale('log')
axs[1,0].set(xlabel='Parent charge',ylabel='Daughter minus parent energy (%)',title='26 sampled two-daughter comparisons')
for Q in [1250.,5000.,100000.]:
    er=[a[Q]['energy_per_charge'],b[Q]['energy_per_charge'],c[Q]['energy_per_charge']]
    axs[1,1].plot(np.array([.4,.2,.1])**2,1e6*(np.array(er)-er[-1]),'o-',label='Q='+format(Q,'.0f'))
axs[1,1].set(xlabel='h squared',ylabel='1e6 x (E/Q - finest E/Q)',title='Three-grid convergence, L=80')
axs[1,1].legend(fontsize=8)
fig.suptitle('Two-complex-scalar model, lambda=12 — exploratory numerical evidence',fontsize=12)
fig.savefig(HERE/'branch_evidence.png',dpi=160)
plt.close(fig)
manifest={}
for file in sorted(HERE.rglob('*')):
    if file.is_file() and file.suffix in ['.py','.json','.csv','.npz','.png','.md','.txt'] and file.name!='manifest.json':
        manifest[str(file.relative_to(HERE)).replace('\\','/')]=hashlib.sha256(file.read_bytes()).hexdigest()
(HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='branch_rows'}),flush=True)

