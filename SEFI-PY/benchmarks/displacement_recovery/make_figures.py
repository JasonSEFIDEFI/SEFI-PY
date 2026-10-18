"""Three data-driven scientific figures, PNG and editable SVG. Jason D. Dutton."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                     'axes.spines.right':False,'svg.fonttype':'none','savefig.facecolor':'white'})
rows=list(csv.DictReader((ROOT/'results/table_1_all_results.csv').open()))
runtime=list(csv.DictReader((ROOT/'results/runtime.csv').open()))
COLORS={'ML':'#596579','Geometry':'#006b8f','Conventional posterior':'#cf8f16','No recovery':'#aa4b64'}

def selected(name):
    return sorted([r for r in rows if r['method']==name],key=lambda r:float(r['sigma']))

def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=220,bbox_inches='tight')
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight');plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(12,4.8),gridspec_kw={'width_ratios':[1.2,1]})
for name in ('No recovery','ML','Geometry','Conventional posterior'):
    s=selected(name);x=np.array([float(r['sigma']) for r in s]);y=np.array([float(r['entanglement_fidelity']) for r in s])
    lo=np.array([float(r['entanglement_fidelity_ci_low']) for r in s]);hi=np.array([float(r['entanglement_fidelity_ci_high']) for r in s])
    style='--' if name=='Conventional posterior' else '-'
    axes[0].plot(x,y,style,marker='o' if name!='Conventional posterior' else None,lw=2,color=COLORS[name],label=name)
    axes[0].fill_between(x,lo,hi,color=COLORS[name],alpha=.12)
s=selected('Geometry');x=np.array([float(r['sigma']) for r in s]);d=np.array([float(r['paired_delta_fe_vs_ml']) for r in s]);lo=np.array([float(r['paired_delta_ci_low']) for r in s]);hi=np.array([float(r['paired_delta_ci_high']) for r in s])
axes[1].errorbar(x,100*d,yerr=np.vstack([d-lo,hi-d])*100,color=COLORS['Geometry'],marker='o',capsize=4,lw=2)
axes[1].axhline(0,color='#666',lw=1)
axes[0].set(xlabel='Displacement noise scale sigma (complex-amplitude units)',ylabel='Unconditional entanglement fidelity',ylim=(0,1.03),title='A  |  Encoded oscillator qubit')
axes[1].set(xlabel='Displacement noise scale sigma',ylabel='Geometry minus ML (percentage points)',title='B  |  Paired effect with 95% seed-bootstrap CI')
axes[0].legend(frameon=False,fontsize=9,loc='lower left')
for ax in axes:ax.grid(alpha=.15);ax.axvline(.2,color='#777',ls=':',alpha=.5)
fig.suptitle('Pilot-assisted recovery: the geometry estimator matches conventional Bayesian inference',fontsize=13,y=1.02)
fig.text(.5,-.035,'20 seeds x 4,000 held-out trials per noise | 16 shared-noise heterodyne probes | leakage retained as failure',ha='center',fontsize=9)
fig.tight_layout();save(fig,'figure_1_noise_performance')

names=['Conventional posterior','No prior distance','No anisotropy','Random metric orientation','Shuffled features','Affine coordinates','Nonlinear coordinates']
s=[next(r for r in rows if float(r['sigma'])==.2 and r['method']==name) for name in names]
d=np.array([float(r['paired_delta_fe_vs_geometry']) for r in s])*100
lo=np.array([float(r['paired_delta_geometry_ci_low']) for r in s])*100;hi=np.array([float(r['paired_delta_geometry_ci_high']) for r in s])*100
fig,axes=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.4,1]})
axes[0].barh(np.arange(len(names)),d,color=['#cf8f16' if i==0 else '#596579' for i in range(len(names))],alpha=.85)
axes[0].errorbar(d,np.arange(len(names)),xerr=np.vstack([np.maximum(d-lo,0),np.maximum(hi-d,0)]),fmt='none',ecolor='#222',capsize=4)
axes[0].set_yticks(np.arange(len(names)),names,fontsize=10);axes[0].invert_yaxis();axes[0].axvline(0,color='#222',lw=1)
axes[0].plot(d,np.arange(len(names)),'o',color='#222',markersize=4)
axes[0].set_xlim(min(lo)-.3,.35)
axes[0].set(xlabel='Change in fidelity versus full geometry (percentage points)',title='A  |  Ablations and representation controls at sigma = 0.2')
stress=list(csv.DictReader((ROOT/'results/distribution_shift.csv').open()))
for i,name in enumerate(('ML','Geometry')):
    matched=next(r for r in rows if float(r['sigma'])==.2 and r['method']==name)
    shifted=next(r for r in stress if r['method']==name)
    axes[1].bar(np.array([0,1])+i*.32,[float(matched['entanglement_fidelity']),float(shifted['fe'])],width=.3,color=COLORS[name],label=name)
axes[1].set_xticks([.16,1.16],['Matched prior','Unmodeled mean shift']);axes[1].set(ylim=(0,1),ylabel='Unconditional entanglement fidelity',title='B  |  Prior mismatch reverses the gain')
axes[1].legend(frameon=False,loc='lower left',fontsize=9);axes[1].grid(axis='y',alpha=.15)
fig.suptitle('Controls separate ordinary prior shrinkage from any geometry-specific claim',fontsize=13,y=1.02)
fig.text(.5,-.025,'Randomizing orientation weakens the gain; it does not eliminate shrinkage. Shuffling trial features destroys their alignment.',ha='center',fontsize=9)
fig.tight_layout();save(fig,'figure_2_ablation_controls')

names=['ML','Geometry','Conventional posterior','No anisotropy','Random metric orientation','Shuffled features','Affine coordinates','Nonlinear coordinates']
s=[next(r for r in runtime if r['method']==n) for n in names]
med=np.array([float(r['median_microseconds_per_trial']) for r in s]);lo=np.array([float(r['timing_p025']) for r in s]);hi=np.array([float(r['timing_p975']) for r in s])
fig,ax=plt.subplots(figsize=(10,5))
ax.errorbar(med,np.arange(len(names)),xerr=np.vstack([med-lo,hi-med]),fmt='o',color='#006b8f',capsize=4)
ax.set_yticks(np.arange(len(names)),names);ax.invert_yaxis();ax.set_xscale('log');ax.grid(axis='x',alpha=.2)
ax.set(xlabel='Batched decoder microseconds per trial (log scale)',title='Figure 3  |  Decode-only runtime on this computer')
fig.text(.5,-.01,'Median and empirical 2.5-97.5% timing range, 51 batches of 4,000 | preparation/scoring excluded | not hardware latency',ha='center',fontsize=9)
fig.tight_layout();save(fig,'figure_3_runtime')
print('Created all three figures as PNG and SVG.')
