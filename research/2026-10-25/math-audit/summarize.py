"""Generate descriptive report, scientific plots and integrity manifest."""
from pathlib import Path
import json,hashlib,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'results.json').read_text())
ref=json.loads((ROOT/'landmark_refinement.json').read_text())
second=json.loads((ROOT/'second_alpha_results.json').read_text())
stability=json.loads((ROOT/'stability_results.json').read_text())
floquet=json.loads((ROOT/'foundations_floquet.json').read_text())
counts={}
for name in ['audit','qec','geometry','nuclear','voxel']:
    log=(ROOT/f'{name}-tests.log').read_text()
    assert not re.search(r'\d+ failed|\d+ error',log),name
    match=re.search(r'(\d+) passed',log)
    assert match,name
    counts[name]=int(match.group(1))
last=data['full'][-1]
improvement=abs(last['q_leading_error']/last['q_first_order_error'])
err0=abs(last['q_leading_error']/last['q_scaled'])
err1=abs(last['q_first_order_error']/last['q_scaled'])
orders=[float(np.log(abs(x['q_first_order_error']/y['q_first_order_error']))/np.log(2))
        for x,y in zip(data['full'][:-1],data['full'][1:])]
fold_error=abs(ref['fold']['q_scaled']/data['fold']['q_scaled']-1)
cross_error=abs(ref['energetic_crossing']['q_scaled']/data['energetic_crossing']['q_scaled']-1)
original=next(r for r in data['full'] if r['delta']==.03125)
full_domain_error=abs(data['full_refinement']['q_scaled']/original['q_scaled']-1)
resonant=[r for r in floquet['runs'] if r['sign']==-1 and r['Omega']==2]
summary=dict(tests=counts,total_tests=sum(counts.values()),
    first_correction=dict(charge_error_reduction_factor=improvement,
        leading_relative_charge_error=err0,corrected_relative_charge_error=err1,
        observed_successive_orders=orders,second_alpha_improvements=[
            abs(r['q_leading_error']/r['q_first_order_error']) for r in second['full']]),
    convergence=dict(fold_scaled_charge_relative_change=fold_error,
                     binding_crossing_scaled_charge_relative_change=cross_error,
                     full_scaled_charge_domain_tolerance_relative_change=full_domain_error),
    scope='Descriptive deterministic convergence; no statistical intervals or new-physics claim')
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2))

fig,axes=plt.subplots(1,3,figsize=(14,4.2))
scan=data['scan']
axes[0].plot([r['alpha'] for r in scan],[r['q_scaled']/1e4 for r in scan],color='#225ea8')
for key,label,color in [('fold','Minimum charge','#d95f02'),('energetic_crossing','Binding crossing','#238b45')]:
    r=data[key]
    axes[0].scatter(r['alpha'],r['q_scaled']/1e4,color=color,label=label,zorder=4)
axes[0].set(xlabel='Reduced frequency parameter alpha',ylabel='Scaled charge / 10,000',title='Two distinct branch landmarks')
axes[0].legend(fontsize=8)
for key,label,color in [('q_leading_error','Leading approximation','#d95f02'),
                         ('q_first_order_error','Derived first correction','#225ea8')]:
    axes[1].loglog([r['delta'] for r in data['full']],
        [abs(r[key]/r['q_scaled']) for r in data['full']],'o-',label=label,color=color)
axes[1].set(xlabel='delta = 16 - lambda',ylabel='Relative scaled-charge error',
            title='Prediction at fixed alpha')
axes[1].legend(fontsize=8)
eps=np.array([r['epsilon'] for r in resonant])
axes[2].plot(eps,[r['growth_exponent'] for r in resonant],'o',label='Monodromy calculation',color='#225ea8')
axes[2].plot(eps,eps/4,'--',label='Mathieu prediction: epsilon / 4',color='#d95f02')
axes[2].set(xlabel='Geometry modulation epsilon',ylabel='Growth exponent',title='Positive Hessian, resonant growth')
axes[2].legend(fontsize=8)
fig.tight_layout();fig.savefig(ROOT/'mathematical_continuation.png',dpi=180);plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(10,4.2))
for branch in stability['branches']:
    grids=branch['grids'];center=branch['reduced_center']
    axes[0].plot([g['h']**2 for g in grids],
        [g['sectors'][0]['real_eigenvalues'][0] for g in grids],'o-',label=f'Center parameter {center}')
    for field,label in [('real_eigenvalues','Translation'),('phi_phase_eigenvalues','Phase')]:
        sector=1 if field=='real_eigenvalues' else 0
        axes[1].plot([g['h']**2 for g in grids],
            [g['sectors'][sector][field][0] for g in grids],'o-',label=f'{label}, center {center}')
axes[0].axhline(0,color='gray',lw=.8)
axes[0].set(xlabel='Squared grid spacing h²',ylabel='Lowest scaled curvature',title='Radial saddle versus positive curvature')
axes[1].axhline(0,color='gray',lw=.8)
axes[1].set(xlabel='Squared grid spacing h²',ylabel='Lowest scaled symmetry curvature',title='Symmetry artifacts approach zero')
for ax in axes:ax.legend(fontsize=8)
fig.tight_layout();fig.savefig(ROOT/'energy_curvature.png',dpi=180);plt.close(fig)

text=f"""# Completed audit and continuation

This is a mathematical/numerical audit of the principal claims developed in this chat, not a certification of every repository manuscript. See MATH_CONTINUED.md for proofs, assumptions, corrections and source context.

## Results

- Exact original and generalized quartic square-completion identities pass. The general threshold is a*lambda=c²; the repository coefficients give lambda=16.
- The reduced charge fold is q={ref['fold']['q_scaled']:.5f}, alpha={ref['fold']['alpha']:.10f}, with positive scaled binding {ref['fold']['binding_scaled']:.6f}.
- The reduced free-particle binding crossing is q={ref['energetic_crossing']['q_scaled']:.5f}, alpha={ref['energetic_crossing']['alpha']:.10f}. These are reduced-limit constants, not measured charges or particle properties.
- Enlarging the domain and tightening tolerance changes these charges relatively by {fold_error:.3g} and {cross_error:.3g}. This estimates numerical sensitivity, not a rigorous error bound or statistical uncertainty.
- At the closest tested delta={last['delta']}, the independent first correction reduces relative charge error from {100*err0:.6f}% to {100*err1:.6f}%, a factor of {improvement:.1f}. Successive remainder orders are {', '.join(f'{x:.3f}' for x in orders)}, consistent with O(delta²).
- A second alpha also improves at every tested delta; the reduction factors are {', '.join(f'{x:.1f}' for x in summary['first_correction']['second_alpha_improvements'])}. No coefficients are fitted to full-model results.
- The full-model domain/tolerance refinement changes scaled charge relatively by {full_domain_error:.3g}. Parameter-derivative refinement separately checks the envelope relation.
- The lower-amplitude test branch has a persistent negative radial curvature. The higher-amplitude test branch has positive radial curvature; phase/translation zero-mode artifacts decrease approximately as h². These are static energy-curvature results, not nonlinear stability.
- Under explicit positivity, monotonicity and boundary assumptions, ground-state identities exclude negative phase and nonradial amplitude energy directions analytically. Radial fixed-charge curvature remains the decisive block.
- In the foundational prescribed-geometry model, all instantaneous stiffnesses are positive but periodic modulation at resonance grows. At epsilon=.02 the calculated exponent is {resonant[0]['growth_exponent']:.10f}, versus the known leading Mathieu prediction .005.

## Verification

{sum(counts.values())} selected tests pass: {json.dumps(counts)}. There are 43 reduced scan profiles, five primary fixed-alpha full solves, one full domain/tolerance refinement, three second-alpha full solves, and two additional full profiles for curvature. Four curvature grids yield 24 lowest-eigenvalue blocks. The Floquet sweep has 24 conditions plus one refinement.

Configurations, profiles, pre-evaluation correction coefficients, timings, environments, raw results and failure logs are retained. The reduced coefficients were computed before full-model evaluation, but this is not an external preregistration. The failed L=320 virial check and inadmissible direct L=480 initialization are not hidden.

## What this leads to

The clearest next question is a controlled near-critical binding transition in a specified classical model. A useful quantitative reduction now exists, including a tested first correction and a separation of branch existence, energetic binding and quadratic stability.

Next: continue full-model charge folds and radial eigenvalue zero crossings across delta; test the first-correction validity near both frequency endpoints; perform charge-conserving time evolution; and compare the specific formulas with prior quartic-extended-model literature.

For the foundational worldline action, a constrained Hamiltonian and a physical observation map are still required. The scalar model has not been derived from it. No new element, planet, experimentally demonstrated matter, unique SEFI advantage, or novelty claim is established.
"""
(ROOT/'REPORT.md').write_text(text,encoding='utf-8')
files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.iterdir()
       if p.is_file() and p.name!='SHA256SUMS.json'}
(ROOT/'SHA256SUMS.json').write_text(json.dumps(files,indent=2))
print(json.dumps(summary,indent=2))
