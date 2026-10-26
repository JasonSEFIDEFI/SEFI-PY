from pathlib import Path
import json,hashlib,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'results.json').read_text())
radial=json.loads((ROOT/'radial_results.json').read_text())
sensitivity=json.loads((ROOT/'sensitivity_results.json').read_text())
assert len(data['folds'])==6 and len(radial)==12
assert len(sensitivity)==2
for r in sensitivity:
    assert abs(r['corrected_tangent']['q_c'])<1e-6
    assert all(a['identity_residual_per_tangent']>b['identity_residual_per_tangent']
               for a,b in zip(r['checks'][:-1],r['checks'][1:]))
assert '3 passed' in (ROOT/'tests.log').read_text()
folds=data['folds'];last=folds[-1];pred=data['predictions']
assert all(f['admissible'] and f['charge_second_derivative_c']>0 for f in folds)
assert all(abs(f['q_first_error'])<abs(f['q_leading_error']) for f in folds)
orders=[float(np.log(abs(a['q_first_error']/b['q_first_error']))/np.log(2))
        for a,b in zip(folds[:-1],folds[1:])]
reference=next(f for f in folds if f['delta']==.0625)
ref_error=abs(data['refinement']['q_scaled']/reference['q_scaled']-1)
improvement=abs(last['q_leading_error']/last['q_first_error'])
crossings=[]
for delta in [.25,.03125]:
    below=next(r for r in radial if r['delta']==delta and r['center_offset']<0 and r['h']==.25)
    above=next(r for r in radial if r['delta']==delta and r['center_offset']>0 and r['h']==.25)
    zero=[next(r for r in radial if r['delta']==delta and r['center_offset']==0 and r['h']==h) for h in [.5,.25]]
    assert below['scaled_lowest_eigenvalues'][0]<0<above['scaled_lowest_eigenvalues'][0]
    crossings.append(dict(delta=delta,below=below['scaled_lowest_eigenvalues'][0],
        above=above['scaled_lowest_eigenvalues'][0],
        at_fold=zero[-1]['scaled_lowest_eigenvalues'][0],
        zero_mode_refinement_ratio=zero[0]['scaled_lowest_eigenvalues'][0]/zero[1]['scaled_lowest_eigenvalues'][0],
        tangent_overlap=zero[-1]['lowest_vector_tangent_overlap'],
        family_identity_residual_ratio=zero[0]['family_identity_residual_per_tangent']/zero[1]['family_identity_residual_per_tangent']))
summary=dict(tests_passed=3,fold_refinement_relative_change=ref_error,
    closest_delta=last['delta'],closest_q_error_reduction=improvement,
    closest_leading_relative_error=abs(last['q_leading_error']/last['q_scaled']),
    closest_corrected_relative_error=abs(last['q_first_error']/last['q_scaled']),
    correction_remainder_orders=orders,effective_Qmin_exponents=data['effective_exponents'],
    radial_crossings=crossings,full_fold_bvp_solves=data['full_solve_count'],
    matching_domain_certificates=sensitivity)
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2))

fig,axes=plt.subplots(2,2,figsize=(11,8))
d=np.array([f['delta'] for f in folds])
axes[0,0].semilogx(d,[f['q_scaled'] for f in folds],'o',label='Original equations')
axes[0,0].semilogx(d,pred['q0']+d*pred['q1'],'--',label='Derived first correction')
axes[0,0].axhline(pred['q0'],color='gray',ls=':',label='Leading critical limit')
axes[0,0].set(xlabel='delta = 16 - lambda',ylabel='Q_min × delta^(3/2)',title='Full-model minimum-charge scaling')
for key,label in [('q_leading_error','Leading approximation'),('q_first_error','First correction')]:
    axes[0,1].loglog(d,[abs(f[key]) for f in folds],'o-',label=label)
axes[0,1].set(xlabel='delta = 16 - lambda',ylabel='Absolute scaled-charge error',title='Prediction made before full-model evaluation')
for delta in [.25,.03125]:
    for h,style in [(.5,'--'),(.25,'-')]:
        selected=sorted([r for r in radial if r['delta']==delta and r['h']==h],key=lambda r:r['center_offset'])
        axes[1,0].plot([r['center_offset'] for r in selected],
                      [r['scaled_lowest_eigenvalues'][0] for r in selected],style+'o',label=f'delta={delta}, h={h}')
axes[1,0].axhline(0,color='gray',lw=.7)
axes[1,0].set(xlabel='Central parameter minus fold parameter',ylabel='Lowest scaled radial curvature',title='Radial energy-curvature sign crossing')
mid=np.sqrt(d[:-1]*d[1:])
axes[1,1].semilogx(mid,data['effective_exponents'],'o-',label='Observed adjacent-pair exponent')
axes[1,1].axhline(1.5,color='gray',ls='--',label='Derived critical exponent 3/2')
axes[1,1].set(xlabel='Geometric-mean delta of adjacent pair',ylabel='Effective exponent of Q_min',title='Approach to the critical exponent')
for ax in axes.flat:ax.legend(fontsize=7)
fig.tight_layout();fig.savefig(ROOT/'full_fold_results.png',dpi=180);plt.close(fig)

table='\n'.join(f"| {f['delta']:.6f} | {f['q_scaled']:.6f} | {f['Q']:.3f} | {f['alpha']:.10f} | {f['q_first_error']:.6f} |" for f in folds)
text=f"""# Full-model fold session — 26 October 2026

This session extends the Oct 25 critical-binding audit. It uses the original dimensionless two-field equations; no physical-matter or novelty claim is made.

## Independently derived prediction

With delta=16-lambda,

Q_min = delta^(-3/2) [q0 + delta*q1 + O(delta²)],

q0 = {pred['q0']:.8f}, q1 = {pred['q1']:.8f}.

The fold frequency parameter satisfies alpha_fold=alpha0+delta*alpha1+O(delta²), with alpha0={pred['alpha0']:.12f}, alpha1={pred['alpha1']:.12f}.
These coefficients were evaluated from reduced solutions before the new full-model sweep. Finite-difference derivative refinements are saved; this is not external preregistration.

## Observed full-model folds

| delta | Scaled minimum charge | Q_min | alpha_fold | Corrected scaled-charge error |
|---|---:|---:|---:|---:|
{table}

Adjacent-pair effective exponents are {', '.join(f'{x:.6f}' for x in data['effective_exponents'])}, approaching the predicted 1.5.
The corrected remainder orders are {', '.join(f'{x:.4f}' for x in orders)}, approximately second order.
At the closest delta the charge-prediction error falls from {100*summary['closest_leading_relative_error']:.6f}% to {100*summary['closest_corrected_relative_error']:.6f}%, a factor of {improvement:.1f}.
The longer-domain/tighter-tolerance fold repeat changes scaled charge relatively by {ref_error:.3g}. This is numerical sensitivity, not a rigorous error bound or statistical interval.

The charge minimum is a local fold of the selected positive monotone spherical branch, not a global minimum over all fields. Its scaled binding remains positive at every tested coupling, so reaching the minimum charge does not by itself give energetic binding below free particles.

## Radial curvature

Two couplings, three central amplitudes each, and two grids produce 12 lowest-radial-eigenvalue blocks. Below each fold the lowest eigenvalue is negative; above it it is positive.
At the fold, small negative grid eigenvalues shrink toward zero; their coarse/fine ratios are {', '.join(f"{r['zero_mode_refinement_ratio']:.4f}" for r in crossings)}.
The corresponding eigenvectors overlap the stationary-family tangents by {', '.join(f"{r['tangent_overlap']:.9f}" for r in crossings)}.
The exact identity H_Q F_c = (2 omega Q_c / I)(f,0) explains why the charge fold has a radial zero mode. The radial matrix uses scaled coordinates, so its eigenvalues are physical curvature divided by delta.

The first family-identity residual used a shorter spectral domain (240) than the stationary boundary-value domain (320); it worsened on one refinement because the unmodified tangent did not obey the shorter-domain outer boundary. Those diagnostics remain in radial_results.json. A new linear differentiated-equation solve on the matching domain (320) supplies the tangent without finite parameter differences. One derivative-based Newton correction reduces absolute charge derivatives to {', '.join(f"{abs(r['corrected_tangent']['q_c']):.3g}" for r in sensitivity)}. On three grids (.5, .25, .125), matching-domain identity residuals decrease roughly fourfold with each halving, reaching {', '.join(f"{r['checks'][-1]['identity_residual_per_tangent']:.3g}" for r in sensitivity)}. This supports the zero-mode identity while making the original boundary mismatch explicit.

These results support a simple local quadratic-energy index crossing. They are not a charge-conserving dynamical simulation or a nonlinear stability theorem.

## Verification and limitations

Three new mathematical tests pass. The fold search used {data['full_solve_count']} boundary-value solves including the domain/tolerance repeat; the radial calculation records its separate profile solves. ODE residuals, off-node derivative checks, virial residuals, monotonicity, frequency bounds, timings and optimization diagnostics are retained.

An initial empty-parameter seed interface error stopped the first sweep before any full fold was evaluated. Its log is retained. Corrected initialization uses a reduced profile with its frequency parameter.

The first-correction coefficient has finite parameter-difference error; do not interpret all printed digits as guaranteed precision. Only one full-fold domain refinement and two eigenvalue grids per sampled profile were run; the matching-domain identity uses three grids but is not an additional eigenvalue computation. Extrapolation to either frequency endpoint, other branches, gauged fields, complex Sigma, quantum behavior or measured nuclei is unsupported.

This strengthens the controlled mathematical description of the model's binding transition. The next discriminating step is charge-conserving time evolution near the radial fold, followed by comparison of these specific formulas with prior quartic-extended soliton literature. The separate foundational worldline-to-field and field-to-observation bridges remain open.
"""
(ROOT/'REPORT.md').write_text(text,encoding='utf-8')
print(json.dumps(summary,indent=2))
