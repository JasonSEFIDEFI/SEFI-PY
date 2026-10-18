"""Write the concise campaign report from finalized convergence diagnostics."""
import argparse
import json
from pathlib import Path


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('root',type=Path)
    args=ap.parse_args()
    d=json.loads((args.root/'convergence.json').read_text())
    groups={}
    for run in d['runs']:
        group=groups.setdefault((run['L'],run['h']),{'targets':0,'complete':0,'pairs':0,'m':set()})
        for m,s in run['sectors'].items():
            group['m'].add(int(m));group['targets']+=s['targets']
            group['complete']+=s['complete_targets'];group['pairs']+=s['returned_pairs']
    lines=['# N=0 full-Z falsification campaign — 30 September 2026','',
        'Starting commit: `bc68d639fe9b175f2e7c4fade05f991c0276157b`, branch',
        '`research/action-audit-2026-09-30`. Existing tracked files, v0.1, retained',
        'backgrounds, and analytical checkpoints remain unchanged. The committed',
        'operator and results were audited before extending the search.','',
        'Parameters: N=0, c=0.6, nu=0.3, lambda=25. Full Z domain; no perturbation parity restriction.','',
        '**Result: no accepted instability; no stability claim.**',
        f'The campaign retained {d["total_returned_pairs"]} returned eigenpairs (including repeats)',
        f'from {d["total_targets"]} recorded target calculations. Nonsymmetry positive-real',
        f'candidates above the 1e-6 screening threshold: **{d["nonsymmetry_candidate_count"]}**.',
        'The positive m=0 pair is the previously calibrated translation splitting.',
        'No genuine unstable mode survived the acceptance gates.','',
        '## Completed scope','',
        '| L | h | m sectors | Complete / recorded targets | Returned pairs |',
        '|---:|---:|:---|---:|---:|']
    for (L,h),g in sorted(groups.items()):
        lines.append(f'| {L:g} | {h:g} | {", ".join(map(str,sorted(g["m"])))} | {g["complete"]} / {g["targets"]} | {g["pairs"]} |')
    lines += ['',
        'The main scans use eleven complex targets per sector, with imaginary',
        'offsets through +/-0.25 and positive real offsets 0.0375 and 0.15, plus',
        '0.0001. The L24 h=0.4 and L32 refinements use 0.01 and 0.12 +/-0.12i.',
        'Independent full-matrix targeting adds m=0 offsets through -1 and +2,',
        'and separately reproduces the m=2,3 near-zero/complex-target calculations.',
        'One unshifted largest-real-part attempt returned no converged eigenpairs.',
        'The -2 offset was interrupted without results; it is excluded from scope.','',
        'Actual returned imaginary ranges (these are not exclusion regions):','',
        '| m | Smallest Im(sigma) | Largest Im(sigma) |',
        '|---:|---:|---:|']
    for m in range(5):
        sectors=[r['sectors'][str(m)] for r in d['runs'] if str(m) in r['sectors'] and r['sectors'][str(m)]['returned_pairs']]
        lines.append(f'| {m} | {min(s["returned_imag_range"][0] for s in sectors):.9g} | {max(s["returned_imag_range"][1] for s in sectors):.9g} |')
    lines += ['','## Symmetry calibration','',
        '| L | h | Longitudinal splitting magnitude | Transverse splitting magnitude |',
        '|---:|---:|---:|---:|']
    for r in d['runs']:
        if r['run'] not in ('L16_h050','L16_h040','L24_h050','L24_h040','L32_h050'):
            continue
        cal=r['calibration'];a=cal['longitudinal_translation']['sigma'];b=cal['transverse_translation']['sigma']
        lines.append(f'| {r["L"]:g} | {r["h"]:g} | {abs(complex(*a)):.9g} | {abs(complex(*b)):.9g} |')
    phase=max(r['calibration']['carrier_phase']['direct_null_relative_l2'] for r in d['runs'] if 'carrier_phase' in r['calibration'])
    agreement=d['independent_strategy_agreement']
    lines += ['',
        f'The largest direct carrier-phase null residual is {phase:.3g}. Translation',
        'splittings shrink with increasing L and show the independent h trend.',
        'Cylindrical overlaps, boundary fractions, interior tangent residuals, and',
        'background checks are saved in convergence.json. Exact infinite-domain',
        'translation zero modes remain unresolved.','',
        'The symmetric nearest-eigenvalue disagreement between Schur and independent',
        'full-matrix targeting at L16 h=0.5 is '+', '.join(f'm={m}: {v["symmetric_max_nearest_eigenvalue_distance"]:.3g}' for m,v in agreement.items())+'.',
        'Both use ARPACK; this verifies different factorizations rather than different',
        'eigensolver libraries.','',
        '## Limits and reproducibility','',
        'Every returned pair has first-order/quadratic residuals and eigenfunction',
        'diagnostics. Representative vectors and every nonsymmetry growing candidate',
        'are retained. No candidate is accepted automatically; all six requested',
        'residual, refinement, localization, symmetry, and reproduction gates remain.',
        'Growth below the screening threshold and modes missed by finite targeting',
        'remain possible. No continuum, nonlinear, or general spectral stability',
        'claim follows from this finite search.','',
        'Six tests pass, including an unstable pair missed by near-zero targeting,',
        'a dense block-inverse oracle, and a dense-spectrum check of the weighted',
        'energy identity. That identity supplies a disk bound when the largest B',
        'eigenvalue is known; the numerical estimate used here is not certified.',
        'A certified bound/inertia or contour count is a useful next completeness',
        'check. Continue simultaneous box/grid refinement before accepting growth.','',
        '[Method, targeting limits, and reproduction commands](broad_scan_strategy.md).',
        '[Machine-readable convergence](broad_runs/convergence.json) and',
        '[preservation/artifact verification](broad_runs/verification.json).',
        'Per-target spectra, retained vectors, backgrounds, and environment/settings',
        'are in broad_runs/. All changes are separate research additions.','']
    Path(__file__).with_name('broad_scan_report.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':
    main()
