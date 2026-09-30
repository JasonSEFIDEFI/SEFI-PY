# N=0 full-Z falsification campaign — 30 September 2026

Starting commit: `bc68d639fe9b175f2e7c4fade05f991c0276157b`, branch
`research/action-audit-2026-09-30`. Existing tracked files, v0.1, retained
backgrounds, and analytical checkpoints remain unchanged. The committed
operator and results were audited before extending the search.

Parameters: N=0, c=0.6, nu=0.3, lambda=25. Full Z domain; no perturbation parity restriction.

**Result: no accepted instability; no stability claim.**
The campaign retained 2412 returned eigenpairs (including repeats)
from 202 recorded target calculations. Nonsymmetry positive-real
candidates above the 1e-6 screening threshold: **0**.
The positive m=0 pair is the previously calibrated translation splitting.
No genuine unstable mode survived the acceptance gates.

## Completed scope

| L | h | m sectors | Complete / recorded targets | Returned pairs |
|---:|---:|:---|---:|---:|
| 16 | 0.4 | 0, 1, 2, 3, 4 | 55 / 55 | 660 |
| 16 | 0.5 | 0, 1, 2, 3, 4 | 67 / 68 | 804 |
| 24 | 0.4 | 0, 1, 2, 3 | 12 / 12 | 144 |
| 24 | 0.5 | 0, 1, 2, 3, 4 | 55 / 55 | 660 |
| 32 | 0.5 | 0, 1, 2, 3 | 12 / 12 | 144 |

The main scans use eleven complex targets per sector, with imaginary
offsets through +/-0.25 and positive real offsets 0.0375 and 0.15, plus
0.0001. The L24 h=0.4 and L32 refinements use 0.01 and 0.12 +/-0.12i.
Independent full-matrix targeting adds m=0 offsets through -1 and +2,
and separately reproduces the m=2,3 near-zero/complex-target calculations.
One unshifted largest-real-part attempt returned no converged eigenpairs.
The -2 offset was interrupted without results; it is excluded from scope.

Actual returned imaginary ranges (these are not exclusion regions):

| m | Smallest Im(sigma) | Largest Im(sigma) |
|---:|---:|---:|
| 0 | -1.02005157 | 2.0071191 |
| 1 | -0.319647823 | 0.319647823 |
| 2 | -0.334722819 | 0.334722819 |
| 3 | -0.348274194 | 0.348274194 |
| 4 | -0.368174767 | 0.368174767 |

## Symmetry calibration

| L | h | Longitudinal splitting magnitude | Transverse splitting magnitude |
|---:|---:|---:|---:|
| 16 | 0.4 | 0.00269583551 | 0.00468127336 |
| 16 | 0.5 | 0.00270686543 | 0.00482831558 |
| 24 | 0.4 | 0.00104666802 | 0.00214632945 |
| 24 | 0.5 | 0.00105073302 | 0.00245837112 |
| 32 | 0.5 | 0.000519913272 | 0.00211569189 |

The largest direct carrier-phase null residual is 2.99e-15. Translation
splittings shrink with increasing L and show the independent h trend.
Cylindrical overlaps, boundary fractions, interior tangent residuals, and
background checks are saved in convergence.json. Exact infinite-domain
translation zero modes remain unresolved.

The symmetric nearest-eigenvalue disagreement between Schur and independent
full-matrix targeting at L16 h=0.5 is m=0: 3.29e-08, m=2: 2.54e-15, m=3: 3.32e-15.
Both use ARPACK; this verifies different factorizations rather than different
eigensolver libraries.

## Limits and reproducibility

Every returned pair has first-order/quadratic residuals and eigenfunction
diagnostics. Representative vectors and every nonsymmetry growing candidate
are retained. No candidate is accepted automatically; all six requested
residual, refinement, localization, symmetry, and reproduction gates remain.
Growth below the screening threshold and modes missed by finite targeting
remain possible. No continuum, nonlinear, or general spectral stability
claim follows from this finite search.

Six tests pass, including an unstable pair missed by near-zero targeting,
a dense block-inverse oracle, and a dense-spectrum check of the weighted
energy identity. That identity supplies a disk bound when the largest B
eigenvalue is known; the numerical estimate used here is not certified.
A certified bound/inertia or contour count is a useful next completeness
check. Continue simultaneous box/grid refinement before accepting growth.

[Method, targeting limits, and reproduction commands](broad_scan_strategy.md).
[Machine-readable convergence](broad_runs/convergence.json) and
[preservation/artifact verification](broad_runs/verification.json).
Per-target spectra, retained vectors, backgrounds, and environment/settings
are in broad_runs/. All changes are separate research additions.
