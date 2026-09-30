# Three-track research checkpoint — 30 September 2026

Started from `9b3978c521f4cde792bc8a8f669dee1991887b9e` on `research/action-audit-2026-09-30`. v0.1 and every preexisting tracked file are preserved. N=0, lambda=25; no matter or gravity premise is adopted.

**Executed:** 86 converged full-Z profiles, 12 box/grid convergence cases, c continuation over 0.50–0.70, nu over 0.26–0.34, 20 centered identity checks and 16 local return checks. Newton failures: 0. Ten tests pass.

## 1. Continuum and domain convergence

Audited physical E, P and both phase charges come directly from the laboratory stress/current definitions. Bare reduced momentum, graph-functional momentum and physical momentum are distinct. The algebraic generator identity is checked before interpreting branch trends.

| L | h | E excess | P physical | Q Sigma | Carrier RMS radius |
|---:|---:|---:|---:|---:|---:|
| 16 | 0.5 | -201.94748 | 822.05411 | 7.1334546 | 5.7384294 |
| 16 | 0.4 | -200.68114 | 823.54649 | 7.1546523 | 5.736908 |
| 16 | 0.333333 | -199.97192 | 824.36655 | 7.1658915 | 5.7360588 |
| 24 | 0.5 | -205.99157 | 764.62648 | 6.9700287 | 5.5553697 |
| 24 | 0.4 | -204.84688 | 766.32614 | 6.9918567 | 5.5552256 |
| 24 | 0.333333 | -204.20401 | 767.25817 | 7.0034491 | 5.5551247 |
| 32 | 0.5 | -204.71454 | 752.71393 | 6.9360433 | 5.5166596 |
| 32 | 0.4 | -203.62291 | 754.4409 | 6.9579756 | 5.5167855 |
| 32 | 0.333333 | -203.00888 | 755.38727 | 6.969627 | 5.5168312 |
| 40 | 0.5 | -203.40477 | 748.63586 | 6.9244055 | 5.503345 |
| 40 | 0.4 | -202.34181 | 750.36873 | 6.9463708 | 5.5035612 |
| 40 | 0.333333 | -201.74335 | 751.31793 | 6.9580409 | 5.5036561 |

Boundary-corrected Pohozaev residuals decrease approximately as h squared. Raw residuals contain finite-box surface terms; their signs can cancel grid errors at an intermediate h. A small raw value at one mesh is therefore not an existence gate. Domain and grid changes remain material in E and P. Three-mesh extrapolations and their fit sensitivity are saved, but do not establish an infinite-domain branch. Negative background-subtracted E is retained, not reinterpreted as particle mass.

## 2. Continuation and first law

Parametric E(P), both charges, radii and all backgrounds are in the branch CSV/JSON files. The c and nu scans use separate outward paths from the same center. Local return solves quantify reproducibility; they do not constitute global fold detection or pseudo-arclength continuation.

| L | h | parameter | First-law residual, step→0 estimate | Graph envelope error, step→0 estimate |
|---:|---:|:---|---:|---:|
| 16 | 0.5 | c | 5.821098 | 0.0006605939 |
| 16 | 0.5 | nu | 0.5692164 | 4.082451e-07 |
| 24 | 0.333333 | c | 1.968959 | 8.489502e-05 |
| 24 | 0.333333 | nu | 0.2462237 | 1.377727e-07 |
| 24 | 0.4 | c | 2.809515 | 0.0002511187 |
| 24 | 0.4 | nu | 0.3521646 | 4.083184e-07 |
| 24 | 0.5 | c | 4.374268 | 0.0002509346 |
| 24 | 0.5 | nu | 0.54402 | 4.096899e-07 |

The discrete envelope error decreases with parameter-step refinement; the physical quadrature first-law error retains a grid-dependent floor. These are distinct tests. The derived finite-box c-law defect is -v R_Z/(2 Omega gamma squared); raw and corrected tests are both saved. The full continuum first law is not yet numerically closed. Dropping charge derivatives gives the wrong slope: raw dE/dP does not equal v on this varying-charge branch.

## 3. Characteristics and universal geometry

The full coupled Cartesian principal symbol is (eta contracted with k twice) times I_4. Potential couplings and rotating-frame gyroscopic terms are lower order; exact microscopic fronts share the already-assumed Minkowski cone. The Phi infrared phase has sound-speed squared 0.2. An admissible homogeneous two-condensate example has two different phase speeds squared, 0.04959603434 and 0.33918901239, independently recovered from exact coupled dispersion. The zero-carrier asymptotic state and ring cores invalidate extending that two-phase elimination everywhere.

**Negative geometry result:** an unqualified universal extension of the phase-only metric to all sectors fails these comparisons. Shared microscopic Minkowski characteristics are assumed structure, not emergent gravitational dynamics. This does not exclude every restricted infrared universality regime. Localized-loop coupling to a varying phase background has not been computed. No matter identification, universal gravitational coupling or Einstein equation follows.

## Reproduction and limits

See [derivation](three_track_derivation.md), [machine analysis](three_track_runs/analysis.json), [artifact hashes](three_track_runs/hashes.json), and [dispersion diagnostics](three_track_runs/characteristics.json). Profiles and failures remain separate research additions.

```text
python -m unittest test_full_z_spectral test_broad_scan test_three_track -v
python run_three_track.py --track convergence
python run_three_track.py --track continuation
python run_three_track.py --track continuation --fine-only
python characteristic_audit.py
python analyze_three_track.py
```

Use the unchanged requirements-spectral.txt. Future work should reduce h at larger L, control the cylinder-shape and far-field charge limits, and test the first law after both limits. No new stability scan was performed in this checkpoint.
