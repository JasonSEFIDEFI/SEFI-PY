# Architecture additions, suggested commits, and pull request description

## Repository structure and integration

All additions live under `SEFI-PY/benchmarks/displacement_recovery/`. No existing import namespace, license, classical core or simulation console is changed. Dependencies for this directory are isolated in its requirements.txt. Local unit tests use unittest, so pytest is not required to run the new physical checks. A future CI job can install this directory's requirements, run the tests and regenerate aggregate values; archive raw NPZ/figures separately if CI size limits matter.

The README lists the complete file structure. `bosonic.py` is the minimal quantum architecture: ladder/coherent state, physical displacement, exact overlaps and erasure-channel representation. `decoders.py` accepts observations and calibration covariance only. `run_benchmark.py` alone owns held-out simulator truth. The metrics layer accepts residuals after prediction. This information boundary is intentionally explicit.

## Suggested commits

1. `audit: inventory pinned SEFI-PY quantum recovery gaps` - technical audit, inventory and static source diagnostics.
2. `feat: add pilot-assisted bosonic displacement recovery reference` - fixed protocol, oscillator channel/metrics, ML and geometry/posterior implementations.
3. `test: verify physical channels and matched-information controls` - CPTP/Bell-pair/six-state checks, coordinate controls and prior limits.
4. `bench: record seeded recovery results and negative controls` - raw and aggregate data, confidence intervals, mismatch stress test and timings.
5. `docs: add QIS benchmark figures and calibrated conclusions` - report, three figures, reproduction instructions and NSF preliminary-work wording.

Stage only this new benchmark directory. Do not inadvertently include existing uncommitted Time-module changes or other local research. Suggestions do not imply that commits have already been created.

## Pull request title

Add reproducible pilot-assisted bosonic recovery benchmark with matched Bayesian control

## Pull request description

The existing repository did not provide a validated bosonic displacement-recovery benchmark with measured quantum fidelity and matched-information decoder controls. This addition simulates an unknown qubit encoded in the first two Fock levels, uses sixteen separate coherent probes to observe a shared displacement, and applies measurement-only recovery. It adds exact unconditional entanglement fidelity, leakage accounting, a maximum-likelihood baseline, geometry-distance regularization, an independently implemented conventional Gaussian posterior, coordinate transformations, feature ablations and randomized negative controls.

Across 560,000 held-out trials and twenty calibration/evaluation seeds, the primary-condition geometry estimator improves Fe from 0.88792 to 0.95212 versus ML; the paired effect is 0.06419 with 95% seed-bootstrap CI [0.06358, 0.06476]. Geometry and conventional posterior predictions agree to 1.78e-15. An unmodeled mean shift reverses the gain. The result supports a reproducible conventional-physics recovery workflow, not a novel geometric decoder advantage, a protecting code, a fault-tolerance threshold, hardware performance or validation of SEFI/DEFI/GWFM.

Validation: nineteen physical/inference tests pass. Independent Bell-pair propagation and six-state averages match analytic fidelity; Kraus/Choi checks establish CPTP erasure decoding; cutoffs 8/16/32/48 converge to exact oscillator overlaps. Every owned source in the pinned repository inventory was read and hash matched; five pre-existing syntax-invalid QEC Stack files are documented and excluded. The new directory does not modify legacy core code or the existing research license.

Review artifacts: editable benchmark report and PDF, three PNG/SVG figures, full raw NPZ and seed/aggregate CSV results, independent Gaussian-quadrature check, protocol/environment hashes, timing samples, negative-control discussion, and suggested NSF preliminary-work paragraph. Local results are computational simulations only.
