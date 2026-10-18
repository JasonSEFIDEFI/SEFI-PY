# Negative controls and limits on interpretation

Jason D. Dutton | October 6, 2026

## Equal-information control

ML, geometry and conventional posterior have access to the same two-dimensional sufficient statistic from sixteen probe heterodyne outcomes, the same measurement variance, and the same independent calibration covariance. The classical simulator's held-out beta is passed to scoring only. No decoder receives the payload state, its density matrix, trial displacement truth, future observations, or a per-trial noise label. The ML rule ignores the available prior because maximum likelihood optimizes the likelihood rather than a posterior. Equality of available information does not make these inference objectives identical; this is why the conventional posterior is essential.

The geometry and independently coded conventional posterior estimates differ by at most 1.78e-15 over all held-out trials. Thus any gain over ML does not survive attribution against the stronger ordinary Bayesian control. This is a null result for a geometry-specific advantage within the tested family.

## Coordinate reparameterization

The affine control uses T = [[1.7, 0.4], [-0.3, 0.65]], transforms both covariance matrices as T Sigma T^T and T R T^T, predicts in that chart, and transforms the correction back. Its maximum estimate discrepancy is 1.78e-15. Metric transformation is essential; applying unchanged Euclidean distances after shear would define a different estimator.

The nonlinear control encodes observations as arcsinh(y), inverts them by sinh, and evaluates the physical-coordinate estimator. Its discrepancy is also 1.78e-15. This is a bijective-information/implementation check, not a hard test of a decoder optimized intrinsically in a curved chart. It does not establish invariant MAP density: density maximizers can change with a coordinate Jacobian. An intrinsic nonlinear manifold estimator and robust chart-transition tests remain future work.

## Feature ablation

The geometry score has two feature groups: distance to the zero-displacement reference under an empirical prior metric, and that metric's calibrated anisotropy/orientation. Removing prior distance gives ML and eliminates its improvement. Replacing the anisotropic metric with trace(Sigma)/2 times I keeps isotropic shrinkage, but reduces primary-condition fidelity by about 0.479 percentage points versus full geometry. There is no recurrence or phase-history feature in this memoryless model, so there is no such feature to ablate.

## Randomized geometric features

Randomizing metric orientation with one independent uniform angle per seed preserves eigenvalues and shrinkage strength. It reduces fidelity by about 0.875 percentage points at the primary condition, but retains much of the gain over ML. This is evidence that ordinary shrinkage explains most of the benefit; a claim that random geometry eliminates every improvement would be false.

Independently permuting probe features across trials destroys their alignment with the payload displacement while preserving marginal feature statistics. Fidelity falls below ML and below full geometry. This is a deliberately invalid correction strategy, not a competitive decoder. It is an information-alignment sanity check, not proof of algorithmic novelty.

## No-recovery control

No recovery avoids acting on measurement noise. It exceeds ML at low displacement levels, including the primary condition, and must appear in every performance figure/table. ML is not uniformly preferable to doing nothing. A geometry-versus-ML comparison alone would overstate the importance of active recovery.

## Distribution-shift negative test

The independent test shift beta -> beta + (0.75, -0.40) is applied to both payload displacement and probe observations while leaving the zero-mean calibration prior unchanged. At sigma=0.2, geometry fidelity falls from 0.95212 to 0.49715; ML remains 0.88792. The paired geometry-minus-ML effect is -0.39077, with 95% seed-bootstrap interval [-0.39375, -0.38776]. This is a synthetic prior-mismatch stress test, not a calibrated model of drift frequency in an instrument.

## Physical and numerical negative checks

Leakage is never conditioned away. The decoded channel retains leaked mass in an orthogonal erasure output. Analytic oscillator scores are compared with independent finite-Fock exponentiation at cutoffs 8, 16, 32 and 48. CPTP Kraus completeness and Choi positivity are checked. The no-displacement state has unit target fidelity; nonzero correction of that state can lower fidelity. A radius-3 matrix grid does not cover all Gaussian tails, but headline scores are exact infinite-oscillator overlaps for every sampled residual, including tails.

All secondary ablation comparisons are exploratory. The reported pointwise intervals are not multiplicity-adjusted; they do not justify a familywise claim about every control. The primary endpoint was designated before the first simulation in protocol.json.
