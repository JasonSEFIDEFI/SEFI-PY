# Full-model fold session — 26 October 2026

This session extends the Oct 25 critical-binding audit. It uses the original dimensionless two-field equations; no physical-matter or novelty claim is made.

## Independently derived prediction

With delta=16-lambda,

Q_min = delta^(-3/2) [q0 + delta*q1 + O(delta²)],

q0 = 10740.73088927, q1 = -825.92568465.

The fold frequency parameter satisfies alpha_fold=alpha0+delta*alpha1+O(delta²), with alpha0=0.002914930803, alpha1=0.000234227196.
These coefficients were evaluated from reduced solutions before the new full-model sweep. Finite-difference derivative refinements are saved; this is not external preregistration.

## Observed full-model folds

| delta | Scaled minimum charge | Q_min | alpha_fold | Corrected scaled-charge error |
|---|---:|---:|---:|---:|
| 0.500000 | 10330.801308 | 29219.919 | 0.0030356901 | 3.033261 |
| 0.250000 | 10535.009752 | 84280.078 | 0.0029743834 | 0.760283 |
| 0.125000 | 10637.680529 | 240703.233 | 0.0029444312 | 0.190351 |
| 0.062500 | 10689.158168 | 684106.123 | 0.0029296252 | 0.047634 |
| 0.031250 | 10714.932630 | 1939609.990 | 0.0029222646 | 0.011919 |
| 0.015625 | 10727.828784 | 5492648.337 | 0.0029185945 | 0.002984 |

Adjacent-pair effective exponents are 1.528239, 1.513992, 1.506965, 1.503475, 1.501735, approaching the predicted 1.5.
The corrected remainder orders are 1.9963, 1.9979, 1.9986, 1.9987, 1.9981, approximately second order.
At the closest delta the charge-prediction error falls from 0.120268% to 0.000028%, a factor of 4324.3.
The longer-domain/tighter-tolerance fold repeat changes scaled charge relatively by 2.89e-10. This is numerical sensitivity, not a rigorous error bound or statistical interval.

The charge minimum is a local fold of the selected positive monotone spherical branch, not a global minimum over all fields. Its scaled binding remains positive at every tested coupling, so reaching the minimum charge does not by itself give energetic binding below free particles.

## Radial curvature

Two couplings, three central amplitudes each, and two grids produce 12 lowest-radial-eigenvalue blocks. Below each fold the lowest eigenvalue is negative; above it it is positive.
At the fold, small negative grid eigenvalues shrink toward zero; their coarse/fine ratios are 4.0049, 4.0069.
The corresponding eigenvectors overlap the stationary-family tangents by 0.999999944, 0.999999944.
The exact identity H_Q F_c = (2 omega Q_c / I)(f,0) explains why the charge fold has a radial zero mode. The radial matrix uses scaled coordinates, so its eigenvalues are physical curvature divided by delta.

The first family-identity residual used a shorter spectral domain (240) than the stationary boundary-value domain (320); it worsened on one refinement because the unmodified tangent did not obey the shorter-domain outer boundary. Those diagnostics remain in radial_results.json. A new linear differentiated-equation solve on the matching domain (320) supplies the tangent without finite parameter differences. One derivative-based Newton correction reduces absolute charge derivatives to 1.16e-09, 3.08e-09. On three grids (.5, .25, .125), matching-domain identity residuals decrease roughly fourfold with each halving, reaching 3.55e-06, 3.53e-06. This supports the zero-mode identity while making the original boundary mismatch explicit.

These results support a simple local quadratic-energy index crossing. They are not a charge-conserving dynamical simulation or a nonlinear stability theorem.

## Verification and limitations

Three new mathematical tests pass. The fold search used 83 boundary-value solves including the domain/tolerance repeat; the radial calculation records its separate profile solves. ODE residuals, off-node derivative checks, virial residuals, monotonicity, frequency bounds, timings and optimization diagnostics are retained.

An initial empty-parameter seed interface error stopped the first sweep before any full fold was evaluated. Its log is retained. Corrected initialization uses a reduced profile with its frequency parameter.

The first-correction coefficient has finite parameter-difference error; do not interpret all printed digits as guaranteed precision. Only one full-fold domain refinement and two eigenvalue grids per sampled profile were run; the matching-domain identity uses three grids but is not an additional eigenvalue computation. Extrapolation to either frequency endpoint, other branches, gauged fields, complex Sigma, quantum behavior or measured nuclei is unsupported.

This strengthens the controlled mathematical description of the model's binding transition. The next discriminating step is charge-conserving time evolution near the radial fold, followed by comparison of these specific formulas with prior quartic-extended soliton literature. The separate foundational worldline-to-field and field-to-observation bridges remain open.
