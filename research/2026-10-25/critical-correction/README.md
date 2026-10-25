# First correction to the stationary critical limit

25 October 2026. **Exploratory / Not Promoted. Working manuscript draft only.**

This extends the [frequency and critical-limit study](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/446e6177ec30b840a0fcd336ff72dd96a1c2df4e/research/2026-10-25/localized-frequency). The first-order coefficients below come from the reduced equations and a transverse-amplitude calculation. No full-model charge data are fitted.

## Main result

For delta=12/lambda-3/sqrt(lambda), the formal predictions are

    Q_turn delta^(3/2) = 38.53902853 + 66.25169133 delta + O(delta^2),
    Q_binding delta^(3/2) = 50.59921950 + 76.66926995 delta + O(delta^2).

The remainder notation is a formal asymptotic claim under the smooth-branch assumptions in [DERIVATION.md](DERIVATION.md), not a proved uniform estimate.

Coefficients were frozen in predictions.json before the declared new full-model couplings were computed. At those couplings:

| lambda | Quantity | Leading relative error | First-correction relative error |
|---|---|---:|---:|
| 15.95 | Charge minimum | 0.20167% | 1.7013 ppm |
| 15.975 | Charge minimum | 0.10078% | 0.4246 ppm |
| 15.95 | Binding crossing charge | 0.17780% | 1.4948 ppm |
| 15.975 | Binding crossing charge | 0.08884% | 0.3728 ppm |

One ppm is a relative error of 1e-6. These are errors against computed full-model stationary branches, not experimental errors or certified continuum bounds.

![Correction evidence](correction_evidence.png)

The scaled charge remainder divided by delta squared is approximately -47.47 and -47.52 for the turn, and -54.74 and -54.78 for binding. Two points are consistent with the next asymptotic order, not proof of its order on an interval.

## Checks completed

- Exact symbolic expansion of the local functional and transverse completion of the square.
- Direct evaluation of the corrected fields in the full finite-delta functional, with remainders consistent with second order.
- Four differentiation step sizes for the on-shell action derivative.
- An independent linear profile-sensitivity BVP reproduces the turn coefficient within about 6.2e-11.
- Secondary position checks, computed after the charge tests, give beta_turn=beta_turn,0+0.01183051 delta+O(delta^2) and beta_binding=beta_binding,0+0.02868270 delta+O(delta^2). Their residuals are also consistent with second order at the two points; these were not part of the frozen protocol.
- The leading envelope identity is checked numerically.
- Full-model domain X=40 to 48 and tolerance 1e-9 to 1e-10 changes the scaled charges by less than 4e-12 relatively, well below the reported prediction errors.
- The frozen prediction file and the source hashes used to generate it are checked before producing the comparison.

The first reduced binding solve attempted a direct jump from the turn seed and converged to the trivial solution; it was rejected before predictions.json was created. Choosing the already available binding seed corrected initialization without changing the equations. No rejected solution is included as evidence.

## Draft and review material

- [MANUSCRIPT_DRAFT.md](MANUSCRIPT_DRAFT.md): an integrated working paper focused on the stationary reduction.
- [DERIVATION.md](DERIVATION.md): correction functional, charge formulas and conditional assumptions.
- [AMPLITUDE_BOUND.md](AMPLITUDE_BOUND.md): an a priori transverse-amplitude error estimate for prescribed f under explicit positivity and regularity assumptions; it does not prove the full coupled expansion.
- [LITERATURE_COMPARISON.md](LITERATURE_COMPARISON.md): prior-work overlap and remaining novelty questions.
- [heldout_comparison.csv](heldout_comparison.csv): machine-readable errors.
- predictions.json: frozen coefficients, numerical derivative checks, source hashes and declared test couplings.
- heldout/results.json and heldout_validation/results.json: full-model solves and refinements.
- sensitivity_results.json and expansion_validation.json: independent mathematical checks.
- NPZ files: seeds, dense correction profile and final full-model profiles.
- manifest.json: file checksums.

## Reproduction

Python 3.14.7 with requirements.txt. Use a fresh copy of source and the three seed NPZ files. Existing output files/directories are protected.

    python -m pip install -r requirements.txt
    python derive_correction.py
    python heldout_test.py
    python validate_expansion.py
    python validate_heldout.py
    python sensitivity_check.py
    python position_checks.py
    python check_amplitude_bound.py
    python summarize_correction.py

The reduced solver is the same critical-curve BVP from the preceding study, with initialization from its saved turn/binding profiles. full_solver.py retains the preceding full-model equations and collocation method, omitting only its old command-line scan. The full-model seeds are the prior lambda=15.9 turn and binding profiles. The method solves the full second-order field equations, not a nonperturbative fourth-order truncation.

The word held-out means these new full-model charges were not used to determine the coefficients. This is an adaptively designed in-session check, not an externally preregistered or independently blinded study.

## Remaining limits

The static chart requires cos(theta)>0; no uniform endpoint/layer analysis is supplied. Smooth dependence and existence of the full branch, rigorous error bounds, global energetic optimality and nonlinear stability remain separate questions. The literature comparison is incomplete. Substantive AI assistance was used throughout; author verification and independent technical review remain required.

No result is approved for promotion and no manuscript is submitted by this record.

