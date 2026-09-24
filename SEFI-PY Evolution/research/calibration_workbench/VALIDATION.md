# Validation record

Executed locally on 24 September 2026 UTC. Data and hardware settings are synthetic. The default experiment uses seed `20260925`, NumPy 2.3.5, 1,024 shots per setting, 12 settings per scan, 45 scans per run, 40 healthy calibration runs, and 40 evaluation runs for each of four scenarios.

## Held-out results

| Method | Healthy false alarms / scans | Healthy vector RMSE | Mixed vector RMSE | Median mixed detection delay | Missed fault runs across 3 scenarios |
|---|---:|---:|---:|---:|---:|
| Independent fit | 20 / 1,800 (1.11%) | 0.00890 | 0.00913 | 1 scan | 0 / 120 |
| Conventional EWMA | 17 / 1,800 (0.94%) | 0.00369 | 0.02170 | 1 scan | 0 / 120 |
| SEFI-inspired tangent EWMA | 22 / 1,800 (1.22%) | 0.00362 | 0.02211 | 1 scan | 0 / 120 |

Amplitude-only and detuning-only median detection delays are 2 scans for every method. After the synthetic reset in the mixed scenario, median time to the first clear alarm is 0 scans for independent fits and 8 scans for both smoothed methods. Alarm clearance is not the report's recovery acceptance test.

The tangent method is not shown to beat conventional EWMA. Differences in false-alarm rate are small, and temporal dependence precludes treating smoothed scans as independent Bernoulli trials. No significance or generalization claim is made.

## Metric definitions

- Vector RMSE: square root of mean squared Euclidean error in `[fractional amplitude error, detuning/Omega0]`, across all scans and evaluation runs. It combines unlike but dimensionless quantities without a hardware cost weighting.
- Healthy false alarm: a score above the separately calibrated method threshold in a held-out healthy scan. Threshold target is the empirical 99th percentile on calibration data; this is a per-scan target, not a per-run guarantee.
- Detection delay: first alert in scans 10–29, minus 10. A value of zero means detection in the first affected scan. Medians exclude missed runs; missed-run counts are reported separately.
- First-clear delay: first below-threshold score after scan 30, minus 30. This is not sustained clearance or proof of recovery.
- Alert scores are reported even when the service workflow abstains on cause attribution. Detection does not mean correct component diagnosis.

## Executed checks

Eleven automated tests passed with runtime warnings treated as errors:

1. Resonant pi-pulse, 2pi-pulse, and zero-drive physical limits.
2. Probability bounds across broad drive and detuning values.
3. Signed detuning and amplitude recovery from noiseless off-resonant scans.
4. Resonant-only sign ambiguity and explicit abstention.
5. Tangent projection agreement for a small physical perturbation.
6. Three-scan recovery acceptance and rejection of an unrecovered fault.
7. Model-mismatch ambiguity and withholding recovery acceptance.
8. Exact CSV replay and separation of measurement columns from simulation truth.
9. Rejection of duplicate, incomplete, and invalid-count CSV input.
10. Rejection of invalid probabilities and configuration.
11. Seed reproducibility and independence of different generated samples.

The complete default benchmark and CSV replay were also executed with warnings treated as errors. Source, observations, evaluation-only truth, metrics, and reports are included. These tests cover this isolated addition, not every pre-existing repository module.

## Important limitations

The generator and estimators use the same physical response and a perfectly known reference/readout model. This is an implementation and workflow test, not an independent physics validation. The geometric projection is local; grid fitting is quantized to 0.005 in both coordinates and bounded to +/-0.20 amplitude error and +/-0.30 detuning/Omega0. At boundaries the report abstains. Local approximate two-sigma bounds omit reference uncertainty and need empirical coverage validation. Real data could expose missing dynamics and invalidate the assumed operating region.
