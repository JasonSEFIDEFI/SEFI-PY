# Microwave calibration and drift diagnostic workbench

**A field engineering demonstration: detect a change, identify what can be inferred, recommend the next measurement, and verify recovery.**

Built for a September 2026 discussion with eleQtron. This is an independent proposal based on public information, not an eleQtron project or endorsed implementation. All included measurements are synthetic.

## Review in five minutes

1. Read the [validation report](VALIDATION.md), including the negative result on algorithmic advantage.
2. Inspect the [service report](sample_output/service_report.md) and [ambiguous case](sample_output/ambiguous_report.md).
3. Download this folder and open `sample_output/index.html` in a browser for the visual demonstration. GitHub displays HTML source rather than hosting this page.
4. Inspect [the implementation](workbench.py), [tests](test_workbench.py), and [numerical results](sample_output/metrics.csv).

## Why this problem

eleQtron researchers publicly describe microwave-chain contributions to gate error and the need for coherent transfer between computational and storage transitions. The company's System Test Engineer role describes RF characterization, debugging, validation, and automation. Those provide a concrete motivation for a measurement-to-diagnosis workflow; they do not establish the company's current internal priorities.

- [eleQtron authors' ECCTI 2026 abstract](https://indico.fysik.su.se/event/9371/contributions/15175/)
- [eleQtron System Test Engineer role](https://eleqtron.com/en/jobs/rf-radio-frequency-test-engineer-f-m-d/)

## What is implemented

- Exact coherent population response for a driven two-level transition, initialized in its ground state.
- Twelve settings per scan: three signed drive-frequency offsets and four pulse durations.
- Binomial measurement counts with a fixed, known 2% symmetric readout error.
- Healthy, amplitude-drift, detuning-drift, and mixed-drift scenarios; faults ramp over scans 10–29 and are synthetically reset at scan 30.
- Independent nonlinear grid fitting; the same fits with exponentially weighted moving average (EWMA); and SEFI-inspired local tangent projection with the same EWMA.
- Method-specific alert thresholds calibrated on separate healthy runs.
- A measurement-only CSV replay interface, observability checks, model-mismatch flags, and a recovery report.

## Physical model and units

Let `a = 1 + amplitude_error`, `d = detuning / Omega0`, `c = configured frequency offset / Omega0`, and `t = physical pulse duration * Omega0`.

`H/(hbar * Omega0) = ((d-c)*sigma_z + a*sigma_x)/2`

`P_excited = a^2/(a^2+(d-c)^2) * sin^2(sqrt(a^2+(d-c)^2)*t/2)`

The measurement probability is `0.02 + 0.96*P_excited`. No nominal frequency in hertz is assumed. Every frequency must use angular-frequency units consistently; convert with `Omega0 = 2*pi*f0`. A scan is an ordered batch, not a specified number of seconds. Twelve settings at 1,024 shots each require 12,288 shots per scan; there is no claim of real-time operation or reduced measurement cost.

The reference is the ideal model at zero drift. Its uncertainty is **not** estimated from hardware. Drift remains fixed within each scan. These are synthetic engineering assumptions, not eleQtron specifications.

## How SEFI informs the proposal

The existing SEFI/DEFI research motivates representing normal behavior geometrically and resolving deviations relative to a reference. Here the reference geometry is defined in measurable population space, not asserted to be a new physical field.

For reference response `p0`, binomial standard deviations `s`, and the model's finite-difference Jacobian `J`:

`z = (observed - p0)/s`

`Jw = J/s`

`estimated deviation = pseudoinverse(Jw) @ z`

The resulting two coordinates estimate amplitude error and detuning. EWMA updates them with `state = 0.7*state + 0.3*estimate`, starting at zero. The orthogonal residual flags measurements that disagree with the local model. Local rank and conditioning checks identify an inadequate scan design.

**Weighted tangent projection, least squares, and EWMA are established methods.** This is a SEFI-inspired engineering adaptation, not a demonstrated novel algorithm or a derivation from the full SEFI field equations. The original geometric engine is not imported: its access to known injected noise is deliberately not used. Hidden simulation labels exist only in the evaluation layer.

All methods share the same measurements. The conventional smoothed baseline distinguishes any benefit of smoothing from a benefit attributable to the geometric representation.

## Run and reproduce

Python 3.10 or newer with NumPy; no SciPy, cloud account, instrument, or API key required. From this folder:

```sh
python -m pip install -r requirements.txt
python -W error -m unittest -v
python workbench.py --output demo_output --runs 40 --seed 20260925
python workbench.py --input demo_output/measurements.csv --output replay_output
```

Open `demo_output/index.html`. The bundled `sample_output` contains the default run. Reproduction was checked with Python 3.12.14 and NumPy 2.3.5; the recorded NumPy version is also in `metrics.json`. Other versions may vary in low-order numerical details.

CSV schema: `scan,setting,offset_omega0,duration_omega0,shots,excited`. Scan IDs start at zero and are contiguous. Each scan includes every setting once. Counts must be valid integers, and settings and shots must match the configured model. Duplicate, incomplete, and incompatible records are rejected. Replay still uses the synthetic reference and thresholds; importing a file does not make the tool hardware-calibrated.

## What the results establish

The default run evaluates 40 independently sampled runs per scenario, using separate random streams for calibration, evaluation, and the displayed example. Every method flags all 120 injected-fault runs at least once during their fault windows. Healthy per-scan false alarms are 1.11%, 0.94%, and 1.22% for independent fit, conventional EWMA, and tangent EWMA respectively. These are observed frequencies from a small synthetic benchmark, not guaranteed operating rates.

**No clear advantage over conventional tracking is established.** Median detection delays match across the three methods in these scenarios. Smoothing reduces healthy-state noise but lags at onset and recovery. In the mixed scenario the tangent tracker has slightly higher overall vector RMSE than conventional EWMA. The value demonstrated is the reproducible diagnostic workflow and explicit handling of uncertainty.

See [VALIDATION.md](VALIDATION.md) for exact values, metrics, limits, and checks.

## Field engineering workflow

1. Record the settings and measurements needed to reproduce an observation.
2. Estimate amplitude and detuning changes; check whether the model and scan support the inference.
3. If ambiguous, acquire positive and negative frequency-offset measurements or investigate reference/readout mismatch.
4. Cross-check with independent RF or resonance measurements before assigning a physical cause. A detuning shift alone cannot distinguish a magnetic change from a microwave-frequency change.
5. After an intervention, require three consecutive scans within the demonstration acceptance bands, accounting for approximate local uncertainty and residual checks. No correction is applied by this software.

Acceptance bands of +/-3% amplitude and +/-0.04 Omega0 detuning are illustrative. Recovery checks use independent fits; smoothed alarms may retain history after an apparent restoration.

## Boundaries and next pilot

This model does not include phase estimation, leakage, multi-level storage transfer, heating, motional dynamics, colored noise, preparation errors beyond the fixed readout abstraction, instrument latency, or correlated drift. Same-model synthetic generation and inference favor the estimators. Grid fitting has finite resolution and range. Local uncertainty estimates are not validated confidence intervals. The residual check can flag nonlinearity as well as a different physical fault. This is not gate-fidelity or logical-error-rate evidence.

A useful next pilot would identify one calibration workflow with an RF validation or operations engineer, obtain an approved anonymized scan format and operating limits, learn a measured reference, test model mismatch, and assess false alarms, diagnosis time, and measurement burden against the team's existing method. Retain whichever method performs best.
