# Reproducible methodology

Jason D. Dutton | October 6, 2026

## Fixed protocol

The checked-in protocol was saved before the first simulation. It fixes sigma=0.2 as the primary comparison; seven noise levels, twenty seeds, calibration size, trial size, probe count, geometry assumptions, bootstrap seed, and stress shift are fixed. Timing requested batches are capped by the actual 4,000-trial dataset. No favorable seeds or noise levels were selected after seeing outcomes. The statistical analysis was expanded to include paired ablation contrasts and independent quadrature after the initial run; decoder definitions and primary trial data were unchanged. No external preregistration was performed.

For each seed and noise level, NumPy SeedSequence spawns separate streams for calibration injection, held-out displacement, probe measurement noise and randomized controls. Calibration uses 2,000 controlled, known displacements to estimate the population covariance; calibration data never enter held-out scoring. Twenty independent fitted priors account for calibration variation in the seed bootstrap. The same 4,000 held-out trials per seed are scored for every decoder/control. A new calibration set is supplied at each noise condition; blind adaptation to unknown noise is not tested.

`raw_trials.npz` contains truth_beta, observations_y, calibration_covariances, sigmas and measurement_covariance. Truth is evaluation-only. The sample mean is a sufficient statistic of the defined Gaussian heterodyne likelihood, so individual probe outcomes are not retained. Calibration sample covariances are retained; underlying calibration draws are regenerated exactly by the specified streams. No wavefunction measurement or noiseless syndrome is substituted for a probe observation.

## Rebuild and verification

Run the commands in README.md. `run_benchmark.py` produces trial data, every seed/aggregate metric, paired intervals, runtime samples, physical checks and a run record. `make_figures.py` reads saved CSVs; figures contain no manually substituted results. `build_report.py` reads the same results and the editable narrative. No online account is required after installing dependencies.

The exact run environment and protocol/dataset hashes are saved in results/run_record.json. A repeat must agree in recovered estimates and numerical metrics within floating-point tolerance. Runtime is system-dependent and will change; a compressed NPZ's byte hash can also depend on serialization versions/timestamps. Compare arrays and numerical CSV columns to verify scientific reproducibility rather than demanding an identical ZIP timestamp or runtime.

Nineteen stdlib-unittest tests passed in the recorded execution. They check coherent-state normalization and energy, independent ladder/displacement reference, inverse/unitarity, finite-cutoff versus exact overlaps, Kraus completeness and Choi partial trace, fidelity of an explicitly propagated Bell pair, six-state fidelity, leakage accounting, heterodyne variance, invalid inputs, ML/quadratic stationarity, independently coded posterior equivalence, coordinate controls, prior limits, truth access boundary, calibration, and paired-bootstrap determinism.

## Statistics

Seed means are independent clusters incorporating both calibration and evaluation variability. For each reported mean, sample 20 seed clusters with replacement 10,000 times; percentile 2.5 and 97.5% quantiles give a 95% interval. For contrasts, bootstrap the paired seed-level differences rather than subtracting two independent intervals. The primary positive-versus-ML criterion is an interval excluding zero at the predesignated noise level. Effect sizes include absolute fidelity difference and relative reduction in 1-Fe. Secondary contrasts have pointwise exploratory intervals, not adjusted global significance.

No hardware shots are simulated for final fidelity tomography: each residual's exact quantum overlap is scored. Probe shot noise is explicitly sampled. The Monte Carlo uncertainty concerns classical channel/probe draws and fitted prior variation, not an assumed experimentally estimated fidelity measurement. Independent 40-by-40 Gauss-Hermite integration checks expected entanglement fidelity of the linear decoders under the same Gaussian model.

## Costs and resources

Each trial uses one unknown payload oscillator and sixteen separate coherent probe modes with alpha=1, mean input probe energy sixteen photons. The maximally mixed encoded payload has input mean number 0.5. Channels share the same classical displacement during one trial; probes have independent heterodyne measurement noise. These are explicit resources for every principal decoder. Payload-probe displacement decorrelation is a critical untested physical risk; feature shuffling is only a synthetic proxy, not a physical decoherence model.

Covariance calibration costs O(Ncal) for fixed two-dimensional features. Geometry/posterior setup uses fixed-size 2-by-2 algebra; inference is O(N) for N trials, four multiply and two add operations per trial aside from array overhead. ML copies two numbers per trial. Exact scores are O(N). Finite-Fock cross-checks use eigendecomposition O(d^3), outside timed inference. Timings exclude calibration, quantum-state preparation, measurements, physical actuation, and fidelity scoring. They are amortized batched CPU throughput, not per-shot instrument latency or a quantum computational speedup.
