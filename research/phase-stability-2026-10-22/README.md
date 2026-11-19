# Phase stability: a bounded correction/degradation model

Jason Duran Dutton, SEFI Research Program | October 22, 2026

Start with **RESEARCH_NOTE.md** for the source audit, explicit new assumptions, dimensional analysis, derived stationary-lock criterion and a possible calibration/experiment protocol. This is classical phase tracking, not quantum error correction or a universal fatigue theory.

## Main result

For gain K(d)=K0/(1+d) and recoverable activity accumulation d_dot=alpha u²-beta d, the principal stationary locked state is locally asymptotically stable if

    |delta| (1 + alpha delta²/beta) < K0.

The note derives this condition from the equations and checks the Jacobian. The added gain and accumulation laws are phenomenological modeling choices, not uniquely implied by the original SEFI source. Establishing actual relevance requires measuring them independently.

With K0=2/s, beta=1/s and alpha=1 s, the stationary-lock threshold is 1 rad/s; the constant-gain limit is 2 rad/s. At detuning 1.3 rad/s, the baseline locks but the degraded-gain trajectory drifts. This predicts a failure mechanism under the model, not better control or new physics.

## Numerical checks

- Five scientific tests passed: cubic threshold, stationary-state equations/stability, independently finite-differenced Jacobian, constant-gain limit and marginal threshold.
- All 36 prescribed trajectories completed: two models, six detunings and three step sizes.
- Every finite-time lock classification agrees across step sizes. Maximum end-phase difference relative to the finest step is 1.49e-7 rad.
- For every test with a strict stationary equilibrium, the trajectory from the specified zero initial state reached the analytical equilibrium within 1e-7. This does not prove global attraction from other initial states.

## Reproduce

Python 3.12+ standard library is sufficient for the model and tests:

```text
python -m unittest -v test_phase_model
python run_study.py
```

To regenerate the scientific SVG figure, install ReportLab 4.4.9 and run `python build_figure.py`. The PNG is a raster rendering of the supplied SVG. Reproduction overwrites the corresponding results; preserve the original package for comparison. A started/completed/interrupted/failed manifest records study status; incomplete stages are not final evidence.

## Boundaries

The local source `Rough/1.tex` is also a transcription; ambiguous original intent is unresolved. No source manuscript or main research-console equation was changed. Assertions of universal biological, gravitational or cognitive equivalence are not carried into this model. Conventional actuator-state models may reproduce the same result exactly. No unique advantage, novelty priority, hardware result or journal readiness is asserted.

The next useful external input is a laboratory's assessment of whether its actuator shows the specified gain-versus-load and recovery behavior. Parameters must be measured on calibration data before predicting held-out phase failures. The author has not commissioned consultancy or contacted anyone through this task.
