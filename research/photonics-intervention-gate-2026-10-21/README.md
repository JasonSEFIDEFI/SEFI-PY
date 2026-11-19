# Photonic correction intervention gate — feasibility study

Jason Duran Dutton, SEFI Research Program | October 21, 2026

## Result

We derived and simulated a conditional correction gate that compares expected physical fidelity with and without actuator disturbance. This is conventional Bayesian decision theory, not a new SEFI law. It supports an engineering effect relative to always applying a noisy correction; it does not establish superiority to a strong conventional gate.

| Actuator variance per real component | Always bypass | Always correct | Conventional MSE gate | Fidelity gate | Fidelity-gate correction fraction |
|---|---:|---:|---:|---:|---:|
| 0 | 0.947291 | 0.965824 | 0.965824 | 0.965824 | 100% |
| 0.0025 | 0.947291 | 0.956465 | 0.959160 | 0.959161 | 51.29% |
| 0.01 | 0.947291 | 0.929195 | 0.951366 | 0.951361 | 15.47% |
| 0.04 | 0.947291 | 0.833747 | 0.947468 | 0.947476 | 0.49% |

At variance 0.01, the gate-minus-always-correct fidelity difference is 0.022166 (2.217 percentage points), descriptive seed-bootstrap 95% interval 0.021722 to 0.022598. The fidelity gate versus conventional MSE gate difference is -0.00000444, with interval -0.00002670 to 0.00001789. None of the four settings shows a resolved fidelity advantage over MSE gating in these exploratory intervals. No favorable setting is a prespecified confirmatory endpoint. At very large correction disturbance, bypass alone explains most of the apparent advantage over always correcting.

All twenty seeds and all four actuator settings completed: 80,000 independent displacement trials, reused across the specified actuator settings. Exact Gaussian expected-fidelity calculations match independent quadrature within 1.12e-16. The conditional gate selects the higher expected-fidelity action under the assumed calibrated model. No real-time timing or hardware test was performed.

## Reproduce

Use Python 3.12+ with NumPy 2.3.5:

```text
python -m pip install -r requirements.txt
python -m unittest -v test_gate
python pilot.py
```

The pilot overwrites its results. Preserve the supplied results before reproducing. `protocol.json` defines all settings and assumptions. `DERIVATION.md` supplies the formulas and scope. `results/simulation_inputs.npz` retains the observations, hidden evaluator truth and actuator realizations; raw audio is not part of this project. The run manifest marks incomplete/interrupted stages.

## What a laboratory could test

The rule is relevant if a laboratory already has calibrated auxiliary sensing, a displacement actuator and a bypass/control option. Compare always-correct, always-bypass, MSE-gated and fidelity-gated operation, counting detector inefficiency, actuator error and bypass loss. A coherent-beam demonstration verifies classical control only. A quantum-retention claim requires a valid quantum-state source, measurement protocol and unconditional scoring. Pilot/payload common noise must be experimentally justified.

The useful question for a potential collaborator is whether actuator-induced disturbance is a significant resource cost in their apparatus. This package does not establish that it is, and is not a reason to purchase apparatus before checking that question.

## Relation to SEFI

This formalizes a possible resource-aware version of the manuscript's restoration idea. It does not operationalize curvature/torsion as a syndrome or prove that geometric projection is a physical quantum channel. It offers a concrete engineering test, not validation of the broad theory. Novelty beyond conventional fidelity-aware control has not been established. A separate source audit records the missing theory-to-observable bridge.

The existing SEFI-PY noncommercial research license applies. No license, submitted manuscript or main research-console equation has been changed. This is an exploratory follow-up to the adaptive-recovery benchmark, not an independent validation of it.
