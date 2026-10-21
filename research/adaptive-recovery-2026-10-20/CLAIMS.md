# Supported claim and exclusions

In a frozen synthetic pilot-assisted displacement-recovery benchmark, a conventional adaptive mixture improved mean entanglement fidelity over a development-selected 128-sample rolling estimator by 0.210 percentage points in the prespecified drifting-noise scenario (paired seed-bootstrap 95% interval 0.177 to 0.244 points; 20 independent streams). It did not outperform the tested EW estimator in that scenario and was slower in the measured Python implementation. This does not establish a SEFI-specific advantage.

## Do not claim

- A new SEFI algorithm, theory validation, quantum advantage, or superiority to all conventional methods.
- Hardware QEC, a protecting-code threshold, measured logical error reduction, or real-time feedback readiness.
- Faster or better at matched computation cost: the candidate uses more measured computation.
- Generalization beyond the five scripted synthetic noise families.
- Independent preregistration or independent reproduction; the protocol freeze and validation were local.

## Useful engineering observation

The EW baseline achieved drift fidelity 0.949550, versus 0.946370 for the rolling baseline. Its median implementation time was 15.37 versus 75.57 microseconds per observation (4.92 times lower). This is a secondary descriptive comparison of implementations, not a new estimator or SEFI-specific speedup.
