# Physical matter research checkpoint

This opt-in research addition preserves the existing SEFI-PY engine. No engine imports, entry points, configuration, dependencies or QEC files are changed. Run it manually from this directory in a separate NumPy environment.

The explicit goal is to investigate whether a consistent field theory can derive physical matter and a geometry universally experienced by it. That goal is not achieved. The supplied model produces finite-domain moving charged-loop candidates, limited disturbance tests, and a rest-state scaling obstruction. Revised rest-loop searches are recorded as unsuccessful, not as solutions.

## Read first

- [Mathematical checkpoint](outputs/research_direction_2026_09_23.md)
- `work/charged_time_dt0.0125_amp0.001.json`: refined nonlinear test
- `outputs/charged_loop_dynamic_comparison.json`: phase alignment and time-step comparison
- `outputs/rest_branch_audit.json`: rest obstruction and proposed parameter window
- `work/rest_loop_lam12_nu085_N1.json`: converged to background
- `work/rest_loop_lam12_nu085_N3.json`: failed to converge

## Reproduction

Python 3 and NumPy 2 or newer are required. No dependency installation is performed by these scripts. The original run's versions and hashes are in `manifest.json`.

From this directory:

```text
python prepare_data.py
python work/rest_branch_audit.py
python work/solve_charged_ring.py --n 32 --h 0.5 --initial work/charged_ring_refined.npz --name checkpoint_smoke --steps 1
python work/charged_ring_time_probe.py --amplitude 0 --name charged_time_control
python work/charged_ring_time_probe.py --amplitude 0.001 --dt 0.025
python work/charged_ring_time_probe.py --amplitude 0.001 --dt 0.0125
python work/charged_ring_azimuthal_probe.py --m 1
python work/charged_ring_azimuthal_probe.py --m 2
```

Reference profiles, including failed-search outputs, are stored in `reference_profiles.zip`. `prepare_data.py` extracts only inside this research directory and refuses to overwrite different data. Profile and evolution runs write to this directory's `work` folder; named reruns can replace research outputs there. Preserve the archive and JSON histories when exploring new parameters.

The profile solver supports `--lam`, `--winding`, `--nu`, and `--c` for explicitly labeled exploratory changes. The time-evolution probes currently implement only the baseline lambda=25, zero carrier-azimuthal-winding model. They reject other parameters and must not be used to claim stability of the lambda=12 searches.

These are finite-difference research solvers. Small residuals and conserved charge are numerical checks, not evidence of quantum matter, electric charge, universal gravity, or global stability. Full three-dimensional nonlinear stability, continuum and domain limits, physical calibration, and independent validation remain open.

## Changes and validation

Added coupled profile and time-evolution solvers, azimuthal linear probes, the algebraic rest-state audit, reference data and a research-direction note. Existing engine behavior is preserved by isolation. The checkpoint includes an unperturbed control, time-step refinement, charge monitoring, phase-aligned diagnostics, and two explicitly unsuccessful rest-loop searches. No new physical result is inferred from software test success alone.
