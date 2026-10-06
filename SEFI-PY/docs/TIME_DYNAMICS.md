# Time manuscript mathematics

Author: **Jason D. Dutton**. Source: *Phase Dynamics, Recurrence, and Relational Temporal Ordering in SEFI-DEFI-GWFM*, submitted manuscript (2026), Sections 2–8 and 13. Submission does not imply acceptance or physical validation.

`core/time_dynamics.py` implements the explicitly specified classical identity reduction. It is independent of the existing vector-blending `WarpDEFI` layer. All quantities are dimensionless; curvature and torsion are prescribed rather than varied. It does not implement the full GWFM spacetime variational problem.

## Implemented results

| Manuscript result | Interface | Conditions and limits |
|---|---|---|
| Positive diagonal metric, base identity equations, energy | `IdentityMetric`, `Coupling`, `IdentityModel.rhs`, `.energy`, `.prescribed_power` | Supply C2 coupling functions and their correct derivatives. Total energy is conserved for constant geometry; changing geometry can supply external work. |
| Origin energy and directed phase | `origin_energy`, `clock_phase`, `phase_history` | Zero amplitude has no phase. Unwrapped history requires a cycle record. Sampled records must match the specified harmonic oscillator and have gaps smaller than half a cycle. |
| Relational chart | `relational_derivative` | Positive phase speed is required. This reparameterizes a solution; it does not remove the evolution parameter from the action. |
| Constant linear spectrum and trajectories | `linear_parameters`, `linear_spectrum`, `linear_modes`, `linear_solution` | Every solution is bounded below r=1. At r=1 drift is possible. Above r=1 growing solutions exist; not every initial state grows. |
| Effective potential and stationary Hessian | `effective_potential`, `effective_derivatives`, `equilibrium`, `equilibrium_spectrum` | Completion of squares is algebraic, not dynamic coordinate elimination. Positive stationary Hessian gives local stability; negative eigenvalues imply instability; zero needs higher-order analysis. |
| Exact linear recurrence | `exact_linear_ratios`, `recurrence_cycles` | Exact rational certificates only. Include excited modes. Finite-precision near returns do not prove exact recurrence. |
| Calibrated lift locking and cycle descent | `locking_residual`, `periodic_descent_residual` | Locking is kinematic, not an attracting synchronization mechanism. A sampled descent residual is a diagnostic, not proof for all phases. |
| Driven soft-mode instability | `floquet_monodromy`, `floquet_diagnostics` | Fixed coupling direction and sinusoidal strength. Repeated multipliers require Jordan analysis; numerical determinant errors are reported as unresolved. |
| DNA norm penalty and helical mismatch | `dna_origin_energies`, `dna_ratio_mismatch`, `dna_linear_spectrum` | Origin is included in the penalty. Base energy is generally not conserved; modified energy is conserved for constant mismatch. The shifted spectrum applies to constant mismatch and geometry only. |

The manuscript's coercivity result is a sufficient global theorem, not a finite numerical test. No generic routine claims to prove coercivity for arbitrary supplied functions. Experimental observability, proper-time calibration, quantum error-correction conditions, and the thermodynamic arrow remain outside this implementation.

## Reproduce the manuscript examples

From the `SEFI-PY` directory:

```console
python demo_time_dynamics.py
python -m pytest tests/test_time_dynamics.py
```

The module and demo require only the Python standard library. Tests require pytest; other pre-existing engine modules have their own dependencies.

The resonant example has r(tau)=0.5+0.1*cos(sqrt(2)*tau), so instantaneous stiffness stays positive while the Floquet growth rate is approximately 0.03531402295. The comparison drive at 2.2 has distinct unit-modulus multipliers. Clock conservation alone therefore does not guarantee system stability.

At r=24/25 the excited linear frequency ratios are exactly 1/5, 1, and 7/5; five Origin cycles yield full-state recurrence at 10*pi. At r=1/2 the soft frequency is irrational relative to the unit Origin clock, and no finite exact full-state recurrence follows when both modes are excited.

For the nonlinear example mu1(IA)=IA^2/2, mu2=0, unit metric, kappa=1, sigma=0, Veff=IA^2/2-IA^4/8. The center is locally stable, saddles occur at +/-sqrt(2), and the barrier is 1/2. This potential is not globally coercive. Rounded DNA geometry gives proximity to pi/phi, not exact locking or a topological invariant.
