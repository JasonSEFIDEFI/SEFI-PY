# Model history and development contract

## 0.1 — Numerical checkpoint console, 23 September 2026

Implements the saved coupled traveling-ring equations and an axisymmetric RK4 evolution for carrier self-coupling λ=25 and carrier azimuthal winding N=0. Other allowed λ and N values are nonlinear profile searches only. The saved refined, coarse, and wider-domain profiles remain immutable comparison evidence in `data/`.

The user requests retaining this console and incorporating subsequent mathematical findings. This means a versioned research instrument, not automatically treating a new conjecture as a validated law.

For each future model change:

1. Record the derivation, definitions, units, approximation domain, and which earlier model version it extends or supersedes.
2. Implement the change in an explicit model version. Preserve the baseline equations and original reference data so old results remain reproducible.
3. Add discriminating tests: limiting cases, conservation or balance laws, residuals, convergence, and an appropriate failed or negative control.
4. Update controls, mathematical explanations, reviewer questions, and claim status together. Keep unimplemented proposals visibly separate from executable models.
5. Write model version, parameters, numerical scheme, and source identifiers into saved runs and checkpoints. Reject incompatible checkpoints rather than silently reinterpret them.
6. Compare old and new predictions. State what changed, what evidence supports it, and what remains unresolved.
7. Review the working local console before publishing the update. Preserve the existing SEFI-PY engine and SEFI-QEC.

No current model version is a derivation of observed physical matter or universally experienced spacetime geometry.
