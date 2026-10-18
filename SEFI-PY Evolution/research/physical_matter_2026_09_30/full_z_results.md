# SEFI-PY full-Z spectral checkpoint — 30 September 2026

Built and ran an independent full-Z quadratic spectral solver from `research/action-audit-2026-09-30`, starting at commit `afed8099d21ed07b19850b4e065ecbf32c3f9115`. The three September 30 analytical checkpoints and the `coupled-ring-axisymmetric-v0.1` implementation are unchanged.

**Result:** the localized carrier phase null mode is recovered. Longitudinal and transverse translation modes are identified and approach zero under domain/grid refinement. Their finite-box splittings remain nonzero. No physical instability or stability conclusion is accepted from this calculation.

## Implementation and scope

The separate research operator uses the full `Z ∈ [-L,L]` cylinder. Every Bogoliubov amplitude has all axial interior nodes; eigenvectors are not constrained to the stationary background's even/odd parity. The retained half-domain background is reflected, interpolated when necessary, and solved again using an independent sparse full-domain Newton method.

With `Y=(ξ,ξ̄,η,η̄)`, the audited equations are assembled as `Y_TT=B Y+G Y_T`, or `(σ²I−σG−B)Y=0`. The first-order matrix is `[[0,I],[B,G]]`. It includes both conjugate couplings, both gyroscopic frequencies, mixed time/axial advection, and the angular centrifugal operator. The m=0 axis uses the regular cylindrical stencil; nonzero m eliminates the axis with vanishing regular amplitudes. Outer perturbations obey homogeneous Dirichlet conditions.

This implementation supports N=0 and explicitly rejects nonzero background winding rather than approximating the different N±m axis spaces. The exercised parameters are c=0.6, ν=0.3, λ=25. The Φ asymptotic phase is fixed and is not incorrectly counted as a localized neutral mode.

The physical longitudinal tangent includes the carrier phase term: `(γ p_Z, γ p_Z*, γ s_Z−iνγv s, γ s_Z+iνγv s)`. Comparisons with its eigenvector remove the independent carrier phase direction. Inner products use cylindrical radial weights, including the axis stencil's discrete weight.

## Calculations completed

Seven profile/spectral runs, independently calculating m=0 and m=1, produced **208 eigenpairs**. Five runs requested 16 eigenvalues per sector; the two larger runs requested 12. Shift-invert targets were σ=0.0001 and, for confirmation, σ=0.01. These are limited searches near zero, not exhaustive searches along the imaginary axis or throughout the right half-plane.

| L | h | Longitudinal pair, σ ≈ ±value | Transverse pair, σ ≈ ±i value |
|---:|---:|---:|---:|
| 16 | 0.5 | 0.002706865 | 0.004828316 |
| 16 | 0.4 | 0.002695835 | 0.004681273 |
| 24 | 0.666667 | 0.001058992 | 0.003028487 |
| 24 | 0.5 | 0.001050733 | 0.002458371 |
| 24 | 0.4 | 0.001046668 | 0.002146329 |
| 32 | 0.5 | 0.000519913 | 0.002115692 |

The L=24, h=0.5 repeat with σ=0.01 reproduces the longitudinal value to about 3×10⁻¹² and the transverse value to about 5×10⁻¹³. It also tightens the stationary solve. The largest first-order eigenpair residual in that repeat is below 9×10⁻¹².

## Symmetry calibration

- **Carrier phase:** direct `B(0,0,is,−is)=0` residuals are at approximately 10⁻¹⁵ after tight Newton solves. The computed phase eigenvalue magnitudes are approximately 5–8×10⁻⁸, with cylindrical overlap indistinguishable from one. Earlier runs stopped Newton at a looser tolerance and give phase eigenvalue magnitudes up to approximately 10⁻⁶; these are recorded, not discarded. Direct null-vector residuals are the primary phase check, since nearly defective zero eigenvalues are more sensitive than the operator identity.
- **Longitudinal translation:** at fixed h=0.5 the real pair shrinks by a factor of about 5.2 from L=16 to L=32. At L=32 its overlap with the physical tangent, modulo carrier phase, is 0.99904. Its norm fraction in the outer 10% boundary strips is about 2.1×10⁻⁵. The shrinking positive real member is not accepted as a continuum instability.
- **Transverse translations:** the m=1 pair shrinks with increasing domain and decreasing grid spacing. At L=32, h=0.5 the tangent overlap is 0.99738. The cosine and sine angular directions give the two real translations; the calculation solves their shared m=1 radial/axial sector.
- **Grid error:** excluding the two outer grid layers, the L=24 longitudinal tangent residual drops from 0.03826 at h=0.5 to 0.02467 at h=0.4; the transverse residual drops from 0.02553 to 0.01637. These ratios are consistent with second-order finite differences. Raw whole-box translation residuals include tangent/boundary incompatibility and must not be mistaken for an exact finite-box symmetry identity.

Thus the exact carrier symmetry and translation convergence trends are recovered at this calibration stage. An exact infinite-domain translation zero mode has not been numerically established. Continued simultaneous domain/grid refinement remains necessary before accepting a physical growth mode.

## Verification and remaining gate

Three independent tests pass: regular angular-axis harmonic behavior; the Bogoliubov operator against a centered finite difference of the nonlinear equations, including gyroscopic/advection terms; and the carrier null identity on the retained profile. The largest stationary residual across the original runs is below 10⁻¹². The solver refuses to compute a spectrum for an unconverged or trivial carrier profile.

No m≥2 spectrum or broad unstable-eigenvalue scan was performed. Future instability acceptance still requires a candidate that persists under independent grid/domain refinement, has a resolved localized eigenfunction, passes eigenpair residual checks and symmetry calibration, and is reproduced without axial parity restrictions. Absence of such a candidate in this limited search does not establish spectral or nonlinear stability, continuum existence, or identification with physical matter.

## Saved artifacts

- `full_z_solver/`: standalone solver, tests, suite driver, pinned requirements, and unchanged retained reference profiles.
- `full_z_runs/`: all seven run directories, each containing the full-domain profile, spectrum metadata, and both sectors' eigenvalues/eigenvectors.
- `symmetry_calibration.json`: cylindrical overlaps, boundary fractions, tangent residuals, and eigenpair residual summaries for all 14 sector calculations.
- `verification.json`: source and output hashes, environment versions, baseline preservation checks, and test outcome.

No remote branch was updated. The new implementation and report are also retained in the isolated local checkout's September 30 research directory.
