# Finite-time moving-ring survival and charge-balance audit

Jason Duran Dutton | 19 October 2026

## Research question and scope

Does the preserved lambda=25, N=0 traveling coupled scalar ring retain a numerically detected core under selected disturbances? Do its discrete energy and charge balances close? How do selected linear angular disturbances behave?

This study does **not** evolve an independent Field Authorship coordinate. See AUTHORSHIP_SCOPE.md. It tests the existing two-complex-field realization, not the complete SEFI/DEFI theory or physical matter.

## Results

- Eight axisymmetric physical-time runs reached dimensionless T=24, retaining at least one detected original-field vortex-core cell at every saved sample. Inputs include unperturbed control; Gaussian amplitude disturbances 0.001, 0.01, 0.02; a 0.02 phase disturbance; half-step refinement; and saved coarse/wide reference profiles.
- At the final time the medium disturbance's full state differs by 7.37e-8 when the time step changes from 0.025 to 0.0125. This is a two-step consistency comparison, not a measured convergence order.
- Maximum carrier-charge relative drift across the eight runs is below 1.88e-9. Carrier density stays mostly localized by the recorded RMS and outer-region diagnostics.
- Four selected linear cos(m theta) perturbations, m=1,2,3,4, have maximum combined field/velocity norm gains about 1.71, 1.44, 1.43, 1.43. These are finite-time responses to one chosen initial direction per m, not eigenvalue bounds or proof that all modes are stable. The field/velocity norm is dimensionless and is not a physical energy norm. Unit-normalized tangent data are interpreted only through linearity, not as a unit-sized nonlinear disturbance.
- For m=2, halving the time step changes the final full tangent state by 2.77e-7. Its final gain is 1.2540 on the refined reference, 1.2586 on the coarse reference and 1.2575 on the wider reference. These are sensitivity checks, not a complete continuum/domain limit.
- For medium amplitude disturbances on the refined and wide references, the discrete rotating/moving-frame energy drifts by less than 2.48e-10 relatively. Original-field excess charge changes by -2.76879 and +0.0185645 respectively, while its independently integrated finite-box balance closes to absolute residual below 1.87e-9.

## What these findings mean

The finite-time tests did not destroy the detected core. Time refinement, selected angular perturbations and balance closure supply stronger evidence than a static picture or residual alone. They do not establish a stable isolated particle.

The fixed outer Phi value is not globally phase-invariant. The corresponding finite-box original-field charge need not be conserved, even though the discrete rotating-frame energy is. Its measured change is accounted for by the derived boundary balance. The energy used here is the rotating/moving-frame functional, not a Lorentz-invariant rest mass.

The coarse and refined references have the same L=16 and different h (2/3 and 1/2). The coarse and wide references share h=2/3 but differ in L (16 and 24). Each is independently solved, not the same interpolated initial state. Their saved starting radii are approximately 5.366, 5.357 and 5.157: domain dependence remains material. Therefore these checks do not establish domain-independent intrinsic ring size.

The rough boundary-influence warning times are 8.6 for L=16 and 15.3 for L=24. Runs extend past those times and have reflecting/fixed boundaries; they cannot be presented as radiation into infinite space. Initial profiles have noncompact tails, so these times are not rigorous isolation windows.

Core detection uses the summed complex phase increments around grid plaquettes, requiring nonzero corner amplitudes. It reports containing cells, not exact continuum zeros. Grid-corner zeros are flagged as ambiguous; finite spatial sampling can miss small structures. Detection at saved times is not a theorem of uninterrupted topology between samples. The older fixed contour winding is also retained in the raw diagnostics, but it is not the sole survival criterion.

## Equations and checks

The production physics.py is unchanged. The angular tangent equations are the real-linear directional derivative of its nonlinear force plus the cylindrical angular Laplacian -m^2/r^2. Since the baseline has carrier azimuthal winding N=0, a cos(m theta) sector closes at linear order; this must not be reused unchanged for N=1. Regular m>0 perturbations vanish on the axis.

For d=Phi-1 and eta=|Sigma|^2, the discrete conserved energy is

    H = integral [ (|Phi_T|^2+|Sigma_T|^2)/2
          - Re(d* Laplacian d + Sigma* Laplacian Sigma)/2
          + (|Phi|^2-1)^2/4 + 25 eta^2/4
          + (4|Phi|^2-3-nu^2)eta/2
          + c Im(d* d_Z)/2 ].

Stars here denote complex conjugation. The integral uses the console's cylindrical weights. The gradient term is a discrete integration-by-parts form, not an extra physical interaction. Subtracting the fixed background before the current term avoids a spurious boundary linear term.

Let v=c/sqrt(8+c^2), gamma=(1-v^2)^(-1/2), Omega=sqrt(2). Excess Phi charge is

    Delta Q = integral [(Omega(|Phi|^2-1)+Im(Phi* Phi_T))/gamma
                        - v Im(Phi* Phi_Z)].

Its derived discrete balance rate is

    Delta Q_T = integral Im(Phi* Laplacian Phi)/gamma
              + 2 Omega v integral Re(Phi* Phi_Z)
              + v integral Im(Phi* (Phi_T)_Z - Phi_T* Phi_Z).

Constitutive forces cancel in this rate. These finite-box terms reduce to boundary contributions under the relevant discrete summation identities. They are not a new source of energy or a new charge violation in the continuum theory. Directional differentiation checks this rate and verifies H_T=0 at a perturbed state with nonzero velocities. RK4 quadrature independently integrates the rate in balance.py.

Four new unit tests pass: core/negative control, corner ambiguity, tangent versus nonlinear derivative, energy/charge differential identities. Existing production console regressions are run separately before publication.

## Reproduction

Requires Python and NumPy. From this directory:

    python -m unittest test_study -v
    python study.py
    python azimuthal.py
    python -c "from azimuthal import run_mode; run_mode(2,reference='trial'); run_mode(2,reference='wide')"
    python balance.py

summary.json and angular_summary.json summarize the runs. Case JSON files retain sampled diagnostics; NPZ files retain initial/final nonlinear fields or final tangent states. They do not contain the complete time-resolved 3D field. Re-running overwrites these study outputs, never the preserved console reference data.

## Next decisive work

1. Specify the Authorship-to-vortex variational bridge from the author's sources before claiming a full SEFI test.
2. Expand/refine independently solved domains and close charge/energy/momentum budgets with open-boundary flux tracking.
3. Resolve the perturbation spectrum and nonlinear 3D outcomes, including symmetry-breaking failure modes.
4. Treat both conserved charges consistently for the separate N=1 charged-torus branch before extending its physical evolution.
5. Only then compare isolated energy/momentum response, scattering and observable signatures against conventional models. Matter identification additionally requires dimensional calibration and experimental evidence.

No particle mass, charge quantization, quantum statistics, universal coupling, or experimental matter identification is established here.
