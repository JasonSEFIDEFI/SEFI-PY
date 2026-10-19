# Vortex benchmark selection and first exact discriminator

19 October 2026

## Choice

Use the **ungauged/global Witten two-complex-scalar model** with a static symmetry-breaking vacuum as the primary comparator. The literature anchor is Tallarita, Peterson, Bolognesi and Bedford, *Vortons with Abelian and non-Abelian currents and their stability*, EPJC 79, 861 (2019), https://arxiv.org/html/1909.01950. Its global limit is taken by setting the gauge coupling to zero and discarding the decoupled gauge field. The paper also studies gauged and non-Abelian variants; those are not silently equated with this benchmark.

This is a close structural comparator, not a claim of the largest difference among all possible vortex theories. There is no defined universal ranking of model difference. A simple single-field vortex would omit the carrier by construction and make a less informative primary comparison. The difference selected here is rotating versus static asymptotic state, measured through linear response.

## Matched normalization and potentials

Use dimensionless coordinates, metric wave speed 1, and canonical real-component kinetic terms: L = (|d_t Phi|^2 - |grad Phi|^2 + |d_t Sigma|^2 - |grad Sigma|^2)/2 - U. Define rho=|Phi|^2 and eta=|Sigma|^2.

The control potential is

    U0 = (rho-1)^2/4 + lambda*eta^2/4 + (4*rho-3)*eta/2.

It is of the global Witten quartic two-field form. Coefficients have been chosen to match the console's effective rest functional; this is not a reproduction of the published paper's numerical parameter runs. With its full complex kinetic convention the paper's potential is 2 U0: lambda_phi=2, lambda_sigma=2 lambda, eta_phi=1, eta_sigma^2=3/lambda, gamma=4. Multiplying the whole Lagrangian by a common constant does not change its classical equations, although charge normalization must be tracked.

The console rotating realization instead corresponds to

    Urot = U0 + Omega^2*rho/2,   Omega=sqrt(2),
    Phi = exp(i Omega t) F,    F -> 1, Sigma -> 0 at infinity.

This reconstruction follows from the implemented rest-frame evolution at c=0, not from a new fundamental SEFI derivation. In that frame, the console has

    F_tt = Laplacian F - 2i Omega F_t + (1-|F|^2-4|Sigma|^2)F.

For the static control, Omega=0 and the same rest effective potential U0 is used. Both have the same effective Phi amplitude stiffness at the selected background and the same dilute carrier mass squared, 1. They have different physical potentials and background Phi charge densities (0 versus Omega in this normalization). The experiment compares these specified physical setups; it is not an equal-charge test of the same isolated system.

## Derivation of the discriminator

Linearize F=1+a+i b with a,b real and Sigma=0. The carrier decouples to first order. For a Fourier mode of wave number k:

    a_tt + (k^2+2)a - 2 Omega b_t = 0,
    b_tt + k^2 b + 2 Omega a_t = 0.

Substitution of exp(-i omega t) gives

    (k^2+2-omega^2)(k^2-omega^2) - 4 Omega^2 omega^2 = 0,
    omega_+/-^2 = k^2 + 1 + 2 Omega^2
                   +/- sqrt[(1+2 Omega^2)^2 + 4 Omega^2 k^2].

Thus the static control has omega_- = k and omega_+^2=k^2+2. The rotating realization has

    omega_+/-^2 = k^2 + 5 +/- sqrt(25+8k^2),
    low-k phase and group speed = 1/sqrt(5) = 0.4472135955.

The static control's low-k speed is 1. The rotating setup therefore has a **55.2786% smaller long-wavelength gapless-mode speed** under this matching convention. Its gapped co-rotating mode starts at sqrt(10), versus sqrt(2). These are small-disturbance eigenfrequencies relative to each selected background. They are not altered causal light cones: the principal wave operator and high-frequency signal-front speed remain 1 in both models. Mixing of real amplitude and phase perturbations in a rotating condensate explains the changed spectrum.

This result describes waves in the homogeneous region outside a vortex. It is a possible input to vortex radiation and stability studies, not a calculated ring lifetime or a direct head-to-head vorton simulation. It does not establish a uniquely SEFI effect: ordinary charged scalar-field theory can describe rotating backgrounds too.

## Verification and artifacts

Run `compare.py` with NumPy and Pillow. It saves `results.json` and `dispersion-comparison.png`. Twelve mode cases compare the exact dispersion against independent eigenvalues of the first-order linear system and against RK4 evolution at two time steps. These are internal mathematical/numerical checks, not experimental validation or spatial-grid convergence tests. There is no spatial discretization in this Fourier-mode test.

## What the existing vortex evidence says

The preserved moving-ring console evolves lambda=25, carrier azimuthal winding N=0; its nonzero traveling parameter is distinct from a stationary spinning vorton. The separate lambda=12, Q=1000, N=1 fixed-charge relaxation candidate loses the original Phi vortex: its minimum amplitude is nonzero on the inspected grid. It must be labeled a charged torus / Q-ring-like candidate, not a vortex-preserving stable particle. Relaxation time is not physical time.

The literature comparator supplies explicit equations and numerical constructions; it does not merely suggest vortons. Distinct equations or background conditions do not themselves show superior predictive accuracy.

## Next fair nonlinear benchmark

1. Identify whether the desired comparison is the preserved traveling-ring branch or the separate spinning charged-torus branch. Do not mix their results.
2. Use the same physical normalization, domain, grid, carrier charge and carrier winding; account explicitly for original-field background charge and energy. Re-solve each model rather than copying an equilibrium from one into the other.
3. Report convergence, field zeros, meridional winding on nonsingular contours, carrier distribution and energy/charge/momentum bookkeeping.
4. Apply equal normalized disturbances and physical-time evolution. A non-axisymmetric three-dimensional test is necessary before asserting full stability. Existing axisymmetric evolution is not validated for the lambda=12, N=1 branch.
5. Compare measured ring survival, radiation spectra and mode growth rates only after those runs exist. Preserve failed and no-vortex outcomes.

The result established here is the exact far-field dispersion difference, not completion of these nonlinear tests or proof of empirical novelty.

## Local evidence inspected

- `SEFI-PY Evolution/sim/vortex_console/physics.py`: implemented dynamics and allowed evolution parameters.
- Same directory, `FIXED_CHARGE_RESEARCH.md`: effective rest functional, numerical branch results, binding bounds and caveats.
- Same directory, `TWO_CHARGE_ACCOUNTING.md`: background charge compensation and its energy cost.
- Same directory, `MODEL_HISTORY.md`: separation of physical evolution and fictitious relaxation.

These are preserved source snapshots, not modified by this comparison.
