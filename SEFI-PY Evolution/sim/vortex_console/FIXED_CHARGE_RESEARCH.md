# Fixed-charge continuation of the vortex research

This advances the rest-state question without changing the published console, the engine, or the saved baseline. The aim remains a physically meaningful localized excitation; no identification with observed matter is established.

## 1. Why change the search constraint?

The previous rest-loop search held the carrier frequency fixed and could converge to zero carrier. At a nonzero fixed carrier charge, that particular endpoint has divergent constrained energy. This changes the numerical variational problem, not the underlying conservation law. It does not remove the analytical obstruction for the original coupling.

Let Phi=exp(i Omega t)F(x), Omega=sqrt(2), and Sigma=exp(i nu t)S(x). For an axisymmetric carrier take S=s(r,z) exp(i N phi). Set x=|F|^2, y=|S|^2, I=integral y d^3x, and Q=nu I.

The rotating-background functional at fixed carrier charge is

    R_Q = K + U_0 + Q^2/(2I),
    K = (1/2) integral (|grad F|^2 + |grad S|^2) d^3x,
    U_0 = integral [(x-1)^2/4 + lambda*y^2/4 + (4x-3)y/2] d^3x.

For carrier winding N, K includes integral N^2*s^2/(2r^2). The frequency is recomputed as nu=Q/I at every search step. Varying R_Q gives the same rest profile equations, with that self-consistent frequency.

IMPORTANT: R_Q is E_relative - Omega*Delta Q_Phi, not the background-subtracted physical energy alone. The excess charge Delta Q_Phi=Omega*integral(|F|^2-1) is not held fixed in this search. Thus even a minimum of R_Q is not yet a demonstrated isolated-particle energy minimum. Full constrained stability must address both fields' conserved quantities and the rotating background.

## 2. Dilation test and its limitations

For F_l(x)=F(x/l), S_l(x)=S(x/l), with carrier charge fixed:

    R_Q(l) = l K + l^3 U_0 + l^(-3) C,  C=Q^2/(2I).

A critical point must satisfy

    K + 3 U_0 - 3 C = 0.

Its second derivative in this restricted dilation direction is 6 U_0 + 12 C. For lambda>9, U_0 is nonnegative for the selected rotating vacuum, so this direction has positive curvature at a nonzero-charge critical point if one exists. This is not an existence proof and does not test other modes. Substituting C=nu^2 I/2 recovers K+3 integral U_nu=0, so it is fully consistent with the earlier no-go result.

## 3. Exact pointwise binding threshold

For y>0, minimize

    2 U_0(x,y)/y = (x-1)^2/(2y) + lambda*y/2 + 4x - 3

over x>=0 and y>0. Minimizing in y gives sqrt(lambda)*|x-1|+4x-3. Therefore, for lambda>9,

    mu_*^2 = min(1, sqrt(lambda)-3).

For 9<lambda<16 the minimum is attained at x=0, y=1/sqrt(lambda); for lambda>16 the value 1 is an infimum approached at the undepleted background with vanishing carrier. At lambda=16 there is a degenerate minimizing family.

Consequently U_0 >= mu_*^2 I/2 and

    R_Q >= K + (mu_*^2 I + Q^2/I)/2 >= mu_* |Q|.

The dilute carrier threshold in this rotating-background functional is |Q| because its far-field mass is 1. At lambda=25, mu_*=1: this bound supplies no below-threshold binding. At lambda=12, mu_*=sqrt(sqrt(12)-3)=0.6812500386: below-threshold binding is energetically permitted by this bound, but not demonstrated by it. Gradient and interface costs, angular momentum, both charge constraints, and actual solution existence remain decisive.

For the stated ansatz, the carrier's axial angular momentum has magnitude |J_z|=|N Q|. With Q=Im(Sigma* d_t Sigma), exp(+i nu t+iN phi), and the usual active spatial-rotation generator, the signed relation is J_z=-NQ. This classical relation does not derive quantized particle spin; Q remains continuous, and F has no azimuthal dependence in this ansatz.

## 4. New numerical method

`work/fixed_charge_rest_search.py` implements an axisymmetric finite-volume-compatible energy and its negative gradient, with Q fixed and nu updated as Q/I. An energy-decreasing line search evolves a fictitious relaxation parameter. This is NOT physical time evolution and must not be shown as such in the console.

The initial original field contains a meridional vortex ring; the carrier has N=1 azimuthal winding. Outer boundaries are F=1, s=0 on a finite cylinder. The original field is regular on the axis, and s=0 there for N>0. No ring-radius or original-field topology constraint is imposed. This allows loss of the vortex to be detected rather than artificially prevented.

Directional finite-difference checks verified the implemented energy gradient to relative discrepancies 1.29e-10 (lambda=12) and 4.71e-11 (lambda=25) at the initial configuration. This checks the discrete variational implementation, not the physical model.

At L=20, h=.5, Q=100, N=1, after 12,000 relaxation steps:

| Coupling | R_Q/|Q| | nu | Maximum residual | Original midplane ring cores | Carrier RMS radius |
|---|---:|---:|---:|---|---:|
| lambda=12 | 1.0165173 | 1.0109154 | 4.23e-5 | none | 12.6419 |
| lambda=25 | 1.0299443 | 1.0377394 | 4.03e-5 | none | 14.6047 |

Neither met the residual tolerance. Both frequencies violate the |nu|<1 requirement for exponential carrier localization in the infinite-domain background. These are unsuccessful rest-vortex searches tending toward diffuse finite-box configurations, not new localized solutions. They do not prove nonexistence at other charges or from other initial data. Energy decreases monotonically by construction; that is not evidence of physical dynamical stability.

## 5. Larger-charge result and checks

At lambda=12, Q=1000 and N=1, the same initial ring relaxed toward a localized toroidal carrier distribution. Continuing that state, enlarging the domain, and refining the grid produced:

| L | h | R_Q | nu | Carrier RMS radius | Maximum residual |
|---:|---:|---:|---:|---:|---:|
| 20 | .5 | 907.870456 | .83960666 | 8.818318 | 9.9973e-6 |
| 24 | .5 | 907.863919 | .83958557 | 8.820394 | 9.9972e-6 |
| 24 | .4 | 907.923579 | .83962216 | 8.821415 | 9.9932e-6 |

All three met the chosen 1e-5 residual tolerance. These are preliminary discrete critical-point candidates, not a continuum existence theorem. The continuum virial residual was -.4873, -.3432 and -.1996 respectively: it is not exactly zero and improves with the domain/grid checks. At L=24,h=.5, the carrier maximum lies at r=5,z=0; only 1.90e-5 of its charge integral lies beyond 80% of the radial or axial boundary extent.

Crucially, the original field has minimum amplitude .20306 in the wider-domain run and no zero anywhere on that grid. Its original vortex has disappeared. Carrier winding N=1 remains part of the ansatz. Thus this is a candidate spinning charged torus / Q-ring-like branch, not a demonstrated vorton retaining the original vortex. A resemblance to known Q-ring branches is a classification suggestion, not proof of model equivalence.

The constrained energy per carrier charge is about .9079, below the free dilute-carrier value 1 in the same rotating-background ensemble. This supports studying a possible binding branch. It does not prove stability against splitting, radiation into other sectors, or general non-axisymmetric perturbations. It does not imply a threshold between Q=100 and Q=1000 has been located; different initial data could find additional branches.

## 6. Background-charge audit: why this is not a mass measurement

For the wider-domain candidate, direct integration gives

    integral (|F|^2-1) = -4480.710656,
    Delta Q_Phi = -6336.681779,
    R_Q = 907.863919,
    E_relative = R_Q + Omega*Delta Q_Phi = -8053.557393.

The negative last value is energy relative to an energetic homogeneous rotating background with a different original-field charge. It is not a negative absolute energy or a negative particle mass. The large background-charge deficit means that quoting R_Q as an isolated rest energy would silently change the physical comparison.

There is a useful next construction. On a finite domain, let I_F=integral |F|^2, and compensate the excess original-field charge in the initial velocity by a uniform extra phase rate deltaOmega=-Delta Q_Phi/I_F. The added energy is -Omega*Delta Q_Phi+(Delta Q_Phi)^2/(2 I_F), so the compensated energy relative to the original background becomes

    E_compensated = R_Q + (Delta Q_Phi)^2/(2 I_F).

This is an initial-data identity, not a stationary solution: changing the phase rate requires a new self-consistent field solution or physical evolution. A localized compensation window is preferable if the outer phase velocity must remain fixed; its quadratic cost depends on the window's charge-storage susceptibility. Increasing a compensation volume can lower that cost, demonstrating the importance of the surrounding medium. A particle interpretation requires separating intrinsic excitation properties from this background dependence.

## 7. Next discriminator

The immediate next test is a two-charge-consistent treatment of this candidate and its non-axisymmetric stability, followed by comparison with a genuinely vortex-preserving branch if one can be found. The new result justifies that test; it does not justify inserting this relaxation as physical evolution into the existing console. Keep the published moving-ring baseline intact, and label this separate branch explicitly if a future console version exposes it.

Reproduction: run `work/fixed_charge_rest_search.py --lam 12 --Q 1000 --name fixed_charge_lambda12_Q1000`, then use `--initial work/fixed_charge_lambda12_Q1000.npz` with `--L 24` and subsequently `--h .4` as recorded in the JSON results. The script uses NumPy only. The 55-second runtime cap can terminate early on slower hardware; inspect the recorded reason and residual rather than assuming convergence. Saved relaxation trajectories are not physical-time trajectories.

Numerical evidence: `fixed_charge_lambda12_Q100.json`, `fixed_charge_lambda25_Q100.json`, `fixed_charge_lambda12_Q1000.json`, `fixed_charge_Q1000_continued.json`, `fixed_charge_Q1000_wide.json`, `fixed_charge_Q1000_refined.json`, and `fixed_charge_background_audit.json` in outputs. Full-precision fields are in corresponding work NPZ files. The analytical bound was also checked on 100,000 sampled (x,y) pairs at each of four couplings; this arithmetic check supplements, not replaces, the derivation.

Related primary research studies vortons and their decay into Q-rings or spinning Q-balls in other models: Tallarita, Peterson, Bolognesi and Bedford, *Vortons with Abelian and non-Abelian currents and their stability* (2019), https://arxiv.org/abs/1909.01950. Those results motivate distinguishing branches; they do not establish existence or stability in this rotating-background model.
