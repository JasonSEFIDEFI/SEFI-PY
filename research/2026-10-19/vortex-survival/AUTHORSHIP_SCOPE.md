# Is Field Authorship included?

19 October 2026

**No independent Authorship coordinate is evolved in this vortex study.** The state is (Phi, Sigma, Phi_T, Sigma_T), not the four-coordinate SEFI/DEFI identity state. The vortex couplings remain fixed. A fitted phase, charge, perturbation amplitude, or displayed identity label is not silently equated with Authorship.

Three author-defined treatments must be kept distinct:

1. SEFI Formal Mathematical Construction, FIELD AUTHORSHIP: layer functionals A_i(r(t), kappa(t)), with their sum defining a scalar identity.
2. DEFI Formal Construction, coupling/action sections: an independent I_A coordinate appears in mu_1(I_A) and mu_2(I_A). Its equation is beta I_A'' = -beta I_A - mu_1'(I_A) kappa^2 I_S - mu_2'(I_A) tau_g^2 I_W. Primes on coupling functions mean differentiation with respect to I_A; time derivatives follow the source's worldline parameter.
3. SEFI-PY/core/field_authorship.py: derived direction, intensity and alignment computed from a FieldOrigin object; this class is not an independent second-order field evolution.

The DEFI equation supplies a candidate worldline internal dynamics. It does not by itself define a spatial field I_A(t,x), an observation map from Phi/Sigma to identity coordinates, a unique vortex-core worldline, coupling functions with calibrated coefficients, or the reciprocal variation of a vortex-field action. Spatial ring curvature is not automatically proper curvature of a timelike worldline. Those missing choices must be specified and justified before the two constructions can be combined.

Required bridge, before coding a full coupled model:

- Decide whether the relevant Authorship is a derived layer functional or an independent coordinate; state which source formulation is being realized.
- Define the object carrying the coordinate (worldline, string worldsheet, or spacetime field). A single center worldline does not capture the full geometry of an extended ring.
- Define its kinetic/potential terms and dimensions, admissible initial data, transformation properties, and maps to observables.
- Supply explicit coupling functions and a single variational law connecting it to the vortex variables, including reciprocal backreaction and conserved/balance laws.
- Recover the preserved vortex model in a stated decoupling limit. Vary Authorship independently with all other initial data controlled to determine whether it changes an observable.

Do not invent these choices merely to make a stable ring. The current results establish properties of the implemented coupled scalar model; they are not a test of the complete SEFI identity/geometry construction.

Source locations: the author's Journal Submissions/SEFI/SEFI_Formal_Mathematical_Construction.tex and Journal Submissions/DEFI/DEFI_Formal_Construction.tex; repository SEFI-PY/core/field_authorship.py. The foundational working manuscript of 18 October distinguishes scalar identity from the independent-coordinate extension.
