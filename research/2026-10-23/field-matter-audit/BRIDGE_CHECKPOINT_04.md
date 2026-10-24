# Five-minute checkpoint 04: is the proposed vector coupling consistent?

23 October 2026. This tests an explicitly proposed extension. Neither a vector
field nor its coupling constant is derived from SEFI here. The existing
worldline and scalar-candidate actions remain unchanged.

## 1. Fix the convention and derive the variation

Consider a real Abelian potential A_mu and constant coupling g, with

    S_int = -g integral d^4x j^mu A_mu
          = -g integral_a^b ds J(s) T^mu A_mu(y(s)).

Under A_mu -> A_mu + partial_mu chi, holding the specified current fixed,

    delta S_int = -g integral_a^b ds J d[chi(y)]/ds
      = g J(a) chi(y(a)) - g J(b) chi(y(b))
        + g integral_a^b ds J_dot chi(y)
      = g <partial_mu j^mu, chi>.

The sign follows from the declared minus sign in S_int. For compactly
supported gauge variations with no contributing endpoints, this interaction
is invariant for arbitrary chi precisely when the source current is
distributionally conserved. Endpoint-balanced currents also meet this
condition. A finite open segment with conserved J generally has a boundary
variation; fixed boundary gauge transformations, a continued history or
appropriately transforming endpoint states require explicit accounting.

A nonzero interior torque is a different obstruction. Choosing chi to vanish
at both endpoints leaves g integral J_dot chi, which endpoint phases cannot
cancel. Removing that term requires a changed model, an additional specified
current/sector or a restricted conserved charge, not a bookkeeping change.
Other gauge-invariant interactions involving field strength are outside this
minimal source-coupling test; failure here does not exclude all vector models.

## 2. A massless Maxwell-type field imposes the same condition

If one additionally postulates

    S_A = -Z/4 integral d^4x F_munu F^munu,  Z>0 constant,
    F_munu = partial_mu A_nu - partial_nu A_mu,

then variation with the source sign above gives

    Z partial_mu F^munu = g j^nu.

Taking a second divergence gives g partial_nu j^nu=0 because distributional
derivatives commute and F is antisymmetric. Thus for nonzero g the proposed
massless field equation itself requires conserved total source. One cannot
repair an arbitrary nonconserved prescribed current by choosing a different
constant g or Z, except by trivial decoupling g=0. Massive-vector or additional
compensating-field models are different hypotheses, not ruled out by this
calculation and not silently substituted here. No Maxwell PDE was solved in
this checkpoint; this is an exact compatibility condition.

## 3. Gauge the internally symmetric sector rather than retaining its old current

Checkpoint 02 identified conditions under which z=q_O+i q_A has an internal
phase symmetry: equal paired squared frequencies and couplings independent
of that pair's phase. For that restricted sector, a possible additional
local-gauge postulate is

    z -> exp(i g chi(y)) z,
    A_mu -> A_mu+partial_mu chi,
    A_s = A_mu T^mu,
    D_s z = z_dot - i g A_s z.

Then D_s z transforms as exp(i g chi) D_s z. A candidate internal action

    S_z = integral ds [ |D_s z|^2/2 - omega^2 |z|^2/2 ]

is invariant. A phase-invariant real potential U(|z|^2) can replace the
quadratic potential without changing this symmetry statement. The already
assumed worldline geometry is kept fixed during this internal gauge test.

Expansion exposes a term that a linear-source-only prescription misses:

    |D_s z|^2/2 = |z_dot|^2/2
      - g A_s J_0 + g^2 A_s^2 |z|^2/2,
    J_0 = Im(z* z_dot).

Variation with respect to A_s yields

    dL_z/dA_s = -g J_cov,
    J_cov = Im(z* D_s z) = J_0-g A_s |z|^2.

J_cov is gauge invariant; J_0 alone is not under a local phase rotation.
Consequently the gauged worldline source is the construction of checkpoint
03 with J_cov, not automatically the old ungauged charge J_0. Do not append
another -g J_0 A_s term to the already covariantized kinetic action; its
linear interaction is already present in the expansion.

For the quadratic example the internal equation is

    D_s^2 z + omega^2 z = 0.

It implies dJ_cov/ds=Im(z* D_s^2 z)=0. A real radial potential gives the same
conclusion. This conservation law is conditional on the gauged equations and
phase-invariant potential. Existing explicit phase-breaking geometric terms
cannot simply be left unchanged and called gauge invariant. The construction
does not settle full worldline geometric dynamics, singular self-fields,
backreaction or a well-posed coupled evolution.

## 4. Actual calculations and checks

Three new tests, together with the preceding 17, gave **20 local passes in
5.41 seconds**. Saved results include these independently evaluated action
differences with g=0.7:

| Case | Direct change of S_int | Predicted change |
| --- | --- | --- |
| Constant J=2, chi vanishes at endpoints | 0 | 0 |
| J(t)=1+t/2, same endpoint-vanishing chi | 0.37333333333333324 | 28/75 |
| Constant J=2, chi=t^2 on 0<=t<=2 | -5.6 | -5.6 |

The first two use the previous timelike line and chi=t^2(2-t)^2. In the
changing-charge case the discrepancy between the directly integrated action
change and the source formula was 5.56e-17. This is a nonzero gauge variation,
not a numerical conservation success.

Thirty seeded complex-state examples additionally verify local covariance
of D_s z, invariance of J_cov, the full quadratic expansion, and the source
obtained by independently finite-differencing the kinetic action with respect
to A_s. These are algebra/variational checks, not coupled field simulations.
The result and validation JSON files record actual inputs, residuals, test
summary and code checksums. The remote CI workflow collects the new test;
remote success is not inferred from the local pass.

## 5. What this establishes and what comes next

The internally symmetric submodel admits a conventional gauge-consistent
extension under explicit additional postulates. Arbitrary phase-breaking
internal torque obstructs using that isolated current as a source for the
specified massless gauge field. This identifies a real restriction on any
attempt to consolidate the models.

It does not show that the vector field is electromagnetism, determine g or Z,
quantize charge, produce spin/statistics, derive the bulk scalar potential or
produce matter. Calling a new source coupling electric charge would require
an independent physical identification and observational comparison.

Next five-minute task: check the field energy of the simplest static point
source in this proposed Maxwell-type model. Determine whether a current
confined to an ideal worldline has finite field energy, and distinguish any
needed finite-size or regularization assumptions from a derived particle mass.
