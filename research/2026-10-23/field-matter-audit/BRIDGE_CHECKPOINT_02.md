# Five-minute checkpoint 02: when do internal phase symmetries exist?

23 October 2026. This is a conditional derivation from the existing diagonal
internal kinetic metric and quadratic stiffnesses. It does not construct
spacetime fields, quantize the model, or settle higher-derivative worldline
stability. No existing action is modified.

## 1. Normalization matters

The internal part of the foundations action is

    L_int = sum_a h_a I_dot_a^2/2 - sum_a k_a I_a^2/2
            - mu1(I_A) kappa^2 I_S - mu2(I_A) tau^2 I_W,

with positive constant h_a,k_a and dots denoting proper-time derivatives.
Consider independent pair rotations (O,A) and (S,W), keeping the worldline
fixed. An ordinary Euclidean rotation of the original I coordinates preserves
the uncoupled kinetic and potential forms precisely when

    h_O=h_A and k_O=k_A;
    h_S=h_W and k_S=k_W.

These equalities are sufficient but unnecessarily restrictive for a physical
phase interpretation with field normalization. Define instead

    q_a = sqrt(h_a) I_a,     omega_a^2 = k_a/h_a.

The kinetic term becomes sum q_dot_a^2/2. Independent ordinary rotations of
(q_O,q_A) and (q_S,q_W) preserve the uncoupled action precisely when

    k_O/h_O = k_A/h_A,
    k_S/h_S = k_W/h_W.

Necessity follows because a diagonal quadratic form invariant under every
rotation in a plane must be proportional to the identity in that plane.
Sufficiency follows directly by writing it as omega^2 times the squared
radius. In the original coordinates these are weighted, generally elliptic
rotations. Unequal h values therefore do not alone obstruct phase symmetry.

The explicit example h=(1,4,9,16), k=(2,8,27,48) has squared frequencies
(2,2,3,3). It fails both ordinary-coordinate equality tests but passes both
normalized symmetry tests. The equal-frequency conditions themselves remain
extra parameter restrictions unless SEFI's premises imply them.

## 2. Derive the charges and their source terms

In normalized coordinates define the candidate internal angular charges

    J_OA = q_O q_dot_A - q_A q_dot_O,
    J_SW = q_S q_dot_W - q_W q_dot_S.

Let mu1' and mu2' denote derivatives with respect to I_A=q_A/sqrt(h_A).
Direct substitution of the actual internal Euler-Lagrange equations gives

    J_dot_OA = (omega_O^2-omega_A^2) q_O q_A
      - q_O/sqrt(h_A) [mu1' kappa^2 I_S + mu2' tau^2 I_W],

    J_dot_SW = (omega_S^2-omega_W^2) q_S q_W
      + q_W mu1 kappa^2/sqrt(h_S)
      - q_S mu2 tau^2/sqrt(h_W).

These are exact torque identities on the stated internal equations, including
along dynamical worldlines where those equations hold. The numerical examples
prescribe geometry only to isolate and test this part of the system.

Consequences:

* Equal paired frequencies and zero geometric couplings give two exact
  conserved internal charges.
* Constant mu1,mu2 remove the Authorship-derivative torque. With equal OA
  frequencies, J_OA is conserved even when geometry varies.
* Nonzero constant couplings generally torque the origin-centered SW pair.
  For a symmetry on all internal states at a fixed geometry, both forcing
  coefficients mu1 kappa^2 and mu2 tau^2 must vanish. An accidental zero
  torque on one trajectory is not an action symmetry.
* Authorship-dependent couplings generally source both candidate charges.
  Special histories or restricted sectors may cancel these terms; the
  general formulas must be used rather than assuming conservation.

For h,k above, q=(1,2,3,4), kappa=2, tau=3:

| mu1(I_A), mu2(I_A) | J_dot_OA | J_dot_SW |
| --- | --- | --- |
| 0, 0 | 0 | 0 |
| 1, 2 | 0 | -8.166666666666668 |
| 1+I_A/2, 2-I_A/4 | 0.125 | -3.8125 |

These are illustrative exact-equation evaluations, not observations or
particle properties. `internal-symmetry-results.json` records the inputs.

## 3. A surviving shifted symmetry limits the negative conclusion

Take constant mu values and prescribed constant geometry, with equal SW
frequencies omega_S^2=omega_W^2=omega_b^2. In normalized coordinates the SW
potential is

    V_SW = omega_b^2 (q_S^2+q_W^2)/2 + B_S q_S + B_W q_W,
    B_S=mu1 kappa^2/sqrt(h_S), B_W=mu2 tau^2/sqrt(h_W).

Its center is q_star=-B/omega_b^2. Completing the square leaves an isotropic
oscillator about that center, so

    J_shifted = (q_S-q_star_S) q_dot_W
                - (q_W-q_star_W) q_dot_S

is conserved. The same example that gives nonzero J_dot_SW passes this shifted
charge check. Thus the prior failed origin-centered pairing does not exclude
every affine identification or every internal phase symmetry.

If geometry varies, q_star generally depends on proper time. Translating to
that moving center introduces derivative/forcing terms. The fixed-background
argument does not establish a symmetry of the full geometric system. Joint
internal/geometric transformations and emergent symmetries are not classified
by this checkpoint.

## 4. What this does and does not connect

There is now an explicit internally symmetric submodel and an exact account
of how the declared couplings break its charges. This gives a sharper bridge
hypothesis than pairing four names into two complex variables. However:

* An internal oscillator on a curve is not yet a spacetime field with spatial
  gradients, a local bulk current, or a quantum particle sector.
* The positive quadratic oscillator potential, including a shifted version,
  does not by itself derive the bulk quartic portal or its symmetry-breaking
  vacuum. Those interactions still require a derivation or declared postulate.
* A change of field normalization cannot create the binding well. Under
  Phi_new=r Phi, Sigma_new=t Sigma with positive r,t, the quartic parameters
  become a_new=a/r^4, b_new=b/t^4, g_new=g/(r^2 t^2), so
  a_new b_new-g_new^2=(ab-g^2)/(r^4 t^4). Its sign is unchanged.
* These internal angular charges have not been identified with electric
  charge, spacetime spin or the bulk candidate's charge.

## 5. Actual checks and next bounded task

Four new tests check normalized versus unnormalized conditions, finite
rotations, independent numerical differentiation of the potential and charge
torques over 20 seeded inputs, and the surviving shifted symmetry. Together
with the prior nine tests, **13 local tests passed in 5.39 seconds**. The exact
ratio conditions use rational arithmetic; gradient/torque checks use NumPy.
No failed numerical outcome was discarded. The remote workflow is updated
to collect all checkpoint tests; its result must be checked separately.

Next five-minute task: construct the spacetime distributional current carried
by a worldline with one of these internal charges. Derive its divergence with
endpoints and nonzero torque included, and test exactly when local charge
conservation follows. This would establish a limited mathematical link to
field-theory language without assuming that it supplies a smooth bulk action.
