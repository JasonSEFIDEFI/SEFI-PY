# Action-to-traveling-wave audit

**Date:** 2026-09-30  
**Status:** analytical checkpoint; no new numerical or experimental claim  
**Baseline implementation preserved:** `coupled-ring-axisymmetric-v0.1`

## Purpose

Re-derive the current coupled-ring equations and traveling variational structure from the four-dimensional action before using continuation data. This checkpoint deliberately separates model assumptions from derived identities.

## Model definition

Metric signature `(+---)` and

[
S=int d^4x,mathcal L,
qquad
mathcal L=
rac12partial_muPhi^*partial^muPhi+
rac12partial_muSigma^*partial^muSigma-V,
]

[
V=
rac12|Phi|^2+rac14|Phi|^4+
rac{lambda}{4}|Sigma|^4+
rac12(4|Phi|^2-3)|Sigma|^2.
]

These are model assumptions/definitions.

## Derived field equations

[
BoxPhi+(1+|Phi|^2+4|Sigma|^2)Phi=0,
]

[
BoxSigma+(lambda|Sigma|^2+4|Phi|^2-3)Sigma=0.
]

The unit-amplitude rotating background (Phi=e^{iOmega t}), (Sigma=0) requires

[
Omega^2=2.
]

## Traveling ansatz

[
Z=gamma(z-vt),qquad gamma=(1-v^2)^{-1/2},
]

[
Phi=e^{iOmega t}psi(r,Z),
qquad
Sigma=e^{i
ugamma(t-vz)}e^{iNarphi}s(r,Z).
]

Direct substitution gives

[

abla^2psi+ic,partial_Zpsi+
(1-|psi|^2-4s^2)psi=0,
]

[

abla^2s-rac{N^2}{r^2}s+
(3+
u^2-4|psi|^2-lambda s^2)s=0,
]

with

[
c=2Omegagamma v,qquad
v=rac{c}{sqrt{8+c^2}}
]

for (Omega=sqrt2). These match the current solver.

## Conserved currents and stress tensor

Independent global phase symmetries give

[
j^mu_Phi=operatorname{Im}(Phi^*partial^muPhi),
qquad
j^mu_Sigma=operatorname{Im}(Sigma^*partial^muSigma).
]

Translation invariance gives

[
T^{mu
u}
=
operatorname{Re}(
partial^muPhi^*partial^
uPhi+
partial^muSigma^*partial^
uSigma)
-eta^{mu
u}mathcal L.
]

At fixed laboratory time, (dz=dZ/gamma). The homogeneous background has

[
T^{00}_{m bg}=rac74.
]

Define the finite relative quantities

[
E_{m exc}=int d^3x,(T^{00}-7/4),
]

[
Delta Q_Phi=int d^3x,(j^0_Phi-Omega),
qquad
Q_Sigma=int d^3x,j^0_Sigma,
]

and physical contravariant momentum

[
P^z=int d^3x,T^{0z}.
]

## Exact constrained traveling functional

The helical/traveling generator is

[
k=partial_t+vpartial_z.
]

For the ansatz,

[
kPhi=iOmegaPhi,
qquad
kSigma=irac{
u}{gamma}Sigma.
]

Therefore the constrained physical combination is

[
H_{m tr}
=
E_{m exc}-vP^z-OmegaDelta Q_Phi-rac{
u}{gamma}Q_Sigma.
]

After changing from (z) to (Z), the reduced functional is

[
mathcal G_c=gamma H_{m tr}
=K_perp+K_Z+U+rac c2P_{m red},
]

where

[
P_{m red}
=int d^3X,operatorname{Im}(psi^*partial_Zpsi),
]

and

[
U=int d^3X,
left[
rac14(|psi|^2-1+4s^2)^2+
rac{lambda-16}{4}s^4+
rac{1-
u^2}{2}s^2
ight].
]

Thus the factor multiplying reduced momentum is (c/2), not (c). An earlier exploratory derivation with (mathcal G=K+U+cP_{m red}) was incorrect by a factor of two and is superseded by this checkpoint.

The reduced momentum (P_{m red}) is not, in general, the physical Minkowski momentum (P^z).

## Pohozaev identities

For localized continuum solutions with boundary terms vanishing, independent transverse and longitudinal scalings give

[
K_Z+U+rac c2P_{m red}=0,
]

[
K_perp-K_Z+U=0.
]

Their isotropic consequence is

[
K+3U+cP_{m red}=0.
]

At (lambda=25) and (|
u|le1), the constrained potential density is nonnegative. At (c=0), the identities force (K_Z=U=K_perp=0), reproducing the rest-state obstruction for this ansatz.

## Branch first law

At each exact critical solution, the first variation of (H_{m tr}) with respect to the physical Cauchy data vanishes. Therefore a tangent variation along any smooth family of such solutions satisfies

[
dE_{m exc}
=
v,dP^z+Omega,d(Delta Q_Phi)
+rac{
u}{gamma},dQ_Sigma.
]

This is the appropriate branch identity when (Omega) and (
u) are the fixed ansatz parameters and (v) varies along the family. Consequently, (dE_{m exc}/dP^z=v) is justified only on a subfamily for which the charge variations vanish (or after the corresponding charge terms are explicitly removed). This prevents a naive Gross-Pitaevskii-style derivative test from being applied to the wrong conserved quantities.

## Falsification tests fixed before continuation data

A candidate continuum branch must simultaneously show:

1. discrete field-equation residual convergence;
2. stable observables under independent grid and domain refinement;
3. convergence of both independent Pohozaev residuals toward zero as boundary/discretization errors decrease;
4. consistency with the branch first law above;
5. explicit separation of (P_{m red}) from (P^z).

Failure of these tests is to be recorded as a negative result, not repaired by redefining observables after inspection.

## Claims explicitly not established

This checkpoint does not establish infinite-domain existence, spectral/nonlinear stability, particle mass, spin, electric charge, fermionic statistics, a universal effective metric, gravity, or identification with observed matter.
