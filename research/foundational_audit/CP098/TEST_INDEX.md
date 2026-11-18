# CP098 — Reproducibility and Test Index

Date: 2026-11-17
Status: OPEN
Parent checkpoint: ../CHECKPOINT_098_STABILITY.md
Branch: research/emergent-consciousness-math

## Evidence policy

- REPORTED PASS: terminal output reviewed in research conversation.
- ARCHIVED: executable source and expected result stored in Git.
- REPRODUCED: archived source rerun successfully.
- PENDING: calculation unfinished or result unavailable.
- CORRECTED: earlier result superseded with explanation.

A reported pass is not an archived reproducible test.
No result is promoted to a general proof merely because
a symbolic test passes on a restricted family.

## Calculation register

| ID | Calculation | Evidence | Archive |
|---|---|---|---|
| CP098-T01 | Normal-secondary bracket reconstruction | REPORTED PASS | PENDING |
| CP098-T02 | Lorentz-boost covariance test | REPORTED PASS | PENDING |
| CP098-T03 | Primary multiplier coefficients | REPORTED PASS | PENDING |
| CP098-T04 | A=0 candidate tertiary independence | REPORTED PASS | PENDING |
| CP098-T05 | A=0 and D_A evolution tests | REPORTED PASS | PENDING |
| CP098-T06 | A=0 constraint compatibility counterexample | REPORTED PASS | PENDING |
| CP098-T07 | Tangent tertiary multiplier brackets | REPORTED PASS | PENDING |
| CP098-T08 | Tangent tertiary isolated-point tests | REPORTED PASS | PENDING |
| CP098-T09 | One-parameter constrained family | REPORTED PASS | PENDING |
| CP098-T10 | Two-parameter constrained family | REPORTED PASS | PENDING |
| CP098-T11 | Changed-geometry constrained family | REPORTED PASS | PENDING |
| CP098-T12 | Exact Hamiltonian decomposition | REPRODUCED | LOCAL FILE |
| CP098-T13 | Constraint-algebra reconstruction | REPORTED PASS | PENDING |
| CP098-T14 | General tangent-secondary bracket | REPRODUCED | LOCAL FILE |
| CP098-T15 | Normal primary-secondary bracket sign | REPRODUCED | LOCAL FILE |


## Established algebraic identity

In the regular domain v>0, kappa>0, C!=0:

H_c = v psi_T + v^2 kappa chi_N - (T.q) chi_T

This identity was checked symbolically without
imposing the constraints.

Executable source: tests/CP098_T12.py
Locally reproduced with Python 3.12.3 and SymPy 1.14.0.
Git archival remains pending.

## Correction register

### CP098-C01: Normal multiplier sign

Original checkpoint reports:

{chi_N,phi_N} ~= +A/32 at v=4

Subsequent calculations report:

{chi_N,phi_N} ~= -A/32 at v=4

The later result supersedes the earlier sign claim,
subject to independent reproduction from archived code.

### CP098-C02: Covector contractions

Canonical covector/vector contractions use p_mu X^mu.

Metric contractions between covectors use eta^mu nu.

A previous two-parameter script incorrectly applied
the Minkowski metric to a covector/vector contraction.
The script was corrected and subsequently passed.

### CP098-C03: Momentum differentiation

Canonical momenta must remain independent symbols
until Poisson-bracket differentiation is complete.
Substitution of zero components before differentiation
invalidated an earlier multiplier-bracket calculation.

### CP098-C04: Longitudinal acceleration

The earlier proposed gamma coefficient depended
incorrectly on a special longitudinal acceleration.

The corrected expression is:

gamma = v^2 kappa (T.q)/2


### CP098-C05: Tangent-secondary bracket candidate

The earlier candidate containing p0.T and P^2 terms
was rejected by symbolic comparison.

Corrected regular-sector expression:

{psi_T,chi_T} = 4 A kappa^2/v

Derivation: homogeneity under u -> alpha*u,
with q and the canonical momenta held fixed.

Direct canonical derivative check at:
A=1, C=1, v=3, kappa=2

Result: 16/3, matching the corrected expression.

Source: tests/CP098_T14.py

Status: locally reproduced, not yet committed.
The full Dirac algorithm remains unresolved.


### CP098-C06: Correction to CP098-C01

CP098-C01 incorrectly identified a sign error by
confusing the two orientations of the Poisson bracket.

The independently reproduced identities are:

{chi_N,phi_N} = +2A/v^3
{phi_N,chi_N} = -2A/v^3

At v=4:

{chi_N,phi_N} = +A/32
{phi_N,chi_N} = -A/32

Therefore, the original checkpoint's positive sign
for {chi_N,phi_N} was correct.

CP098-C01 is superseded by this correction.

Evidence: tests/CP098_T15.py
Status: symbolic verification passed locally.
Git archival remains pending.


## CP098-T16–T23 extended audit

These results were reported from local terminal executions.
The scripts exist locally but have not been committed.
They are not yet independently reproduced from Git.

| ID | Calculation | Evidence | Archive |
|---|---|---|---|
| CP098-T16 | Tangent-tertiary closure requirement | ALGEBRAIC NECESSARY CONDITION | LOCAL FILE |
| CP098-T17 | Direct canonical bracket at one constrained point | REPORTED PASS | LOCAL FILE |
| CP098-T18 | Three constrained momentum samples | REPORTED PASS | LOCAL FILE |
| CP098-T19 | Symbolic two-parameter momentum family | REPORTED PASS | LOCAL FILE |
| CP098-T20 | Three geometry-dependent symbolic families | REPORTED PASS | LOCAL FILE |
| CP098-T21 | Five-parameter direct symbolic bracket | TIMEOUT | LOCAL FILE |
| CP098-T21A | Optimized directional derivative | REPORTED PASS; CANDIDATE REJECTED | LOCAL FILE |
| CP098-T21B | Reduced five-parameter full bracket | REPORTED PASS | LOCAL FILE |
| CP098-T22 | Transverse canonical contribution | DIAGNOSTIC MISMATCH | LOCAL FILE |
| CP098-T22A | Independent transverse-factor correction | REPORTED PASS | LOCAL FILE |
| CP098-T23 | Nine-parameter adapted-family bracket | REPORTED PASS | LOCAL FILE |

### CP098-T16

The tangent-tertiary consistency requirement was
derived algebraically using the exact Hamiltonian
decomposition and the previously tested bracket
{psi_T,chi_T}.

This is a necessary condition, not a general proof
of constraint preservation.

### CP098-T17 through T20

Exact canonical calculations passed on constrained
points and restricted symbolic families.

The tests did not establish unrestricted closure.

### CP098-T21

The direct five-parameter symbolic calculation
exceeded its 120-second time limit while calculating
the (u,p1) bracket.

Result: INCOMPLETE. No mathematical failure inferred.

### CP098-T21A and T21B

An optimized directional derivative agreed with
the previously tested T19 geometry.

An incorrect intermediate diagnostic candidate was
rejected in T21A.

T21B combined the directional derivative, reduced
(q,P) contribution, and internal canonical bracket.

The resulting five-parameter identity passed.

The T21B calculation reused the T21A contribution
and the previously reconstructed F_N expression;
it is not a wholly independent derivation.

### CP098-T22 and T22A

T22 found the transverse contribution:

{psi_T,chi_N}_{q,P}
    = kappa*v^2*(r*x+s*y)/(2*C)

The original diagnostic candidate incorrectly
contained v instead of v^2.

T22 therefore recorded a DIAGNOSTIC MISMATCH,
not a successful comparison with its candidate.

T22A independently checked the corrected v^2
coefficient and rejected the old candidate.

### CP098-T23

The constrained nine-parameter adapted family used:

v>0, kappa>0;
a,t,n,r,s,x,y real.

The toy internal specialization was:

A=I, C=1, B=p_I^2/2.

All five constraints vanished identically.

The complete tested canonical bracket matched
the T16 requirement with zero symbolic residual.

The transverse contribution was:

{psi_T,chi_N}_{q,P}
    = kappa*v^2*(r*x+s*y)/2

The calculation still depends on an adapted-frame
specialization and a previously reconstructed
normal-preservation expression F_N.

It is not unrestricted Dirac-Bergmann closure.

### Outstanding mathematical dependencies

1. Establish the general regular constraint-surface
   identities without the toy internal specialization.
2. Independently derive and verify F_N on that surface.
3. Incorporate both primary multipliers and all
   constraint-preservation equations.
4. Determine whether additional constraints arise.
5. Classify first-class and second-class constraints.
6. Establish the physical degree-of-freedom count.
7. Analyze A=0 crossings and the singular C=0 sector.
8. Investigate well-posedness and dynamical stability.
9. Test correspondence with established physics.
10. Seek independent physical and experimental evidence.

DYNAMICAL CONSISTENCY: OPEN.

## CP098-T24–T26 independent normal-preservation audit

Evidence: locally executed symbolic calculations.
Git archival and independent reproduction remain pending.

| ID | Calculation | Evidence | Archive |
|---|---|---|---|
| CP098-T24 | Independent normal-preservation bracket, nine-parameter toy family | REPORTED PASS | LOCAL FILE |
| CP098-T25 | General internal canonical contribution for two internal coordinates | REPORTED PASS | LOCAL FILE |
| CP098-T26 | Symbolic-coupling external normal-preservation bracket | REPORTED PASS | LOCAL FILE |

### Findings

T24 independently computed {chi_N,H_c} from the
canonical Hamiltonian on the nine-parameter adapted
family with A=I, C=1 and B=p_I^2/2.

T25 derived the internal contribution 2*kappa*D_A
for two internal coordinates with quadratic kinetic
terms and coordinate-dependent couplings.

T26 independently computed the external contribution
with symbolic A, B and nonzero C on an adapted
regular constrained family.

These calculations support the reconstructed F_N
expression within their stated domains.

They do not establish an unrestricted covariant
identity, complete constraint closure, constraint
classification or physical stability.

DYNAMICAL CONSISTENCY: OPEN.

## Open mathematical dependencies

- General constraint-surface bracket identities
- Complete Dirac-Bergmann closure
- First-class and second-class classification
- Physical degrees of freedom
- Singular C=0 sector
- Well-posedness and dynamical stability
- Correspondence with established physics
- Independent experimental validation

## Reproduction environment

Python 3.12.3
SymPy 1.14.0

## Publication status

This index records research progress.
Several executable tests exist as local files,
but the CP098 test collection has not yet been
committed or independently reproduced from Git.

Reported symbolic passes are not physical validation.

## CP098-T27 through T30 milestone

### T27 — Non-adapted normal bracket

Initial test failed because the supplied acceleration vector
did not have the expected curvature.

Original fixture preserved:
`tests/CP098_T27_initial_fixture.py.bak`

Correction:
The acceleration vector was changed to the exact
Lorentz transform of the intended rest-frame vector.

Corrected result:
- v = 4
- kappa = 1
- Primary constraints satisfied
- Direct (q,P) bracket = -64/3
- Predicted bracket = -64/3
- Difference = 0
- Normal-vector derivative tensor verified

### T28 — Boosted external bracket

- All five constraints satisfied
- Symbolic transverse momenta retained
- Direct and candidate (u,p1) brackets agree
- Difference = 0

### T29 — Symbolic regular family

- Symbolic v, kappa, tangential acceleration,
  A, B, C, and transverse momenta
- All five constraints satisfied
- Direct and candidate brackets agree
- Difference = 0

### T30 — Complete external momentum components

- Restored two transverse p0 components
- All five constraints satisfied
- Direct and candidate brackets agree
- Difference = 0

### Interpretation

The regular constrained-family calculation is compatible
with a general Lorentz-frame reduction argument.

The frame-exhaustiveness argument must be explicitly
documented and checked before treating the resulting
normal-preservation formula as a fully audited
regular-surface identity.

No unrestricted off-surface identity is claimed.

DYNAMICAL CONSISTENCY: OPEN.
GENERAL DIRAC CLOSURE: NOT ESTABLISHED.
PHYSICAL VALIDATION: NOT ESTABLISHED.

## CP098 T31-T39 REGULAR CLOSURE MILESTONE

Status: PASS for the specified regular-sector Hamiltonian model.

Regularity assumptions:
- v > 0
- kappa > 0
- A != 0
- C != 0
- Specified quadratic internal kinetic model
- Regular five-constraint surface and adapted Lorentz-frame coverage

Evidence:
- T31: regular constrained-family parameterization
- T32: external tertiary bracket
- T33: internal tertiary bracket
- T34: combined identity and normal multiplier
- T35: external tangential preservation
- T36: complete tangential cancellation
- T37: chi_T preservation
- T38R: primary preservation
- T39: five-equation consistency assembly

T38 development failures:
- Original derivative implementation failed.
- First correction left a second derivative defect.
- Unrestricted symbolic implementation timed out.
- T38R replaced it with efficient directional derivatives.
- Preserve failed versions and diagnostic logs.

Conclusion:
The five constraints have mutually consistent preservation
equations on the specified regular constraint surface.

The normal multiplier is fixed:
lambda_N = -v^3 F_N/(2A).

The tangential multiplier lambda_T remains undetermined
by these preservation equations.

Limitations:
- Constraint classification remains open.
- Gauge interpretation remains open.
- Singular sectors remain open.
- General physical validation remains open.
- T39 assembles earlier results; it is not an independent
  canonical derivation.

No experimental or consciousness-related claim follows
from this checkpoint.
