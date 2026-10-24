# Five-minute checkpoint 01: which bridge is actually missing?

23 October 2026. Scope: source tracing and three exact-algebra checks, not a
complete review of all SEFI manuscripts or a new field theory.

## Source correction and map

There is a SEFI/DEFI action. The missing step is not absence of every action;
it is a derived relationship between the worldline/internal-coordinate action
and the bulk two-complex-scalar action used for the localized candidate.

| Source inspected | Mathematical content | Status relative to the candidate |
| --- | --- | --- |
| Public SEFI Formal Mathematical Construction, section 6.2 | Scalar sum of four layer functionals along a worldline | Does not alone introduce four independent dynamical fields |
| Public DEFI Formal Construction, introduction A | Explicit extension to four independent internal coordinates and a coupled variational model | A genuine additional modeling postulate; not merely a change of notation |
| Local Foundations of a Geometric Identity Theory draft, 17 October 2026, action/coupling and field-extension sections | Proper-time worldline action, four internal oscillators, curvature/torsion couplings, distributional occupation construction | Provides a concrete classical model; does not derive the bulk candidate action or quantum field amplitudes |
| Existing scalar candidate derivation | Two complex spacetime fields, canonical kinetic terms, independent phase symmetries, quartic potential | The fields, coefficients and symmetry choice remain assumptions relative to the inspected worldline framework |
| Existing source-coupling working note | An explicitly additional scalar interaction for source/momentum accounting | Already labeled a candidate addition, not a derivation of this two-field soliton theory |

Public sources:
[SEFI manuscript](https://github.com/JasonSEFIDEFI/PhD/blob/main/SEFI_Formal_Mathematical_Construction.pdf),
[DEFI manuscript](https://github.com/JasonSEFIDEFI/PhD/blob/main/DEFI_Formal_Construction.pdf).
Local drafts were inspected as LaTeX at `output/theory-foundations/foundations.tex`
and `output/sefi-source-coupling/sefi-source-coupling.tex` in the shared
workspace, not republished or assumed to be accepted literature. The public
PDFs were read by text extraction; no layout or visual-equation fidelity audit
was performed. The conclusions below have independent explicit algebra.

## Exact result 1: symmetry permits a family, not the selected constants

Conditionally assume two complex scalar fields in 3+1 dimensions, independent
global U(1) phases, and a local potential polynomial through fourth field
degree. The complete nonderivative invariant basis is

    1, rho, eta, rho^2, rho eta, eta^2,
    rho=|Phi|^2, eta=|Sigma|^2.

Independent phases force equal powers of each field and its conjugate. Direct
enumeration in `bridge_check.py` checks the complete degree-four list. This
does not classify higher-order or derivative EFT operators, and the symmetry
itself is an assumption to be justified from SEFI.

Write a useful positive-parameter family as

    U = a rho^2/4 + b(eta-v^2)^2/4 + c rho/2 + g rho eta/2,
    a,b,c,g,v^2 > 0.

It has vacuum rho=0, eta=v^2 and m_Phi^2=c+g v^2. Exact polynomial expansion
gives

    U - m_Phi^2 rho/2
      = b/4 (eta-v^2+g rho/b)^2 + (ab-g^2)rho^2/(4b).

Thus ab>=g^2 obstructs below-free-Phi-threshold binding as in the prior audit.
For g^2>ab, sufficiently small positive rho with
eta=v^2-g rho/b >=0 makes this difference negative: a pointwise subthreshold
well exists. This remains permission for binding, not a solution or stability
proof. The earlier parameters are a=c=1, g=4, b=lambda, v^2=3/lambda.
Both lambda=12 and lambda=25 satisfy the same phase symmetry; they lie on
opposite sides of ab-g^2=0. The symmetry therefore cannot by itself select
the binding regime. A SEFI derivation must restrict parameters beyond it.

## Exact result 2: four real coordinates do not supply the required symmetry

The tempting pairing Phi=I_O+i I_A, Sigma=I_S+i I_W matches a count of real
variables. It does not establish matching dynamics, spacetime dependence,
kinetic coefficients or independent phase symmetries.

One allowed example from the foundations action has all internal stiffnesses
one, constant mu1=1, prescribed curvature kappa=1 and zero torsion:

    V_internal = (I_O^2+I_A^2+I_S^2+I_W^2)/2 + I_S.

At (0,0,1,0), V=3/2. The proposed Sigma phase rotation by pi gives
(0,0,-1,0), where V=-1/2. The bulk portal potential would be invariant under
that rotation. Therefore the direct pairing does not generally identify the
two models. This counterexample does not rule out every coarse-graining,
restricted sector, emergent symmetry or alternative bridge.

## Concrete next five-minute task

Test the simplest proposed bridge before adding new fields: turn off the
curvature/torsion couplings and ask which equalities among the internal kinetic
weights and stiffnesses are required for rotations within (O,A) and (S,W).
Then restore the existing couplings and identify which symmetries they break.
Record those equalities as extra assumptions unless the premises force them.

After that, choose a clearly labeled bridge hypothesis: an ensemble/coarse
graining of worldlines with a specified measure, or a new bulk-field postulate.
Simply replacing I(s) by I(x) changes the model and does not derive spatial
gradient terms. No choice is made silently here. This checkpoint narrows the
missing derivation and provides a falsifiable starting calculation.

## Validation

`test_bridge.py` checks the complete invariant monomial list, the exact
rational polynomial identity with both parameter regimes, and the explicit
failed symmetry pairing. Actual test output and machine-readable algebra are
recorded alongside this note. No new particle, quantum or empirical claim is
made.

Actual local outcome: all nine tests passed in 4.03 seconds (three bridge
checks plus the six prior action/dynamics checks). `bridge-results.json`
contains the exact-algebra outputs. The source manifest records the examined
public PDF and local draft hashes, and the captured test summary. The remote
workflow now includes the bridge checks; its outcome is not inferred from the
local run. Earlier scientific output files are unchanged.
