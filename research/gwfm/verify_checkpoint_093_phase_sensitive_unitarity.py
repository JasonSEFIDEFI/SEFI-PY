#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 093

PROBABILITY CONSERVATION AND
PHASE-SENSITIVE TRANSFORMATIONS

Established baseline:
- Classical probability distributions
- Stochastic observation maps
- Classical density equivalence

Conditional extension:
- Complex amplitudes
- Norm-preserving unitary transformation
- Output probabilities defined by modulus squared

Purpose:
Determine whether norm conservation alone
excludes or permits phase-sensitive outcomes.

No physical quantum interpretation is assumed.
"""

import sympy as sp

I = sp.I
theta = sp.symbols("theta", real=True)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 093")
print("PROBABILITY CONSERVATION AND PHASE SENSITIVITY")
print("=" * 64)

def eq(a, b):
    if isinstance(a, sp.MatrixBase):
        return all(
            sp.simplify(z) == 0
            for z in a-b
        )
    return sp.simplify(a-b) == 0

def density(z):
    return sp.Matrix([
        sp.simplify(z[j]*sp.conjugate(z[j]))
        for j in range(z.rows)
    ])

def total_probability(rho):
    return sp.simplify(sum(rho))

# 1. Normalized classical input.
rho = sp.Matrix([
    sp.Rational(1,3),
    sp.Rational(2,3)
])

assert eq(total_probability(rho), 1)

print("PASS: classical input normalized")

# 2. Complex representatives differing
# only by relative phase.
psi = sp.Matrix([
    sp.sqrt(sp.Rational(1,3)),
    sp.sqrt(sp.Rational(2,3))*sp.exp(I*theta)
])

assert eq(density(psi), rho)

print("PASS: relative phase leaves")
print("      input density unchanged")

# 3. Classical stochastic mixing.
K = sp.Matrix([
    [sp.Rational(1,2), sp.Rational(1,2)],
    [sp.Rational(1,2), sp.Rational(1,2)]
])

for j in range(2):
    assert eq(sum(K[:,j]), 1)

assert all(entry >= 0 for entry in K)

print("PASS: classical mixing matrix")
print("      is column-stochastic")

classical_output = K*rho

assert eq(
    classical_output,
    sp.Matrix([
        sp.Rational(1,2),
        sp.Rational(1,2)
    ])
)

print("PASS: classical output independent")
print("      of relative phase")

# 4. Conditional unitary mixing.
#
# The matrix below is a mathematical
# unitary operator. No physical
# beam splitter is assumed.

H = sp.Matrix([
    [1, 1],
    [1,-1]
])/sp.sqrt(2)

assert eq(
    H.conjugate().T*H,
    sp.eye(2)
)

print("PASS: amplitude mixing operator")
print("      is exactly unitary")

# 5. Apply unitary transformation.
amplitude_output = H*psi
coherent_output = density(amplitude_output)

cross = sp.sqrt(2)*sp.cos(theta)/3

expected = sp.Matrix([
    sp.Rational(1,2)+cross,
    sp.Rational(1,2)-cross
])

assert eq(coherent_output, expected)

print("PASS: coherent output contains")
print("      exact relative-phase cross terms")

# 6. Total probability conserved.
assert eq(
    total_probability(coherent_output),
    1
)

print("PASS: unitary amplitude mixing")
print("      preserves total probability")

# 7. Output probabilities are
# nonnegative because they are
# squared complex moduli.
#
# Also verify explicit extrema.
for phase_value in (0, sp.pi/2, sp.pi):
    result = [
        sp.simplify(z.subs(theta,phase_value))
        for z in coherent_output
    ]
    assert all(
        sp.simplify(z).is_nonnegative
        for z in result
    )
    assert eq(sum(result),1)

print("PASS: tested coherent outputs")
print("      are nonnegative and normalized")

# 8. Different phases produce
# different output distributions.
output_zero = coherent_output.subs(theta,0)
output_pi = coherent_output.subs(theta,sp.pi)

assert not eq(output_zero,output_pi)

print("PASS: phase-sensitive outcomes")
print("      coexist with norm conservation")

# 9. Classical inputs remain identical
# at these two phases.
assert eq(
    density(psi.subs(theta,0)),
    density(psi.subs(theta,sp.pi))
)

print("PASS: distinguishable coherent outputs")
print("      arise from identical input densities")

# 10. No fixed classical stochastic
# matrix acting only on rho can
# reproduce both phase-dependent
# output distributions.
#
# Both inputs have the same rho.
# Every function of rho alone must
# assign the same output to both.
#
# The distinct outputs give an
# explicit contradiction.

assert eq(K*density(psi.subs(theta,0)),
          K*density(psi.subs(theta,sp.pi)))

assert not eq(output_zero,output_pi)

print("PASS: fixed density-only stochastic")
print("      maps cannot reproduce both outputs")

# 11. Dephasing / phase averaging.
#
# Uniformly average the conditional
# output over the relative phase.

averaged = sp.Matrix([
    sp.simplify(
        sp.integrate(
            coherent_output[j],
            (theta,0,2*sp.pi)
        )/(2*sp.pi)
    )
    for j in range(2)
])

assert eq(averaged,classical_output)

print("PASS: uniform phase averaging")
print("      recovers classical mixing")

# 12. The classical mixing matrix
# is the entrywise squared modulus
# of the chosen unitary matrix.

K_from_H = sp.Matrix([
    [
        sp.simplify(
            H[i,j]*sp.conjugate(H[i,j])
        )
        for j in range(2)
    ]
    for i in range(2)
])

assert eq(K_from_H,K)

print("PASS: classical transition matrix")
print("      equals squared unitary entries")

# 13. For a single populated input
# component, the output has no
# relative-phase interference.
basis0 = sp.Matrix([1,0])
basis1 = sp.Matrix([0,1])

assert eq(
    density(H*basis0),
    K*basis0
)

assert eq(
    density(H*basis1),
    K*basis1
)

print("PASS: basis-state probabilities")
print("      agree with classical mixing")

# 14. For coherent two-component
# input, the classical and conditional
# outputs differ at a fixed phase.

assert not eq(
    coherent_output.subs(theta,0),
    classical_output
)

print("PASS: coherent mixing differs")
print("      from classical mixing")

# 15. The conditional unitary does
# not define a map on the classical
# density quotient.
#
# Equivalent inputs have different
# output densities, violating the
# descent criterion of Checkpoint 089.

assert eq(
    density(psi.subs(theta,0)),
    density(psi.subs(theta,sp.pi))
)

assert not eq(
    density(H*psi.subs(theta,0)),
    density(H*psi.subs(theta,sp.pi))
)

print("PASS: conditional unitary mixing")
print("      fails classical quotient descent")

print("\n" + "=" * 64)
print("CHECKPOINT 093: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Classical stochastic mixing preserves")
print("  probability and ignores local phase.")
print("- The tested unitary transformation")
print("  preserves the complex-vector norm.")
print("- Modulus-squared outputs remain")
print("  normalized and phase-sensitive.")
print("- Identical classical input densities")
print("  can produce different coherent outputs.")
print("- No fixed density-only stochastic map")
print("  reproduces both phase-dependent outputs.")
print("- Uniform phase averaging recovers")
print("  the corresponding classical mixing.")
print("- The tested coherent transformation")
print("  does not descend to the classical")
print("  density-equivalence quotient.")

print("\nMATHEMATICAL QUALIFICATIONS:")
print("- The example is finite-dimensional.")
print("- The unitary and output rule were")
print("  explicitly introduced as conditional")
print("  mathematical structures.")
print("- Norm conservation does not itself")
print("  select a physical state space.")
print("- Phase averaging is a mathematical")
print("  operation, not a derived decoherence law.")
print("- No physical apparatus was modeled.")

print("\nNOT ESTABLISHED:")
print("- A physical coherent mixing mechanism.")
print("- Quantum measurement from GWFM.")
print("- A derivation of the Born rule.")
print("- A preferred action scale.")
print("- Emergent quantum mechanics.")
