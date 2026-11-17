#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 092

CLASSICAL COARSE-GRAINING VERSUS
RELATIVE-PHASE INTERFERENCE

Established inputs:
- Classical probability distributions
- Deterministic, invertible evolution
- Density-equivalent complex representatives
- Classical coarse-grained measurements

Conditional comparison:
- Coherent addition of complex amplitudes

Purpose:
Determine whether ordinary classical
coarse-graining introduces dependence on
relative phase.

No quantum measurement rule is assumed.
"""

import sympy as sp

I = sp.I
theta = sp.symbols("theta", real=True)

print("="*64)
print("GWFM / SEFI — CHECKPOINT 092")
print("CLASSICAL COARSE-GRAINING AND INTERFERENCE")
print("="*64)

def eq(a,b):
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

# Three classical microscopic states.
rho = sp.Matrix([
    sp.Rational(1,2),
    sp.Rational(1,3),
    sp.Rational(1,6)
])

assert eq(sum(rho),1)

print("PASS: classical microscopic")
print("      distribution normalized")

# A deterministic invertible evolution.
# This permutation preserves total
# probability and the discrete measure.
P = sp.Matrix([
    [0,1,0],
    [0,0,1],
    [1,0,0]
])

assert eq(P.T*P,sp.eye(3))
assert eq(sum(P*rho),1)

print("PASS: deterministic classical")
print("      evolution preserves probability")

# Classical detector with two bins.
# First bin collects microscopic
# states 0 and 1; second collects 2.
C = sp.Matrix([
    [1,1,0],
    [0,0,1]
])

detected = C*P*rho

assert eq(
    detected,
    sp.Matrix([
        sp.Rational(1,2),
        sp.Rational(1,2)
    ])
)

assert eq(sum(detected),1)

print("PASS: classical coarse-graining")
print("      produces normalized outcomes")

# Complex Koopman representatives.
psi = sp.Matrix([
    sp.sqrt(rho[j]) for j in range(3)
])

D = sp.diag(
    1,
    sp.exp(I*theta),
    sp.exp(2*I*theta)
)

phi = D*psi

assert eq(density(psi),rho)
assert eq(density(phi),rho)

print("PASS: distinct local phases")
print("      preserve microscopic density")

# The classical evolution preserves
# density equivalence.
assert eq(
    density(P*psi),
    density(P*phi)
)

print("PASS: classical evolution preserves")
print("      density equivalence")

# Coarse-grained measurements depend
# only on microscopic probabilities.
assert eq(
    C*density(P*psi),
    C*density(P*phi)
)

assert eq(
    C*density(P*phi),
    detected
)

print("PASS: classical detector outcomes")
print("      are independent of phase")

# Check several relative phases.
for value in (
    0,
    sp.pi/3,
    sp.pi/2,
    sp.pi
):
    assert eq(
        (C*density(P*phi)).subs(theta,value),
        detected
    )

print("PASS: multiple relative phases")
print("      produce identical classical outcomes")

# Classical coarse-graining is linear
# on probability distributions.
rho2 = sp.Matrix([
    sp.Rational(1,6),
    sp.Rational(1,6),
    sp.Rational(2,3)
])

weight = sp.Rational(2,5)

mixture = weight*rho+(1-weight)*rho2

assert eq(
    C*P*mixture,
    weight*C*P*rho
    +(1-weight)*C*P*rho2
)

print("PASS: classical probability mixtures")
print("      commute with coarse-graining")

# Conditional comparison:
#
# Instead of adding probabilities,
# add complex amplitudes in each bin.
#
# This is a DIFFERENT mathematical map.
# It is not the classical detector rule.
z = P*phi

coherent_amplitudes = sp.Matrix([
    z[0]+z[1],
    z[2]
])

coherent_intensities = density(
    coherent_amplitudes
)

expected_first = (
    rho[1]+rho[2]
    +2*sp.sqrt(rho[1]*rho[2])
    *sp.cos(theta)
)

assert eq(
    coherent_intensities[0],
    expected_first
)

print("PASS: conditional amplitude addition")
print("      produces a cross term")

# The cross term depends on phase.
at_zero = sp.simplify(
    coherent_intensities[0].subs(theta,0)
)

at_pi = sp.simplify(
    coherent_intensities[0].subs(theta,sp.pi)
)

assert not eq(at_zero,at_pi)

print("PASS: conditional interference term")
print("      depends on relative phase")

# The classical probability sum does not
# contain that cross term.
classical_first = detected[0]

assert eq(
    coherent_intensities[0]
    -classical_first,
    2*sp.sqrt(rho[1]*rho[2])
    *sp.cos(theta)
)

print("PASS: classical probability addition")
print("      contains no interference term")

# Amplitude addition is not generally
# probability preserving without an
# additional normalization or structure.
coherent_total = sp.simplify(
    sum(coherent_intensities)
)

expected_total = (
    1
    +2*sp.sqrt(rho[1]*rho[2])
    *sp.cos(theta)
)

assert eq(coherent_total,expected_total)

assert not eq(
    coherent_total.subs(theta,0),
    1
)

print("PASS: coherent bin summation is not")
print("      automatically probability preserving")

# General finite-dimensional statement:
#
# For a nonnegative, column-stochastic
# classical observation matrix K,
# observed probabilities are K*rho.
#
# Local phases do not enter rho.
#
# Therefore no such fixed classical
# observation map can distinguish
# density-equivalent representatives.

u,v,w = sp.symbols(
    "u v w",
    real=True,
    nonnegative=True
)

generic_rho = sp.Matrix([u,v,w])

assert eq(
    sum(C*generic_rho),
    sum(generic_rho)
)

print("PASS: probability-preserving")
print("      classical observation map")

print("\n"+"="*64)
print("CHECKPOINT 092: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Classical deterministic evolution")
print("  preserves probability distributions.")
print("- Classical coarse-graining acts")
print("  on probabilities, not amplitudes.")
print("- Density-equivalent representatives")
print("  give identical classical outcomes.")
print("- Classical mixtures remain linear")
print("  under coarse-graining.")
print("- Coherent amplitude summation")
print("  introduces a relative-phase cross term.")
print("- That amplitude rule is not implied")
print("  by classical coarse-graining.")
print("- The tested coherent bin sum is not")
print("  automatically probability preserving.")

print("\nSCOPE LIMITATIONS:")
print("- This is an exact finite-state model.")
print("- The general conclusion applies to")
print("  fixed classical observation maps")
print("  that depend only on probability density.")
print("- It does not rule out phase-sensitive")
print("  models with additional structure.")
print("- The amplitude comparison is conditional.")
print("- No physical detector model was derived.")

print("\nNOT ESTABLISHED:")
print("- A physical amplitude-combination rule.")
print("- Quantum interference emerging from GWFM.")
print("- A quantum measurement law.")
print("- The Born rule.")
print("- A preferred nonzero action scale.")
print("- Emergent quantum mechanics.")
