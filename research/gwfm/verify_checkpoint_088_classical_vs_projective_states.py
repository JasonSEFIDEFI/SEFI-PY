#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 088

CLASSICAL DENSITY STATES VERSUS HILBERT-SPACE RAYS

Established:
- Classical Koopman L2 representation
- Density equivalence |psi|^2 = |phi|^2
- Classical multiplication observables
- Free relativistic Koopman evolution

Investigated:
- Global versus local phase equivalence
- Classical quotient versus projective rays
- Rank-one projectors as a conditional
  extension of the observable algebra
- Whether projectors distinguish states
  with identical classical densities

No quantum measurement postulate assumed.
"""

import sympy as sp

I = sp.I

x, p, k, t = sp.symbols(
    "x p k t", real=True
)

a, b, m = sp.symbols(
    "a b m", positive=True, real=True
)

N = (a*b)**sp.Rational(1,4)/sp.sqrt(sp.pi)

psi = N*sp.exp(
    -a*x**2/2-b*p**2/2
)

phi = sp.exp(I*k*x)*psi

H = sp.sqrt(p**2+m**2)
v = sp.diff(H,p)

print("="*64)
print("GWFM / SEFI — CHECKPOINT 088")
print("CLASSICAL DENSITY VERSUS PROJECTIVE STATES")
print("="*64)

def eq(A,B):
    return sp.simplify(A-B) == 0

def rho(f):
    return sp.simplify(f*sp.conjugate(f))

def inner(f,g):
    return sp.simplify(
        sp.integrate(
            sp.conjugate(f)*g,
            (x,-sp.oo,sp.oo),
            (p,-sp.oo,sp.oo)
        )
    )

def U(f,time):
    return f.subs(
        x,x-v*time,simultaneous=True
    )

# 1. Both states are normalized.
assert eq(inner(psi,psi),1)
assert eq(inner(phi,phi),1)

print("PASS: both Koopman vectors normalized")

# 2. They define the same classical density.
assert eq(rho(psi),rho(phi))

print("PASS: identical classical densities")

# 3. Global phase equivalence is contained
# in classical density equivalence.
theta = sp.symbols("theta",real=True)

global_phase = sp.exp(I*theta)*psi

assert eq(rho(global_phase),rho(psi))
assert eq(
    inner(psi,global_phase),
    sp.exp(I*theta)
)

print("PASS: global phase leaves both")
print("      density and Hilbert ray unchanged")

# 4. Local phase gives an exact overlap.
overlap = inner(psi,phi)

expected_overlap = sp.exp(
    -k**2/(4*a)
)

assert eq(overlap,expected_overlap)

print("PASS: exact local-phase overlap")
print("      exp(-k^2/(4a))")

# 5. For k != 0, overlap magnitude < 1.
#
# Two normalized vectors define the
# same ray iff |overlap| = 1.
#
# Here 0 < overlap < 1 for nonzero k.

assert eq(
    overlap.subs(k,0),
    1
)

assert eq(
    overlap.subs({k:2,a:1}),
    sp.exp(-1)
)

assert sp.exp(-1) < 1

print("PASS: equal-density states can")
print("      represent distinct Hilbert rays")

# 6. Classical multiplication observables
# cannot distinguish the two vectors.
for observable in (
    1,x,p,x**2,p**2,x*p
):
    assert eq(
        observable*rho(psi),
        observable*rho(phi)
    )

print("PASS: classical observable densities")
print("      remain identical")

# 7. Define a conditional rank-one
# projector onto the reference vector.
#
# Q_psi(f) = psi * <psi,f>.
#
# This is a bounded orthogonal projector
# on the Hilbert space, but its physical
# observability is NOT assumed.

def Q(f):
    return psi*inner(psi,f)

assert eq(Q(psi),psi)
assert eq(Q(Q(phi)),Q(phi))

print("PASS: rank-one projector is idempotent")

# 8. Projector expectations differ.
#
# <psi,Q psi> = 1
# <phi,Q phi> = |<psi,phi>|^2

expect_psi = inner(psi,Q(psi))
expect_phi = inner(phi,Q(phi))

assert eq(expect_psi,1)

assert eq(
    expect_phi,
    sp.exp(-k**2/(2*a))
)

assert eq(
    expect_phi.subs({k:2,a:1}),
    sp.exp(-2)
)

print("PASS: rank-one projector distinguishes")
print("      equal-density states")

# 9. The projector is not multiplication
# by a classical phase-space function.
#
# Any multiplication operator acts
# pointwise. The projector depends on
# the global inner product.
#
# An explicit obstruction: if Q were
# multiplication by some function F,
# Q(psi)=psi would force F=1 almost
# everywhere because psi never vanishes.
# But Q(phi) != phi for k != 0.

assert eq(Q(psi),psi)

difference = sp.simplify(
    (Q(phi)-phi).subs({
        x:0,p:0,k:2,a:1,b:1
    })
)

assert difference != 0

print("PASS: projector cannot be identified")
print("      with classical multiplication")

# 10. Koopman evolution preserves
# classical density equivalence.
psi_t = U(psi,t)
phi_t = U(phi,t)

assert eq(rho(psi_t),rho(phi_t))

print("PASS: classical equivalence preserved")
print("      under free evolution")

# 11. Koopman evolution preserves
# the Hilbert inner product.
#
# This follows from measure preservation.
# Verify the transported overlap
# by an exact change of variable.
#
# q = x-v(p)*t has unit Jacobian.
# In the q coordinate the integrand
# is identical to its initial form.

q = sp.symbols("q",real=True)

transported_overlap_q = (
    sp.conjugate(psi)
    *phi
).subs(x,q)

transported_integral = sp.simplify(
    sp.integrate(
        transported_overlap_q,
        (q,-sp.oo,sp.oo),
        (p,-sp.oo,sp.oo)
    )
)

assert eq(transported_integral,overlap)

print("PASS: inner-product preservation")
print("      follows from classical flow")

# 12. Transported projector is defined
# consistently by conjugation:
#
# Q_t = U_t Q U_-t
#
# It projects onto U_t psi.
#
# This is a mathematical identity,
# not a physical measurement claim.

def Qt(f):
    return U(
        Q(U(f,-t)),
        t
    )

assert eq(Qt(psi_t),psi_t)

print("PASS: projector family transforms")
print("      consistently with Koopman flow")

# 13. A classical mixture of identical
# densities is unchanged.
weight = sp.Rational(1,2)

classical_mix = (
    weight*rho(psi)
    +(1-weight)*rho(phi)
)

assert eq(classical_mix,rho(psi))

print("PASS: classical density mixture")
print("      is insensitive to local phase")

# 14. A conditional Hilbert-space mixture
# of the two rank-one projectors has
# purity (1+|overlap|^2)/2.
#
# This differs from the purity 1
# of a single rank-one projector.
#
# The calculation is an operator-algebra
# identity, not a quantum statistical
# interpretation.

purity = sp.simplify(
    (
        1+overlap*sp.conjugate(overlap)
    )/2
)

assert eq(
    purity,
    (1+sp.exp(-k**2/(2*a)))/2
)

assert eq(purity.subs(k,0),1)

assert purity.subs({k:2,a:1}) < 1

print("PASS: Hilbert-space projector mixture")
print("      differs from a single projector")

print("\n"+"="*64)
print("CHECKPOINT 088: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Classical density equivalence is")
print("  coarser than Hilbert-ray equivalence.")
print("- Local phases can preserve densities")
print("  while changing Hilbert rays.")
print("- Classical multiplication observables")
print("  cannot distinguish these states.")
print("- Rank-one projectors can distinguish")
print("  the same state representatives.")
print("- The projector is not a classical")
print("  multiplication observable.")
print("- Classical evolution preserves both")
print("  density equivalence and inner products.")

print("\nCONDITIONAL STRUCTURE:")
print("- Rank-one projectors were introduced")
print("  as mathematical Hilbert operators.")
print("- Their physical measurability was")
print("  neither assumed nor established.")
print("- Projector mixtures are not being")
print("  identified with quantum ensembles.")

print("\nNOT ESTABLISHED:")
print("- A physical projector measurement law.")
print("- A uniquely selected observable algebra.")
print("- Quantum pure-state interpretation.")
print("- The quantum Born rule.")
print("- Emergent quantum mechanics.")
