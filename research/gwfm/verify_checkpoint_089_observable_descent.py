#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 089

OBSERVABLE DESCENT TO THE CLASSICAL STATE QUOTIENT

Established:
- Classical states represented by densities
- Koopman complex L2 representation
- Classical multiplication observables
- Density equivalence of phase representatives

Conditional:
- Differential operators
- Rank-one Hilbert-space projectors

Question:
Which operator expectations are well-defined
on classical density-equivalence classes?

No quantum measurement postulate is assumed.
"""

import sympy as sp

I = sp.I

x, p, k = sp.symbols("x p k", real=True)
a, b = sp.symbols("a b", positive=True)

N = (a*b)**sp.Rational(1,4)/sp.sqrt(sp.pi)

psi = N*sp.exp(-a*x**2/2-b*p**2/2)
phi = sp.exp(I*k*x)*psi

print("="*64)
print("GWFM / SEFI — CHECKPOINT 089")
print("OBSERVABLE DESCENT TO CLASSICAL STATE QUOTIENT")
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

def Dx(f):
    return -I*sp.diff(f,x)

def Q(f):
    return psi*inner(psi,f)

def W(f):
    return sp.exp(I*k*x)*f

def Winv(f):
    return sp.exp(-I*k*x)*f

def Q_transformed(f):
    return W(Q(Winv(f)))

# 1. Normalization.
assert eq(inner(psi,psi),1)
assert eq(inner(phi,phi),1)

print("PASS: both representatives normalized")

# 2. Classical equivalence.
assert eq(rho(psi),rho(phi))

print("PASS: phase-related states have")
print("      identical classical density")

# 3. Classical multiplication observables
# give the same expectation integrands.
for f in (1,x,p,x**2,p**2,x*p):
    assert eq(
        f*rho(psi),
        f*rho(phi)
    )

print("PASS: multiplication expectations")
print("      descend to density quotient")

# 4. Exact moments.
assert eq(inner(psi,x**2*psi),1/(2*a))
assert eq(inner(phi,x**2*phi),1/(2*a))

assert eq(inner(psi,p**2*psi),1/(2*b))
assert eq(inner(phi,p**2*phi),1/(2*b))

print("PASS: exact representative-independent")
print("      classical second moments")

# 5. A differential operator distinguishes
# these classically equivalent states.
assert eq(inner(psi,Dx(psi)),0)
assert eq(inner(phi,Dx(phi)),k)

print("PASS: Dx expectation depends")
print("      on phase representative")

# 6. Rank-one projector onto psi.
assert eq(Q(psi),psi)
assert eq(Q(Q(phi)),Q(phi))

print("PASS: rank-one projector is idempotent")

# 7. The fixed projector distinguishes
# representatives of the same density.
assert eq(
    inner(psi,Q(psi)),
    1
)

assert eq(
    inner(phi,Q(phi)),
    sp.exp(-k**2/(2*a))
)

print("PASS: fixed projector expectation")
print("      depends on phase representative")

# 8. For a specific nonzero phase gradient
# the expectation values differ.
value = sp.exp(-2)

assert eq(
    inner(phi,Q(phi)).subs({k:2,a:1}),
    value
)

assert value != 1

print("PASS: explicit counterexample")
print("      to projector descent")

# 9. A simultaneous unitary change of
# representation transforms the projector.
assert eq(
    Q_transformed(phi),
    phi
)

assert eq(
    inner(phi,Q_transformed(phi)),
    1
)

print("PASS: simultaneous state/projector")
print("      transformation preserves expectation")

# 10. Fixed versus transformed projector
# are distinct operations.
assert not eq(
    Q(phi).subs({k:2,a:1,b:1}),
    Q_transformed(phi).subs({k:2,a:1,b:1})
)

print("PASS: changing the state alone")
print("      differs from conjugating the observable")

# 11. The local phase also transforms Dx.
#
# W Dx W_inverse = Dx - k.
#
# This is a representation transformation,
# not a new physical gauge law.
assert eq(
    W(Dx(Winv(psi))),
    Dx(psi)-k*psi
)

print("PASS: derivative operator transforms")
print("      with the phase gradient")

# 12. Simultaneous transformation preserves
# the derivative expectation.
def Dx_transformed(f):
    return W(Dx(Winv(f)))

assert eq(
    inner(phi,Dx_transformed(phi)),
    inner(psi,Dx(psi))
)

print("PASS: simultaneous transformation")
print("      preserves derivative expectation")

# 13. The density map is invariant under W.
assert eq(rho(W(psi)),rho(psi))

print("PASS: local phase leaves classical")
print("      density map unchanged")

print("\n"+"="*64)
print("CHECKPOINT 089: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Classical multiplication expectations")
print("  are independent of phase representatives.")
print("- Dx and the tested rank-one projector")
print("  fail this quotient compatibility test.")
print("- The rank-one projector therefore cannot")
print("  be assigned a representative-independent")
print("  expectation on the classical quotient.")
print("- Simultaneously transforming states and")
print("  operators preserves Hilbert expectations.")
print("- Such a representation transformation")
print("  differs from changing only the state.")

print("\nSCOPE LIMITATIONS:")
print("- Counterexamples establish failure for")
print("  the specific operators tested.")
print("- No full classification of all bounded")
print("  quotient-compatible operators is claimed.")
print("- Differential identities use smooth states.")
print("- Projector measurability is not assumed.")
print("- No fundamental gauge principle is assumed.")

print("\nNOT ESTABLISHED:")
print("- A unique physical observable algebra.")
print("- A quantum state interpretation.")
print("- A physical projector measurement law.")
print("- A preferred nonzero action scale.")
print("- The quantum Born rule.")
print("- Emergent quantum mechanics.")
