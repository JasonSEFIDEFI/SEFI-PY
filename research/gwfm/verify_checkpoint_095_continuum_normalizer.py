#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 095

CONTINUUM DENSITY-QUOTIENT NORMALIZERS

Concrete model:
H = L2(R2, dx dp)

Coordinate transformation:
T_c(x,p) = (c*x,p), c > 0.

Unitary pullback:
U_c psi(x,p) = sqrt(c)*psi(c*x,p).

The Jacobian factor makes this unitary
even though T_c does not preserve
phase-space volume when c != 1.

Scope:
Exact symbolic verification of one
family of continuum transformations.

No general normalizer theorem is
claimed from these computations alone.
"""

import sympy as sp

I = sp.I

x, p = sp.symbols("x p", real=True)
c, d = sp.symbols("c d", positive=True)
a, b = sp.symbols("a b", positive=True)
k, ell = sp.symbols("k ell", real=True)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 095")
print("CONTINUUM DENSITY-QUOTIENT NORMALIZERS")
print("=" * 64)

def eq(left, right):
    difference = sp.expand_complex(left-right)
    difference = sp.trigsimp(difference)
    return sp.simplify(difference) == 0

def check(left, right):
    assert eq(left, right)

def U(f, scale):
    return sp.sqrt(scale)*f.subs(
        x, scale*x, simultaneous=True
    )

def Uinv(f, scale):
    return f.subs(
        x, x/scale, simultaneous=True
    )/sp.sqrt(scale)

def rho(f):
    return sp.simplify(f*sp.conjugate(f))

def M(obs, f):
    return obs*f

N = (a*b)**sp.Rational(1,4)/sp.sqrt(sp.pi)

psi = N*sp.exp(
    -a*x**2/2-b*p**2/2
)

phi = sp.exp(I*(k*x+ell*p))*psi

# 1. Classical density equivalence.
check(rho(psi), rho(phi))

print("PASS: local phases preserve density")

# 2. Coordinate map and inverse.
check(
    (c*x).subs(x,x/c),
    x
)

print("PASS: coordinate transformation invertible")

# 3. Jacobian determinant.
J = sp.Matrix([
    [c,0],
    [0,1]
])

check(J.det(),c)

print("PASS: coordinate Jacobian calculated")

# 4. Unitary representation has
# the required Jacobian factor.
check(
    rho(U(psi,c)),
    c*rho(psi).subs(x,c*x)
)

print("PASS: transformed density includes Jacobian")

# 5. Exact normalization by Gaussian
# integration in the continuum.
norm = sp.integrate(
    rho(psi),
    (x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

norm_transformed = sp.integrate(
    rho(U(psi,c)),
    (x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

check(norm,1)
check(norm_transformed,1)

print("PASS: original and transformed states normalized")

# 6. Inverse representation.
check(Uinv(U(psi,c),c),psi)
check(U(Uinv(psi,c),c),psi)

print("PASS: operator inverse verified")

# 7. Group composition.
check(
    U(U(psi,c),d),
    U(psi,c*d)
)

print("PASS: scaling operators form a group")

# 8. Density equivalence descends
# through the transformation.
check(
    rho(U(psi,c)),
    rho(U(phi,c))
)

print("PASS: transformed equivalent states")
print("      remain density-equivalent")

# 9. Observable conjugation.
f = x**2+x*p+2*p**2

left = Uinv(
    M(f,U(psi,c)),c
)

right = M(
    f.subs(x,x/c),psi
)

check(left,right)

print("PASS: multiplication observables")
print("      remain multiplication operators")

# 10. Conjugation preserves products.
g = x+3*p

check(
    (f*g).subs(x,x/c),
    f.subs(x,x/c)*g.subs(x,x/c)
)

print("PASS: observable products preserved")

# 11. Conjugation preserves adjoints.
h = f+I*g

check(
    sp.conjugate(h).subs(x,x/c),
    sp.conjugate(h.subs(x,x/c))
)

print("PASS: observable adjoints preserved")

# 12. Local phase operators remain
# local phase operators.
chi = k*x+ell*p

left = Uinv(
    sp.exp(I*chi)*U(psi,c),
    c
)

right = sp.exp(
    I*chi.subs(x,x/c)
)*psi

check(left,right)

print("PASS: local phase algebra normalized")

# 13. Position transforms.
check(
    Uinv(x*U(psi,c),c),
    x*psi/c
)

print("PASS: position multiplication")
print("      transforms by inverse scaling")

# 14. Momentum multiplication is
# unchanged by this coordinate map.
check(
    Uinv(p*U(psi,c),c),
    p*psi
)

print("PASS: momentum multiplication unchanged")

# 15. The map is generally not
# volume preserving.
assert not eq(
    J.det().subs(c,2),1
)

print("PASS: example includes")
print("      non-volume-preserving maps")

# 16. Omitting the Jacobian factor
# fails normalization.
bare = psi.subs(x,c*x)

bare_norm = sp.integrate(
    rho(bare),
    (x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

check(bare_norm,1/c)

assert not eq(
    bare_norm.subs(c,2),1
)

print("PASS: Jacobian amplitude is necessary")
print("      for this unitary pullback")

# 17. Arbitrary phase after scaling
# cannot change output density.
eta = k*x**2+ell*p

check(
    rho(sp.exp(I*eta)*U(psi,c)),
    rho(U(psi,c))
)

print("PASS: additional local phase")
print("      preserves transformed density")

# 18. The induced probability map
# depends only on input density.
r = rho(psi)

check(
    rho(U(psi,c)),
    c*r.subs(x,c*x)
)

print("PASS: continuum transformation")
print("      descends to density quotient")

print("\n"+"="*64)
print("CHECKPOINT 095: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED IN THIS MODEL:")
print("- Reversible coordinate scaling admits")
print("  a Jacobian-corrected unitary pullback.")
print("- It preserves total L2 probability.")
print("- It preserves density equivalence.")
print("- It normalizes multiplication and")
print("  local phase operator algebras.")
print("- Its induced density transformation")
print("  is independent of local phase.")
print("- Volume preservation is not necessary")
print("  for the corrected unitary pullback.")

print("\nCONTINUUM THEOREM CONTEXT:")
print("- On standard sigma-finite measure spaces,")
print("  normalizers of the multiplication algebra")
print("  are weighted composition unitaries,")
print("  modulo measure-theoretic identifications.")
print("- General transformations require")
print("  measurable invertibility, measure-class")
print("  preservation, and Radon-Nikodym weights.")
print("- These symbolic checks establish only")
print("  the concrete scaling example.")
print("- They do not prove the general theorem.")

print("\nNOT ESTABLISHED:")
print("- A physically selected transformation.")
print("- A unique continuum normalizer.")
print("- A new interaction law.")
print("- Quantum measurement or Born rule.")
print("- A preferred action scale.")
print("- Emergent quantum mechanics.")
