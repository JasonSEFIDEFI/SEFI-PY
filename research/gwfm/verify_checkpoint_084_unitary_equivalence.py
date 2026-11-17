#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 084

Unitary equivalence of classical
phase-space translation representations.

Assumptions:
- Classical phase space R^2
- Complex L2 representation
- Real nonzero kappa
- Free relativistic Hamiltonian

No quantum physical interpretation assumed.
"""

import sympy as sp

x, p, a = sp.symbols("x p a", real=True)
m, alpha, beta = sp.symbols(
    "m alpha beta", positive=True
)
kappa = sp.symbols(
    "kappa", real=True, nonzero=True
)

I = sp.I
psi = sp.exp(-alpha*x**2-beta*p**2)

H = sp.sqrt(p**2+m**2)
v = sp.diff(H,p)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 084")
print("UNITARY EQUIVALENCE AND REPRESENTATION FREEDOM")
print("=" * 64)

def eq(A,B):
    return sp.simplify(A-B) == 0

def W(f):
    return sp.exp(-I*x*p/kappa)*f

def Winv(f):
    return sp.exp(I*x*p/kappa)*f

def T(f,d):
    return f.subs(x,x-d,simultaneous=True)

def V(f,d):
    return sp.exp(-I*d*p/kappa)*T(f,d)

def Dx(f):
    return -I*sp.diff(f,x)

def P(f):
    return p*f

def Pk(f):
    return P(f)+kappa*Dx(f)

def L(f):
    return -v*sp.diff(f,x)

# 1. W has a pointwise inverse.
assert eq(W(Winv(psi)),psi)
assert eq(Winv(W(psi)),psi)

print("PASS: W has the stated inverse")

# 2. W preserves pointwise probability density.
rho = psi*sp.conjugate(psi)
rho_W = W(psi)*sp.conjugate(W(psi))

assert eq(rho_W,rho)

print("PASS: W preserves classical density")

# 3. Because |W|=1, it preserves L2 norms.
# This is a pointwise proof, not merely
# a numerical integration check.
assert eq(
    sp.conjugate(
        sp.exp(-I*x*p/kappa)
    )*sp.exp(-I*x*p/kappa),
    1
)

print("PASS: W is unitary on L2")

# 4. Translation conjugation identity.
assert eq(
    W(T(Winv(psi),a)),
    V(psi,a)
)

print("PASS: V_a = W T_a W_inverse")

# 5. Generator conjugation identity.
assert eq(
    kappa*W(Dx(Winv(psi))),
    Pk(psi)
)

print("PASS: P_kappa = kappa W Dx W_inverse")

# 6. Distinct nonzero kappa values
# admit the same construction.
for K in (
    sp.Integer(1),
    sp.Integer(2),
    -sp.Integer(1),
    sp.Rational(1,2)
):
    lhs = W(T(Winv(psi),a)).subs(
        kappa,K
    )
    rhs = V(psi,a).subs(kappa,K)

    assert eq(lhs,rhs)

print("PASS: equivalence exists for multiple kappa")

# 7. Classical multiplication observables
# commute with W.
assert eq(W(x*Winv(psi)),x*psi)
assert eq(W(p*Winv(psi)),p*psi)

print("PASS: classical X and P unchanged")
print("      by unitary conjugation")

# 8. The classical transport generator
# DOES change under W.
#
# L_transformed = W L W_inverse.
#
# For the free relativistic flow:
#
# L_transformed = L - i*v*p/kappa.

L_transformed = W(L(Winv(psi)))

expected = (
    L(psi)-I*v*p*psi/kappa
)

assert eq(L_transformed,expected)

print("PASS: transformed classical generator")
print("      includes a phase term")

# 9. This extra term is nonzero
# away from p=0.
difference = sp.simplify(
    L_transformed-L(psi)
)

assert eq(
    difference,
    -I*p**2*psi/
    (kappa*sp.sqrt(p**2+m**2))
)

assert not eq(
    difference.subs({
        p:1,
        x:0,
        kappa:1,
        m:1
    }),
    0
)

print("PASS: W does not preserve the")
print("      original generator L")

# 10. The transformed evolution is
# consistent with the transformed state.
#
# If phi = W psi, then
# L_transformed phi = W L psi.

phi = W(psi)

def Lt(f):
    return (
        L(f)-I*v*p*f/kappa
    )

assert eq(
    Lt(phi),
    W(L(psi))
)

print("PASS: transformed states and")
print("      dynamics remain consistent")

# 11. The transformed generator still
# commutes with the transformed
# translation representation.
#
# Since L commutes with T, conjugation
# preserves that commutation.

def comm_transformed(f):
    return sp.simplify(
        Lt(V(f,a))-V(Lt(f),a)
    )

assert eq(comm_transformed(psi),0)

print("PASS: transformed dynamics remains")
print("      translation symmetric")

# 12. W is not a symmetry of the
# original L when L is held fixed.
assert not eq(
    W(L(psi)),
    L(W(psi))
)

print("PASS: representation change is")
print("      distinguished from a symmetry of L")

print("\n" + "=" * 64)
print("CHECKPOINT 084: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED MATHEMATICALLY:")
print("- W is unitary for real nonzero kappa.")
print("- V and T are unitarily equivalent.")
print("- Their translation generators are")
print("  related by unitary conjugation.")
print("- Classical multiplication X and P")
print("  are unchanged by W.")
print("- The classical transport generator")
print("  acquires an additional phase term.")
print("- Consistent transformation of states,")
print("  operators and dynamics preserves")
print("  the translated mathematical structure.")

print("\nIMPORTANT DISTINCTION:")
print("- W is a representation change.")
print("- W is not a symmetry of the fixed")
print("  original Koopman generator L.")
print("- The result applies for kappa != 0.")
print("- It does not identify a physical")
print("  quantum momentum operator.")

print("\nNOT ESTABLISHED:")
print("- A unique physical action scale.")
print("- Quantum dynamics from GWFM.")
print("- A quantum measurement theory.")
print("- The Born rule.")
print("- Emergent quantum mechanics.")
