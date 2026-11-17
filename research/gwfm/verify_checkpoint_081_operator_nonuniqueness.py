#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 081

Operator non-uniqueness in the classical
Koopman Hilbert-space representation.

Classical phase space: R^2 with (x,p).

Test:
- Classical multiplication operators
- Translation generators
- Noncommuting candidate operators
- Freedom in their commutator scale
- Classical relativistic transport generator

No physical quantum interpretation is assumed.
"""

import sympy as sp

x, p, t = sp.symbols("x p t", real=True)
m, a, b = sp.symbols(
    "m a b", positive=True, real=True
)
kappa, eta = sp.symbols(
    "kappa eta", real=True
)

I = sp.I

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 081")
print("CLASSICAL HILBERT-SPACE OPERATOR NON-UNIQUENESS")
print("=" * 64)

# A smooth, rapidly decreasing test function.
psi = sp.exp(-a*x**2-b*p**2)

def eq(A, B):
    return sp.simplify(sp.expand(A-B)) == 0

def comm(A, B, f):
    return sp.simplify(A(B(f))-B(A(f)))

# Classical multiplication observables.
def X(f):
    return x*f

def P(f):
    return p*f

# Translation generators on phase space.
def Dx(f):
    return -I*sp.diff(f,x)

def Dp(f):
    return -I*sp.diff(f,p)

# 1. Classical observables commute.
assert eq(comm(X,P,psi),0)

print("PASS: classical [X,P] = 0")

# 2. Translation generators have
# nonzero commutators with coordinates.
assert eq(comm(X,Dx,psi),I*psi)
assert eq(comm(P,Dp,psi),I*psi)

print("PASS: [X,Dx] = i*Identity")
print("PASS: [P,Dp] = i*Identity")

# 3. Cross commutators vanish.
assert eq(comm(X,Dp,psi),0)
assert eq(comm(P,Dx,psi),0)
assert eq(comm(Dx,Dp,psi),0)

print("PASS: independent phase-space")
print("      translation directions commute")

# 4. Candidate momentum family:
#
# P_kappa = P + kappa*Dx.
#
# This is a mathematical operator family,
# not an identification of physical
# quantum momentum.

def P_candidate(f):
    return P(f)+kappa*Dx(f)

assert eq(
    comm(X,P_candidate,psi),
    I*kappa*psi
)

print("PASS: [X,P_kappa] = i*kappa*Identity")

# 5. Every real kappa satisfies the
# corresponding commutator identity.
for value in (0,1,2,sp.Rational(1,2)):
    lhs = comm(
        X,
        lambda f: P(f)+value*Dx(f),
        psi
    )
    assert eq(lhs,I*value*psi)

print("PASS: multiple commutator scales admitted")

# 6. Candidate position family:
#
# X_eta = X + eta*Dp.
#
# Combined candidates have
# [X_eta,P_kappa] = i*(kappa-eta)*I.

def X_candidate(f):
    return X(f)+eta*Dp(f)

assert eq(
    comm(X_candidate,P_candidate,psi),
    I*(kappa-eta)*psi
)

print("PASS: two-parameter operator family")

# 7. Distinct choices can produce the
# same commutator scale.
#
# (kappa,eta) = (1,0) and (2,1)
# both give i*Identity.

pair_A = sp.simplify(
    (kappa-eta).subs({kappa:1,eta:0})
)

pair_B = sp.simplify(
    (kappa-eta).subs({kappa:2,eta:1})
)

assert pair_A == pair_B == 1

print("PASS: commutator scale does not")
print("      uniquely select operators")

# 8. The established classical Hamiltonian
# determines a transport generator.
H = sp.sqrt(p**2+m**2)
v = sp.diff(H,p)

def L(f):
    return -v*sp.diff(f,x)

# Its action is fixed by the classical flow.
assert eq(
    L(psi),
    2*a*x*v*psi
)

print("PASS: classical transport generator recovered")

# 9. The classical generator does not
# contain kappa or eta.
assert sp.diff(L(psi),kappa) == 0
assert sp.diff(L(psi),eta) == 0

print("PASS: classical transport is independent")
print("      of candidate commutator scales")

# 10. The physical classical momentum
# observable remains multiplication by p.
#
# The candidate P+kappa*Dx is a different
# operator unless kappa=0.

difference = sp.simplify(
    P_candidate(psi)-P(psi)
)

assert eq(
    difference,
    2*I*a*kappa*x*psi
)

assert eq(difference.subs(kappa,0),0)
assert not eq(difference.subs(kappa,1),0)

print("PASS: candidate momentum differs")
print("      from classical momentum")

print("\n" + "=" * 64)
print("CHECKPOINT 081: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED ALGEBRAICALLY:")
print("- Classical X and P multiplication commute.")
print("- Translation generators are noncommuting")
print("  with their corresponding coordinates.")
print("- Candidate operators can realize")
print("  [X,P_kappa] = i*kappa*Identity.")
print("- The parameter kappa remains arbitrary.")
print("- Different operator choices can produce")
print("  the same commutator scale.")
print("- Classical transport does not fix kappa.")

print("\nMATHEMATICAL QUALIFICATIONS:")
print("- Identities were checked on a Gaussian")
print("  test function; the differential-operator")
print("  identities extend to suitable common")
print("  smooth domains.")
print("- Self-adjointness and physical observable")
print("  interpretation were not established.")

print("\nNOT ESTABLISHED:")
print("- A physically preferred quantum momentum.")
print("- A uniquely determined Planck constant.")
print("- Quantum dynamics from GWFM.")
print("- The Born rule.")
print("- Emergent quantum mechanics.")
