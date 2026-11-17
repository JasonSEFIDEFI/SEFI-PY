#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 082

Does classical relativistic dynamics select
a noncommuting operator representation?

Established:
H(p) = sqrt(p**2 + m**2)
L = -v(p)*d/dx

Conditional candidate operators:
X_eta = X + eta*Dp
P_kappa = P + kappa*Dx

No quantum interpretation is assumed.
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

H = sp.sqrt(p**2 + m**2)
v = sp.diff(H, p)
vp = sp.diff(v, p)

psi = sp.exp(-a*x**2 - b*p**2)

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 082")
print("DYNAMICAL SELECTION OF OPERATORS")
print("=" * 64)

def equal(A, B):
    return sp.simplify(A-B) == 0

def X(f):
    return x*f

def P(f):
    return p*f

def Dx(f):
    return -I*sp.diff(f, x)

def Dp(f):
    return -I*sp.diff(f, p)

def L(f):
    return -v*sp.diff(f, x)

def Xeta(f):
    return X(f) + eta*Dp(f)

def Pkappa(f):
    return P(f) + kappa*Dx(f)

def comm(A, B, f):
    return sp.simplify(A(B(f))-B(A(f)))

def U(f, time):
    return f.subs(
        x, x-v*time, simultaneous=True
    )

def evolved_operator(A, f):
    # A_t = U_-t A U_t
    return U(A(U(f, t)), -t)

# 1. Verify the established transport equation.
assert equal(
    sp.diff(U(psi, t), t),
    L(U(psi, t))
)
print("PASS: Koopman evolution equation")

# 2. Classical position evolution.
assert equal(
    comm(X, L, psi),
    v*psi
)

assert equal(
    evolved_operator(X, psi),
    (x+v*t)*psi
)

print("PASS: classical position evolves as x+v*t")

# 3. Classical momentum conservation.
assert equal(comm(P, L, psi), 0)

assert equal(
    evolved_operator(P, psi),
    p*psi
)

print("PASS: classical momentum is conserved")

# 4. Spatial translation generator.
assert equal(comm(Dx, L, psi), 0)

print("PASS: Dx commutes with classical evolution")

# 5. Momentum translation generator.
assert equal(
    comm(Dp, L, psi),
    -vp*Dx(psi)
)

print("PASS: Dp has a velocity-gradient correction")

# 6. Candidate momentum conservation.
assert equal(
    comm(Pkappa, L, psi),
    0
)

assert equal(
    evolved_operator(Pkappa, psi),
    Pkappa(psi)
)

print("PASS: P_kappa conserved for arbitrary kappa")

# 7. Candidate position evolution.
expected_rate = v*psi - eta*vp*Dx(psi)

assert equal(
    comm(Xeta, L, psi),
    expected_rate
)

print("PASS: X_eta has an eta-dependent correction")

# 8. Exact finite-time evolution.
expected_finite = (
    Xeta(psi) + t*expected_rate
)

assert equal(
    evolved_operator(Xeta, psi),
    expected_finite
)

print("PASS: exact finite-time X_eta evolution")

# 9. Velocity gradient.
assert equal(
    vp,
    m**2/(p**2+m**2)**sp.Rational(3, 2)
)

assert equal(
    vp.subs(p, 0),
    1/m
)

print("PASS: velocity gradient is nonzero")

# 10. Conditional classical position law.
correction = sp.simplify(
    comm(Xeta, L, psi)-v*psi
)

assert equal(correction.subs(eta, 0), 0)

coefficient = sp.simplify(
    sp.diff(correction, eta).subs({
        x: 1,
        p: 0
    })
)

assert coefficient != 0

print("PASS: imposing classical position")
print("      evolution selects eta=0")

# 11. Classical evolution leaves kappa free.
for K in (0, 1, 2, sp.Rational(1, 2)):
    def candidate(f, K=K):
        return P(f)+K*Dx(f)

    assert equal(
        comm(candidate, L, psi),
        0
    )

    assert equal(
        comm(X, candidate, psi),
        I*K*psi
    )

print("PASS: dynamical compatibility")
print("      does not select kappa")

# 12. Exact classical momentum identification.
difference = sp.simplify(
    Pkappa(psi)-P(psi)
)

assert equal(
    difference.subs(kappa, 0),
    0
)

assert not equal(
    difference.subs(kappa, 1),
    0
)

print("PASS: identification with classical")
print("      momentum requires kappa=0")

print("\n" + "=" * 64)
print("CHECKPOINT 082: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED ALGEBRAICALLY:")
print("- Classical X and P obey Hamiltonian motion.")
print("- Dx is conserved under free evolution.")
print("- P_kappa is conserved for every kappa.")
print("- X_eta has an eta-dependent correction.")
print("- Requiring the classical position law")
print("  selects eta=0.")
print("- Momentum conservation leaves kappa free.")

print("\nCONDITIONAL IDENTIFICATION:")
print("- Requiring P_kappa to equal classical")
print("  momentum multiplication selects kappa=0.")

print("\nNOT ESTABLISHED:")
print("- A preferred nonzero commutator scale.")
print("- Quantum operators selected by GWFM.")
print("- A physical quantum state space.")
print("- The Born rule.")
print("- Emergent quantum mechanics.")
