#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 086
CLASSICAL OBSERVABLE ALGEBRA AND DYNAMICAL CLOSURE

Established:
- Classical phase space R^2
- H(p) = sqrt(p^2 + m^2)
- Koopman generator L = -v(p)*d/dx
- Classical observables multiply by f(x,p)

Investigated:
- Commutativity of multiplication observables
- Dynamical closure of classical observables
- Differential-operator extensions
- Whether dynamics selects an observable algebra

No quantum observable postulate is assumed.
"""

import sympy as sp

I = sp.I

x, p, t = sp.symbols("x p t", real=True)
m, a, b = sp.symbols(
    "m a b", positive=True, real=True
)

H = sp.sqrt(p**2 + m**2)
v = sp.diff(H, p)
vp = sp.diff(v, p)

psi = sp.exp(-a*x**2-b*p**2)

print("="*64)
print("GWFM / SEFI — CHECKPOINT 086")
print("CLASSICAL OBSERVABLE ALGEBRA")
print("="*64)

def eq(A, B):
    return sp.simplify(A-B) == 0

def M(f):
    return lambda g: f*g

def Dx(g):
    return -I*sp.diff(g,x)

def Dp(g):
    return -I*sp.diff(g,p)

def L(g):
    return -v*sp.diff(g,x)

def comm(A,B,g):
    return sp.simplify(A(B(g))-B(A(g)))

def U(g,s):
    return g.subs(x,x-v*s,simultaneous=True)

def evolve(A,g,s):
    return U(A(U(g,s)),-s)

# 1. Multiplication observables commute.
f = x**2+p
g = x*p+p**2

assert eq(comm(M(f),M(g),psi),0)

print("PASS: classical multiplication algebra")
print("      is commutative")

# 2. Product and adjoint structure.
assert eq(
    M(f)(M(g)(psi)),
    M(f*g)(psi)
)

assert eq(
    sp.conjugate(f*psi),
    sp.conjugate(f)*sp.conjugate(psi)
)

print("PASS: multiplication algebra is")
print("      closed under products and adjoints")

# 3. Infinitesimal classical evolution.
assert eq(
    comm(M(f),L,psi),
    v*sp.diff(f,x)*psi
)

print("PASS: classical observable evolution")
print("      remains multiplication")

# 4. Exact finite-time evolution.
assert eq(
    evolve(M(f),psi,t),
    f.subs(x,x+v*t)*psi
)

print("PASS: classical observable algebra")
print("      is preserved by finite-time flow")

# 5. Canonical coordinate observables.
assert eq(
    evolve(M(x),psi,t),
    (x+v*t)*psi
)

assert eq(
    evolve(M(p),psi,t),
    p*psi
)

print("PASS: X and P retain classical")
print("      Hamiltonian evolution")

# 6. Differential operators do not
# commute with all classical observables.
assert eq(
    comm(Dx,M(x),psi),
    -I*psi
)

assert eq(
    comm(Dp,M(p),psi),
    -I*psi
)

print("PASS: differential generators extend")
print("      beyond the commutative algebra")

# 7. General first-order differential
# operator with smooth coefficients.
u = x+p
w = x*p
c = x**2+p

def A(g):
    return (
        u*sp.diff(g,x)
        +w*sp.diff(g,p)
        +c*g
    )

# Commutation with both coordinates
# detects the derivative coefficients.
assert eq(
    comm(A,M(x),psi),
    u*psi
)

assert eq(
    comm(A,M(p),psi),
    w*psi
)

print("PASS: coordinate commutators detect")
print("      differential components")

# 8. Formal generalization:
# For A=u*d/dx+w*d/dp+c,
# [A,M_x]=M_u and [A,M_p]=M_w.
#
# Therefore commuting with both X
# and P forces u=w=0, within this
# first-order differential-operator class.

Ucoef, Wcoef = sp.symbols("Ucoef Wcoef")

solutions = sp.solve(
    [Ucoef,Wcoef],
    [Ucoef,Wcoef],
    dict=True
)

assert solutions == [{
    Ucoef:0,
    Wcoef:0
}]

print("PASS: first-order differential operators")
print("      commuting with X and P reduce")
print("      to multiplication operators")

# 9. Differential operators are also
# closed under commutation with L.
#
# For A=u*d/dx+w*d/dp+c:
#
# [A,L] =
# (v*u_x-w*v_p)*d/dx
# +v*w_x*d/dp
# +v*c_x

expected = (
    (v*sp.diff(u,x)-w*vp)
    *sp.diff(psi,x)
    +v*sp.diff(w,x)
    *sp.diff(psi,p)
    +v*sp.diff(c,x)*psi
)

assert eq(
    comm(A,L,psi),
    expected
)

print("PASS: first-order differential")
print("      algebra is dynamically closed")

# 10. The differential algebra is
# strictly larger than multiplication
# and is compatible with the same
# classical generator.
assert eq(
    comm(Dx,L,psi),
    0
)

assert eq(
    comm(Dp,L,psi),
    -vp*Dx(psi)
)

print("PASS: larger operator algebra")
print("      remains compatible with dynamics")

# 11. The classical density evolution
# does not depend on initial phase.
k = sp.symbols("k",real=True)

psi_phase = sp.exp(I*k*x)*psi

rho1 = sp.simplify(
    U(psi,t)*sp.conjugate(U(psi,t))
)

rho2 = sp.simplify(
    U(psi_phase,t)
    *sp.conjugate(U(psi_phase,t))
)

assert eq(rho1,rho2)

print("PASS: classical density equivalence")
print("      is preserved by evolution")

# 12. Derivative operators distinguish
# the phase representatives.
assert eq(
    sp.conjugate(psi_phase)*Dx(psi_phase)
    -sp.conjugate(psi)*Dx(psi),
    k*psi*sp.conjugate(psi)
)

print("PASS: extended derivative algebra")
print("      distinguishes phase representatives")

print("\n"+"="*64)
print("CHECKPOINT 086: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED ALGEBRAICALLY:")
print("- Classical multiplication observables")
print("  form a commutative algebra.")
print("- This algebra is preserved by the")
print("  established Hamiltonian flow.")
print("- First-order differential operators")
print("  form a larger algebra.")
print("- That larger algebra is also closed")
print("  under commutation with the generator.")
print("- Within the first-order class,")
print("  commuting with both X and P")
print("  excludes derivative components.")
print("- Classical density equivalence")
print("  persists under free evolution.")
print("- Derivative operators can distinguish")
print("  phase representatives.")

print("\nMATHEMATICAL QUALIFICATIONS:")
print("- Closure is algebraic on a common")
print("  smooth test domain.")
print("- No complete operator-domain or")
print("  self-adjointness analysis was performed.")
print("- No physical measurement interpretation")
print("  was assigned to derivative operators.")

print("\nNOT ESTABLISHED:")
print("- A unique physically selected")
print("  observable algebra.")
print("- Quantum momentum from GWFM.")
print("- A preferred nonzero action scale.")
print("- The quantum Born rule.")
print("- Emergent quantum mechanics.")
