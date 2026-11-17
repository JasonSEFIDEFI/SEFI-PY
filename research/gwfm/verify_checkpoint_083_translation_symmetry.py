#!/usr/bin/env python3
"""
GWFM / SEFI — Checkpoint 083

Translation symmetry versus physical momentum.

Established classical model:
Phase space R^2 with coordinates x,p.
Poisson bracket {x,p}=1.
Free relativistic Hamiltonian H(p).

Conditional representation:
Complex L2 functions on phase space.

Goal:
Determine whether spatial translation symmetry
uniquely identifies quantum momentum or an
action scale.

No quantum emergence is assumed.
"""

import sympy as sp

x, p, a, b = sp.symbols(
    "x p a b", real=True
)
m, alpha, beta, kappa = sp.symbols(
    "m alpha beta kappa",
    positive=True, real=True
)

I = sp.I

print("=" * 64)
print("GWFM / SEFI — CHECKPOINT 083")
print("TRANSLATION SYMMETRY AND MOMENTUM")
print("=" * 64)

def equal(A, B):
    return sp.simplify(A-B) == 0

def PB(f, g):
    return sp.simplify(
        sp.diff(f,x)*sp.diff(g,p)
        -sp.diff(f,p)*sp.diff(g,x)
    )

def X(f):
    return x*f

def P(f):
    return p*f

def Dx(f):
    return -I*sp.diff(f,x)

def comm(A, B, f):
    return sp.simplify(A(B(f))-B(A(f)))

psi = sp.exp(
    -alpha*x**2-beta*p**2
)

# 1. Classical momentum generates
# spatial translations through PB.
assert PB(x,p) == 1
assert PB(p,p) == 0

print("PASS: classical momentum generates")
print("      canonical spatial translations")

# 2. Spatial translation of phase-space
# functions by inverse-flow pullback.
def T(f, distance):
    return f.subs(
        x, x-distance,
        simultaneous=True
    )

assert equal(T(T(psi,a),b),T(psi,a+b))
assert equal(T(psi,0),psi)

print("PASS: translation group composition")

# 3. Its infinitesimal generator is
# -d/dx = -i*Dx.
assert equal(
    sp.diff(T(psi,a),a).subs(a,0),
    -I*Dx(psi)
)

print("PASS: translation generator is Dx=-i*d/dx")

# 4. Classical momentum multiplication
# is distinct from Dx.
assert equal(comm(X,P,psi),0)
assert equal(comm(X,Dx,psi),I*psi)

assert not equal(P(psi),Dx(psi))

print("PASS: classical momentum and")
print("      translation generator are distinct")

# 5. Translation covariance of X and P.
# T_-a X T_a = X+a.
def conjugate_T(A,f,distance):
    return T(A(T(f,distance)),-distance)

assert equal(
    conjugate_T(X,psi,a),
    (x+a)*psi
)

assert equal(
    conjugate_T(P,psi,a),
    p*psi
)

print("PASS: X shifts and P remains invariant")

# 6. Free relativistic dynamics
# commutes with spatial translations.
H = sp.sqrt(p**2+m**2)
v = sp.diff(H,p)

def L(f):
    return -v*sp.diff(f,x)

assert equal(
    T(L(psi),a),
    L(T(psi,a))
)

print("PASS: translation symmetry of")
print("      free relativistic dynamics")

# 7. Candidate momentum operator.
def Pk(f):
    return P(f)+kappa*Dx(f)

assert equal(
    comm(X,Pk,psi),
    I*kappa*psi
)

print("PASS: candidate commutator scale")
print("      remains kappa")

# 8. Candidate translation family:
#
# V_a = exp(-i*a*P_kappa/kappa)
#
# Since P and Dx commute:
#
# V_a psi = exp(-i*a*p/kappa)
#           psi(x-a,p).
#
# This differs from T_a by a phase.

def V(f,distance):
    return (
        sp.exp(-I*distance*p/kappa)
        *T(f,distance)
    )

assert equal(
    V(V(psi,a),b),
    V(psi,a+b)
)

assert equal(V(psi,0),psi)

print("PASS: candidate unitary translation")
print("      family satisfies group law")

# 9. Infinitesimal generator of V.
assert equal(
    sp.diff(V(psi,a),a).subs(a,0),
    -I*Pk(psi)/kappa
)

print("PASS: candidate momentum generates V")

# 10. Candidate translation covariance.
def conjugate_V(A,f,distance):
    return V(A(V(f,distance)),-distance)

assert equal(
    conjugate_V(X,psi,a),
    (x+a)*psi
)

assert equal(
    conjugate_V(P,psi,a),
    p*psi
)

print("PASS: V implements the same")
print("      classical coordinate translation")

# 11. V and T differ on complex states.
ratio = sp.simplify(V(psi,a)/T(psi,a))

assert equal(
    ratio,
    sp.exp(-I*a*p/kappa)
)

assert sp.simplify(
    ratio.subs({
        a:sp.pi,
        p:1,
        kappa:1
    })
) == -1

print("PASS: translation implementations")
print("      differ by momentum-dependent phase")

# 12. Both produce identical
# classical phase-space densities.
rho_T = T(psi,a)*sp.conjugate(T(psi,a))
rho_V = V(psi,a)*sp.conjugate(V(psi,a))

assert equal(rho_T,rho_V)

print("PASS: both induce the same")
print("      classical density translation")

# 13. The candidate representation
# commutes with free transport.
assert equal(
    L(V(psi,a)),
    V(L(psi),a)
)

print("PASS: candidate translations")
print("      commute with free dynamics")

# 14. Different positive kappa values
# satisfy the same tested covariance.
for K in (sp.Integer(1),
          sp.Integer(2),
          sp.Rational(1,2)):

    Vk = V(psi,a).subs(kappa,K)

    assert equal(
        Vk*sp.conjugate(Vk),
        rho_T
    )

    assert equal(
        comm(
            X,
            lambda f,K=K: P(f)+K*Dx(f),
            psi
        ),
        I*K*psi
    )

print("PASS: multiple action scales")
print("      satisfy classical covariance")

# 15. Requiring exactly T rather than V
# is an additional representation choice.
#
# V has a p-dependent phase for finite
# nonzero kappa. Classical densities
# cannot distinguish the two actions.

print("\n" + "=" * 64)
print("CHECKPOINT 083: ALL TESTS PASSED")
print("=" * 64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Classical p generates translations")
print("  through the Poisson bracket.")
print("- Dx generates translations on L2.")
print("- Multiplication by p differs from Dx.")
print("- T and V implement the same coordinate")
print("  translation on classical densities.")
print("- Both commute with free relativistic flow.")
print("- Multiple kappa values remain admissible.")

print("\nCONDITIONAL:")
print("- V introduces a momentum-dependent phase.")
print("- Interpreting P_kappa as physical momentum")
print("  requires an additional identification.")

print("\nNOT ESTABLISHED:")
print("- A preferred nonzero kappa.")
print("- Quantum momentum from GWFM.")
print("- Physical quantum amplitudes.")
print("- The Born rule.")
print("- Emergent quantum mechanics.")
