#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 087

CLASSICAL STATE EQUIVALENCE AND QUOTIENT STRUCTURE

Established inputs:
- Classical phase space R^2
- Free relativistic Hamiltonian
- Koopman inverse-flow pullback
- Classical multiplication observables

Investigated:
- Equivalence by identical densities
- Preservation of equivalence by dynamics
- Classical expectation-value equivalence
- Nonlinearity of the density map
- Why Hilbert-space addition does not
  descend to the classical state quotient

No quantum postulates are introduced.
"""

import sympy as sp

I = sp.I

x, p, t, k = sp.symbols(
    "x p t k", real=True
)

m, a, b = sp.symbols(
    "m a b", positive=True, real=True
)

H = sp.sqrt(p**2 + m**2)
v = sp.diff(H, p)

N = (a*b)**sp.Rational(1, 4)/sp.sqrt(sp.pi)

psi = N*sp.exp(
    -a*x**2/2-b*p**2/2
)

phi = sp.exp(I*k*x)*psi

print("="*64)
print("GWFM / SEFI — CHECKPOINT 087")
print("CLASSICAL STATE EQUIVALENCE AND QUOTIENT")
print("="*64)

def eq(A, B):
    return sp.simplify(A-B) == 0

def density(f):
    return sp.simplify(f*sp.conjugate(f))

def U(f, time):
    return f.subs(
        x, x-v*time,
        simultaneous=True
    )

def L(f):
    return -v*sp.diff(f, x)

def M(observable, f):
    return observable*f

def Dx(f):
    return -I*sp.diff(f, x)

def expectation(f, observable):
    return sp.integrate(
        observable*density(f),
        (x,-sp.oo,sp.oo),
        (p,-sp.oo,sp.oo)
    )

# 1. Normalized reference state.
assert eq(expectation(psi,1),1)

print("PASS: reference classical density normalized")

# 2. Distinct complex functions
# can represent identical densities.
assert eq(density(phi),density(psi))

assert not eq(
    phi.subs(k,1),
    psi
)

print("PASS: distinct Koopman functions")
print("      can represent the same density")

# 3. Equivalence is reflexive,
# symmetric and transitive for
# equality of densities.
#
# Demonstrate a three-state example.
theta = sp.exp(I*2*k*p)*psi

assert eq(density(psi),density(phi))
assert eq(density(phi),density(theta))
assert eq(density(psi),density(theta))

print("PASS: density equality defines")
print("      an equivalence relation")

# 4. Classical observables depend
# only on the density.
observables = (
    1, x, p, x**2, p**2, x*p
)

for f in observables:
    assert eq(
        f*density(psi),
        f*density(phi)
    )

print("PASS: all tested classical")
print("      observable densities agree")

# 5. Exact Gaussian moments.
assert eq(expectation(psi,x),0)
assert eq(expectation(psi,p),0)

assert eq(
    expectation(psi,x**2),
    1/(2*a)
)

assert eq(
    expectation(psi,p**2),
    1/(2*b)
)

print("PASS: exact classical expectations")

# 6. Koopman evolution preserves
# density equivalence.
psi_t = U(psi,t)
phi_t = U(phi,t)

assert eq(
    density(psi_t),
    density(phi_t)
)

print("PASS: density equivalence")
print("      preserved by Koopman evolution")

# 7. Density evolves by classical
# Liouville transport.
rho_t = density(psi_t)

assert eq(
    sp.diff(rho_t,t)
    +v*sp.diff(rho_t,x),
    0
)

print("PASS: quotient density obeys")
print("      classical Liouville equation")

# 8. Density evolution depends
# only on initial density.
rho_0 = density(psi)

assert eq(
    density(U(psi,t)),
    U(rho_0,t)
)

print("PASS: classical density dynamics")
print("      is well-defined on equivalence classes")

# 9. Invertibility of the flow.
assert eq(
    U(U(psi,t),-t),
    psi
)

assert eq(
    U(U(phi,t),-t),
    phi
)

print("PASS: inverse classical flow exists")

# 10. Multiplication observables
# preserve density equivalence
# when applied to equivalent states.
for f in (x,p,x+p,x**2+p**2):
    assert eq(
        density(M(f,psi)),
        density(M(f,phi))
    )

print("PASS: multiplication operators")
print("      preserve density equivalence")

# 11. Differential operators need not
# preserve equivalence.
#
# At k=1, Dx(phi) and Dx(psi)
# generally have different densities.

difference_Dx = sp.simplify(
    density(Dx(phi))-density(Dx(psi))
)

test_difference = sp.simplify(
    difference_Dx.subs({
        x:0,
        p:0,
        k:1
    })
)

assert test_difference != 0

print("PASS: derivative operators need not")
print("      preserve density equivalence")

# 12. Hilbert-space addition is not
# well-defined on density equivalence
# classes.
#
# psi and -psi represent the same
# classical density, but:
#
# psi+psi = 2psi
# psi+(-psi) = 0

assert eq(
    density(-psi),
    density(psi)
)

sum_same = density(psi+psi)
sum_opposite = density(psi-psi)

assert eq(
    sum_same,
    4*density(psi)
)

assert eq(sum_opposite,0)

assert not eq(
    sum_same,
    sum_opposite
)

print("PASS: Hilbert-space addition does not")
print("      descend to the density quotient")

# 13. Classical probability mixtures
# are defined directly on densities.
#
# A convex mixture is independent
# of phase representatives.

q = sp.Rational(1,3)

rho_A = density(psi)
rho_B = density(U(psi,1))

mixture = q*rho_A+(1-q)*rho_B

assert eq(
    mixture,
    q*density(phi)
    +(1-q)*density(U(phi,1))
)

print("PASS: classical probability mixtures")
print("      are phase-independent")

# 14. Mixtures evolve consistently.
assert eq(
    U(mixture,t),
    q*U(rho_A,t)
    +(1-q)*U(rho_B,t)
)

print("PASS: classical mixture evolution")
print("      is linear on densities")

print("\n"+"="*64)
print("CHECKPOINT 087: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Density equality defines a")
print("  classical state equivalence relation.")
print("- Classical multiplication observables")
print("  cannot distinguish equivalent states.")
print("- Koopman evolution preserves equivalence.")
print("- The density quotient follows Liouville")
print("  evolution independently of phase.")
print("- Differential operators need not preserve")
print("  classical density equivalence.")
print("- Hilbert-space vector addition does not")
print("  define addition on the quotient.")
print("- Classical convex probability mixtures")
print("  remain well-defined on densities.")

print("\nMATHEMATICAL QUALIFICATIONS:")
print("- Explicit tests use smooth Gaussians.")
print("- General density equivalence is understood")
print("  almost everywhere in L2.")
print("- The quotient is a classical state-space")
print("  construction, not a quantum derivation.")
print("- No claim of fundamental physical")
print("  gauge symmetry has been made.")

print("\nNOT ESTABLISHED:")
print("- Quantum superposition probabilities.")
print("- A noncommutative physical observable algebra.")
print("- A preferred action scale.")
print("- The quantum Born rule.")
print("- Quantum mechanics emerging from GWFM.")
