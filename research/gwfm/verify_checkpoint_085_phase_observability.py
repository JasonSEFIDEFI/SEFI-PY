#!/usr/bin/env python3
"""
GWFM / SEFI — CHECKPOINT 085

LOCAL PHASE FREEDOM AND CLASSICAL OBSERVABILITY

Established input:
- Classical phase space R^2
- Free relativistic Hamiltonian
- Koopman transport on complex L2
- Classical multiplication observables

Conditional input:
- Real smooth phase chi(x,p)
- Unitary multiplication W=exp(i*chi)

No physical gauge principle or quantum
measurement postulate is assumed.
"""

import sympy as sp

I = sp.I

x, p = sp.symbols("x p", real=True)
m, alpha, beta = sp.symbols(
    "m alpha beta", positive=True
)
c, k, r = sp.symbols(
    "c k r", real=True
)

H = sp.sqrt(p**2+m**2)
v = sp.diff(H,p)

chi = c*x*p+k*x+r*p**2

N = (alpha*beta)**sp.Rational(1,4)/sp.sqrt(sp.pi)

psi = N*sp.exp(
    -alpha*x**2/2-beta*p**2/2
)

print("="*64)
print("GWFM / SEFI — CHECKPOINT 085")
print("LOCAL PHASE FREEDOM AND CLASSICAL OBSERVABILITY")
print("="*64)

def eq(A,B):
    return sp.simplify(A-B) == 0

def W(f):
    return sp.exp(I*chi)*f

def Winv(f):
    return sp.exp(-I*chi)*f

def Dx(f):
    return -I*sp.diff(f,x)

def Dp(f):
    return -I*sp.diff(f,p)

def L(f):
    return -v*sp.diff(f,x)

def M(f,g):
    return f*g

def rho(f):
    return f*sp.conjugate(f)

# 1. Unitary inverse.
assert eq(W(Winv(psi)),psi)
assert eq(Winv(W(psi)),psi)

print("PASS: local phase operator has inverse")

# 2. Pointwise density preservation.
assert eq(rho(W(psi)),rho(psi))

print("PASS: arbitrary tested real phase")
print("      preserves classical density")

# 3. Classical observables are unchanged.
for observable in (
    x,p,x**2,p**2,x*p,x**2+p**2
):
    transformed = W(
        M(observable,Winv(psi))
    )

    assert eq(
        transformed,
        M(observable,psi)
    )

print("PASS: classical multiplication")
print("      observables remain invariant")

# 4. Differential operators transform
# with gradients of the phase.
assert eq(
    W(Dx(Winv(psi))),
    Dx(psi)-sp.diff(chi,x)*psi
)

assert eq(
    W(Dp(Winv(psi))),
    Dp(psi)-sp.diff(chi,p)*psi
)

print("PASS: differential operators acquire")
print("      phase-gradient corrections")

# 5. Classical transport transforms.
transformed_L = W(L(Winv(psi)))

expected_L = (
    L(psi)+I*v*sp.diff(chi,x)*psi
)

assert eq(transformed_L,expected_L)

print("PASS: transport generator transforms")
print("      with a phase-gradient term")

# 6. Consistency of transformed states.
phi = W(psi)

def Lprime(f):
    return (
        L(f)+I*v*sp.diff(chi,x)*f
    )

assert eq(
    Lprime(phi),
    W(L(psi))
)

print("PASS: transformed state and")
print("      transformed dynamics agree")

# 7. Phase functions depending only
# on p commute with free transport.
chi_p_only = r*p**2

def Wp(f):
    return sp.exp(I*chi_p_only)*f

assert eq(
    L(Wp(psi)),
    Wp(L(psi))
)

print("PASS: momentum-only phases commute")
print("      with free transport")

# 8. Position-dependent phases generally
# do not commute with fixed L.
difference = sp.simplify(
    L(W(psi))-W(L(psi))
)

assert eq(
    difference,
    -I*v*sp.diff(chi,x)*W(psi)
)

assert not eq(
    difference.subs({
        x:0,p:1,c:0,k:1,r:0,
        alpha:1,beta:1,m:1
    }),
    0
)

print("PASS: general local phases are not")
print("      symmetries of fixed L")

# 9. Normalization of the Gaussian.
norm = sp.integrate(
    rho(psi),(x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

assert eq(norm,1)

print("PASS: Gaussian classical density")
print("      is exactly normalized")

# 10. Classical observable expectations
# depend only on the density.
#
# Evaluate representative moments.
density = sp.simplify(rho(psi))

Ex = sp.integrate(
    x*density,(x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

Ep = sp.integrate(
    p*density,(x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

Ex2 = sp.integrate(
    x**2*density,(x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

Ep2 = sp.integrate(
    p**2*density,(x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

assert eq(Ex,0)
assert eq(Ep,0)
assert eq(Ex2,1/(2*alpha))
assert eq(Ep2,1/(2*beta))

# The equality of phase-transformed
# expectations follows pointwise from
# rho(W psi)=rho(psi), already proved.

print("PASS: exact classical moments")
print("      depend only on density")

# 11. Differential operators can distinguish
# phases when treated as observables.
#
# Choose the simple phase chi=k*x.
# Then Dx acting on e^(ikx)*psi has
# expectation k instead of zero.

phi_k = sp.exp(I*k*x)*psi

Dx_expectation_density = sp.simplify(
    sp.conjugate(phi_k)*Dx(phi_k)
)

base_Dx_density = sp.simplify(
    sp.conjugate(psi)*Dx(psi)
)

assert eq(
    Dx_expectation_density-base_Dx_density,
    k*density
)

# The base expectation is zero by parity:
base_expectation = sp.integrate(
    base_Dx_density,
    (x,-sp.oo,sp.oo),
    (p,-sp.oo,sp.oo)
)

assert eq(base_expectation,0)

# The shifted expectation equals k.
assert eq(
    base_expectation+k*norm,
    k
)

print("PASS: derivative operator Dx")
print("      distinguishes phase choices")

# 12. But classical multiplication
# observables cannot distinguish them.
assert eq(rho(phi_k),rho(psi))

assert eq(
    x*rho(phi_k),
    x*rho(psi)
)

assert eq(
    p*rho(phi_k),
    p*rho(psi)
)

print("PASS: classical observables remain")
print("      insensitive to that phase")

print("\n"+"="*64)
print("CHECKPOINT 085: ALL TESTS PASSED")
print("="*64)

print("ESTABLISHED MATHEMATICALLY:")
print("- Real phase multiplication is unitary.")
print("- Classical phase-space densities are")
print("  invariant under such transformations.")
print("- Classical multiplication observables")
print("  remain unchanged.")
print("- Differential operators acquire")
print("  phase-gradient terms.")
print("- Consistently transformed dynamics")
print("  preserves the representation structure.")
print("- Momentum-only phases commute with")
print("  this free classical transport.")
print("- General position-dependent phases")
print("  need not commute with fixed transport.")
print("- Derivative operators can distinguish")
print("  phases that classical observables cannot.")

print("\nMATHEMATICAL QUALIFICATIONS:")
print("- The local-phase identities are")
print("  formal operator identities on suitable")
print("  smooth domains.")
print("- The explicit symbolic tests use")
print("  polynomial real phase functions.")
print("- The Gaussian moment identities are exact.")
print("- Physical observability of derivative")
print("  operators has not been assumed.")

print("\nNOT ESTABLISHED:")
print("- Local phase freedom as a physical gauge law.")
print("- A quantum observable algebra.")
print("- A preferred action scale.")
print("- The quantum Born rule.")
print("- Quantum mechanics emerging from GWFM.")
