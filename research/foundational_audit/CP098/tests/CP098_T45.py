"""CP098-T45: Off-surface Hamiltonian identity.

Checks Hc = v*psi_T + v^2*k*chi_N + a*chi_T
in the regular adapted frame, with all canonical
momentum components unrestricted.

Also checks a second, non-adapted Lorentz frame.

This does not yet calculate off-surface Poisson
brackets or establish a gauge generator.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)

A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

P0, P1, P2, P3 = sp.symbols(
    "P0 P1 P2 P3", real=True
)

p10, p11, p12, p13 = sp.symbols(
    "p10 p11 p12 p13", real=True
)

p00, p01, p02, p03 = sp.symbols(
    "p00 p01 p02 p03", real=True
)

P = sp.Matrix([P0, P1, P2, P3])
p1 = sp.Matrix([p10, p11, p12, p13])
p0 = sp.Matrix([p00, p01, p02, p03])

def dot(X, Y):
    return (X.T * eta * Y)[0]

def pair(p, X):
    return (p.T * X)[0]

def quantities(u, q):
    speed = sp.sqrt(-dot(u, u))
    T = u/speed
    tq = dot(T, q)
    w = q + T*tq
    w2 = dot(w, w)
    kap = sp.sqrt(w2)/speed**2
    N = w/sp.sqrt(w2)

    phi_T = pair(P, T)
    phi_N = pair(P, N)

    chi_T = pair(p1, T)
    chi_N = pair(p1, N) + 2*A*kap/speed

    psi_T = (
        pair(p0, T)
        + B
        - A*kap**2
        - w2*dot(P, P)/(4*C)
    )

    Hc = (
        pair(p0, u)
        + pair(p1, q)
        + speed*(B + A*kap**2)
        - speed*w2*dot(P, P)/(4*C)
    )

    reconstructed = (
        speed*psi_T
        + speed**2*kap*chi_N
        - tq*chi_T
    )

    return {
        "speed": speed,
        "kap": kap,
        "tq": tq,
        "phi_T": phi_T,
        "phi_N": phi_N,
        "chi_T": chi_T,
        "chi_N": chi_N,
        "psi_T": psi_T,
        "Hc": Hc,
        "reconstructed": reconstructed
    }

print("=== CP098-T45 ===", flush=True)

# Test 1: adapted frame, unrestricted momenta.
u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

F = quantities(u, q)

print("Checking adapted frame...", flush=True)

assert sp.simplify(F["speed"] - v) == 0
assert sp.simplify(F["kap"] - k) == 0
assert sp.simplify(F["tq"] + a) == 0

residual = sp.factor(
    F["Hc"] - F["reconstructed"]
)

print("Adapted identity residual =", residual,
      flush=True)
assert residual == 0

print("PASS: adapted off-surface identity",
      flush=True)

# Verify the constraints have NOT been imposed.
# In particular, phi_T=P0 and phi_N=P1.
assert sp.simplify(F["phi_T"]-P0) == 0
assert sp.simplify(F["phi_N"]-P1) == 0

print("PASS: primary momenta unrestricted",
      flush=True)

# Verify independence from every momentum
# component by checking the identity residual
# before imposing any momentum constraints.
for z in list(P)+list(p1)+list(p0):
    assert sp.diff(residual, z) == 0

print("PASS: identity holds for arbitrary momenta",
      flush=True)

# Test 2: non-adapted timelike u and
# non-collinear q, with exact rational data.
#
# u^2 = -16, speed=4.
# q has a nonzero transverse component.
u2 = sp.Matrix([5, 3, 0, 0])
q2 = sp.Matrix([2, 1, 4, 0])

assert -dot(u2, u2) == 16

G = quantities(u2, q2)

print("Checking non-adapted frame...", flush=True)

assert sp.simplify(G["kap"]) > 0

residual2 = sp.simplify(
    G["Hc"] - G["reconstructed"]
)

print("Non-adapted identity residual =",
      residual2, flush=True)

assert residual2 == 0

print("PASS: non-adapted off-surface identity",
      flush=True)

print("OFF-SURFACE HAMILTONIAN IDENTITY: VERIFIED")
print("OFF-SURFACE POISSON CHAINS: OPEN")
print("GAUGE GENERATOR: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
