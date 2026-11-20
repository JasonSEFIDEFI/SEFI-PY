import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26A ===", flush=True)

v, m = S.symbols("v m", positive=True)

hA, hS, hW = S.symbols(
    "hA hS hW", positive=True
)
kA, kS, kW = S.symbols(
    "kA kS kW", real=True
)

IA, IS, IW = S.symbols("IA IS IW", real=True)
dA, dS, dW = S.symbols("dA dS dW", real=True)

hs = (hA, hS, hW)
ks = (kA, kS, kW)
positions = (IA, IS, IW)
velocities = (dA, dS, dW)

K = sum(
    hs[i]*velocities[i]**2/2
    for i in range(3)
)

V = sum(
    ks[i]*positions[i]**2/2
    for i in range(3)
)

Lextra = -m*v + K/v - v*V

print("STAGE 1: identity Legendre map", flush=True)

momenta = [
    S.diff(Lextra, velocity)
    for velocity in velocities
]

expected = [
    hs[i]*velocities[i]/v
    for i in range(3)
]

assert all(
    S.simplify(momenta[i]-expected[i]) == 0
    for i in range(3)
)

print("IDENTITY MOMENTA =", momenta)

velocity_hessian = S.hessian(
    Lextra, velocities
)

expected_hessian = S.diag(
    hA/v, hS/v, hW/v
)

assert velocity_hessian == expected_hessian

determinant = S.factor(velocity_hessian.det())

print("VELOCITY HESSIAN DETERMINANT =", determinant)
assert determinant == hA*hS*hW/v**3

print("STAGE 1 PASS", flush=True)

print("STAGE 2: partial identity Hamiltonian", flush=True)

pA, pS, pW = S.symbols("pA pS pW", real=True)
ps = (pA, pS, pW)

inverse = {
    velocities[i]: v*ps[i]/hs[i]
    for i in range(3)
}

for i in range(3):
    assert S.simplify(
        momenta[i].subs(inverse)-ps[i]
    ) == 0

Hidentity = sum(
    momenta[i]*velocities[i]
    for i in range(3)
) - (K/v-v*V)

Hidentity = S.factor(Hidentity.subs(inverse))

expected_Hidentity = v*(
    sum(ps[i]**2/(2*hs[i]) for i in range(3))
    + V
)

assert S.simplify(
    Hidentity-expected_Hidentity
) == 0

print("PARTIAL IDENTITY HAMILTONIAN =", Hidentity)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: full p0 contribution", flush=True)

# From R25A:
# delta_p0 = (m + V + K/v**2)*eta*u/v
#
# Since u dot_eta u = -v**2,
# delta_p0 dot u = -v*(m+V+K/v**2).

delta_p0_dot_u = -v*(m+V+K/v**2)

# Extra-sector contribution to the COMPLETE
# Ostrogradsky canonical energy.
#
# The geometric contribution was checked in R24.
Hextra = (
    delta_p0_dot_u
    + sum(
        momenta[i]*velocities[i]
        for i in range(3)
    )
    - Lextra
)

Hextra = S.factor(Hextra)

print("EXTRA CANONICAL ENERGY =", Hextra)
assert Hextra == 0

print("STAGE 3 PASS", flush=True)

print("STAGE 4: independent algebraic identity",
      flush=True)

# Sum pI*I' = 2K/v.
kinetic_identity = S.factor(
    sum(
        momenta[i]*velocities[i]
        for i in range(3)
    ) - 2*K/v
)

assert kinetic_identity == 0

# Independent decomposition of the extra
# canonical energy.
decomposition = S.factor(
    delta_p0_dot_u
    + 2*K/v
    + m*v
    - K/v
    + v*V
)

assert decomposition == 0

print("KINETIC IDENTITY RESIDUAL =", kinetic_identity)
print("FULL CANCELLATION RESIDUAL =", decomposition)

print("STAGE 4 PASS", flush=True)

print()
print("R26A COMPLETE")
print("IDENTITY LEGENDRE MAP: INVERTIBLE")
print("PARTIAL IDENTITY HAMILTONIAN: VERIFIED")
print("EXTRA CANONICAL ENERGY: IDENTICALLY ZERO")
print("FULL GEOMETRIC ENERGY: PREVIOUSLY CHECKED")
print("FULL CANONICAL ENERGY: ZERO CONDITIONAL ON R24")
print("PRIMARY CONSTRAINT CLASSIFICATION: OPEN")
print("FULL DIRAC ALGORITHM: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-start, 3))
