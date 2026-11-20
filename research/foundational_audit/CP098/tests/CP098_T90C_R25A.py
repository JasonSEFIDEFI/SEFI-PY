import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R25A ===", flush=True)

# Two generic identity coordinates are sufficient to
# test the structure; the formulas extend to any finite
# number of independent identity coordinates.

u = S.symbols("u0:4", real=True)
I = S.symbols("I0:2", real=True)
Id = S.symbols("Id0:2", real=True)
h = S.symbols("h0:2", positive=True)
k = S.symbols("k0:2", real=True)
m = S.symbols("m", real=True)

eta = S.diag(-1, 1, 1, 1)
U = S.Matrix(u)
v2 = -(U.T * eta * U)[0]
v = S.sqrt(v2)

K = sum(h[a]*Id[a]**2/2 for a in range(2))
V = sum(k[a]*I[a]**2/2 for a in range(2))

Lextra = -m*v + K/v - v*V

print("STAGE 1: zeroth-momentum contribution", flush=True)

computed = S.Matrix([
    S.diff(Lextra, x) for x in u
])

expected = (
    (m + V + K/v2) * eta*U/v
)

residual = (computed-expected).applyfunc(S.simplify)
print("DELTA p0 RESIDUAL =", residual, flush=True)
assert residual == S.zeros(4,1)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: identity momenta", flush=True)

pI = [S.diff(Lextra, x) for x in Id]
expected_pI = [h[a]*Id[a]/v for a in range(2)]

for a in range(2):
    assert S.simplify(pI[a]-expected_pI[a]) == 0

print("IDENTITY MOMENTUM RESIDUALS = [0, 0]",
      flush=True)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: canonical energy", flush=True)

# Euler contraction for the extra sector.
Hextra = (
    (computed.T*U)[0]
    + sum(pI[a]*Id[a] for a in range(2))
    - Lextra
)

Hresidual = S.simplify(Hextra)
print("EXTRA-SECTOR CANONICAL ENERGY =", Hresidual)
assert Hresidual == 0

print("STAGE 3 PASS", flush=True)

print("STAGE 4: geometric coupling dependence",
      flush=True)

A, C = S.symbols("A C", real=True)
W, Z, ss = S.symbols(
    "W Z ss", positive=True
)

Lgeom = -A*W/ss**S.Rational(3,2) \
        -C*Z/(S.sqrt(ss)*W)

dLdA = S.diff(Lgeom,A)
dLdC = S.diff(Lgeom,C)

assert S.simplify(
    dLdA + W/ss**S.Rational(3,2)
) == 0
assert S.simplify(
    dLdC + Z/(S.sqrt(ss)*W)
) == 0

print("dL/dA =", dLdA)
print("dL/dC =", dLdC)
print("STAGE 4 PASS", flush=True)

print()
print("R25A COMPLETE")
print("EXTRA p0: SYMBOLICALLY VERIFIED")
print("EXTRA p1 AND P: ZERO BY DERIVATIVE ORDER")
print("IDENTITY MOMENTA: SYMBOLICALLY VERIFIED")
print("EXTRA CANONICAL ENERGY: IDENTICALLY ZERO")
print("GEOMETRIC COUPLING DERIVATIVES: VERIFIED")
print("COUPLING JET CONSISTENCY: OPEN")
print("FULL HAMILTONIAN RECONSTRUCTION: OPEN")
print("FULL DIRAC SYSTEM: OPEN")
print("DEFI MANUSCRIPT CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =", round(monotonic()-start,3))
