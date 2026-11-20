import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26C ===", flush=True)

# Canonical Ostrogradsky coordinates.
y = S.symbols("y0:4", real=True)
u = S.symbols("u0:4", real=True)
q = S.symbols("q0:4", real=True)

p0 = S.symbols("p0_0:4", real=True)
p1 = S.symbols("p1_0:4", real=True)
P = S.symbols("P0:4", real=True)

# Regular identity canonical coordinates.
I = S.symbols("IA IS IW", real=True)
pI = S.symbols("pIA pIS pIW", real=True)

coordinates = y + u + q + I
momenta = p0 + p1 + P + pI

assert len(coordinates) == len(momenta) == 15

def PB(f, g):
    result = S.Integer(0)
    for x, p in zip(coordinates, momenta):
        result += (
            S.diff(f, x)*S.diff(g, p)
            - S.diff(f, p)*S.diff(g, x)
        )
    return S.expand(result)

phi1 = sum(u[i]*P[i] for i in range(4))
phi2 = sum(q[i]*P[i] for i in range(4))

print("STAGE 1: canonical bracket definition",
      flush=True)

for i in range(4):
    assert PB(q[i], P[i]) == 1
    assert PB(u[i], p1[i]) == 1
    assert PB(y[i], p0[i]) == 1

for i in range(3):
    assert PB(I[i], pI[i]) == 1

assert PB(q[0], p1[0]) == 0
assert PB(u[0], P[0]) == 0

print("CANONICAL PAIRS: PASS")

print("STAGE 2: primary-constraint brackets",
      flush=True)

b11 = PB(phi1, phi1)
b12 = PB(phi1, phi2)
b21 = PB(phi2, phi1)
b22 = PB(phi2, phi2)

bracket_matrix = S.Matrix([
    [b11, b12],
    [b21, b22]
])

print("BRACKET MATRIX =")
S.pprint(bracket_matrix)

assert b11 == 0
assert b22 == 0
assert S.expand(b12 + phi1) == 0
assert S.expand(b21 - phi1) == 0

print("{phi1,phi2} = -phi1")
print("{phi2,phi1} = phi1")
print("STAGE 2 PASS")

print("STAGE 3: constraint-surface closure",
      flush=True)

weak_matrix = bracket_matrix.subs(phi1, 0)

assert weak_matrix == S.zeros(2)

print("BRACKET MATRIX ON phi1=phi2=0 =")
S.pprint(weak_matrix)

print("PRIMARY-SUBALGEBRA CLOSURE: PASS")

print("STAGE 4: generated canonical flows",
      flush=True)

# Infinitesimal transformation:
# delta F = epsilon * {F, phi}.

for i in range(4):
    # phi1 generates a shift of q along u.
    assert PB(q[i], phi1) == u[i]
    assert PB(P[i], phi1) == 0
    assert PB(p1[i], phi1) == -P[i]

    # phi2 generates reciprocal scaling of q,P.
    assert PB(q[i], phi2) == q[i]
    assert PB(P[i], phi2) == -P[i]
    assert PB(p1[i], phi2) == 0

    # Neither generator changes y or p0.
    for phi in (phi1, phi2):
        assert PB(y[i], phi) == 0
        assert PB(p0[i], phi) == 0

for i in range(3):
    for phi in (phi1, phi2):
        assert PB(I[i], phi) == 0
        assert PB(pI[i], phi) == 0

print("phi1: delta q=u, delta P=0, delta p1=-P")
print("phi2: delta q=q, delta P=-P")
print("IDENTITY CANONICAL PAIRS: UNCHANGED")
print("STAGE 4 PASS")

print()
print("R26C COMPLETE")
print("CANONICAL POISSON BRACKETS: VERIFIED")
print("PRIMARY CONSTRAINT SUBALGEBRA: CLOSED")
print("PRIMARY BRACKET MATRIX: WEAKLY ZERO")
print("GENERATED CANONICAL FLOWS: VERIFIED")
print("FULL FIRST-CLASS CLASSIFICATION: OPEN")
print("SECONDARY CONSTRAINTS: NOT DERIVED")
print("DIRAC CONSISTENCY: OPEN")
print("PHYSICAL DEGREES OF FREEDOM: NOT COUNTED")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =",
      round(monotonic()-start, 3))
