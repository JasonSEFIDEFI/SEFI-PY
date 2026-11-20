import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26L ===", flush=True)

b = S.symbols("b", nonzero=True)
B = S.symbols("B", real=True)

names = ["phi1", "phi2", "chi", "psi", "theta"]

print("STAGE 1: construct weak matrix", flush=True)

C = S.Matrix([
    [0,  0,   0,   0,    0],
    [0,  0,   0,  -b, -2*b],
    [0,  0,   0, 2*b,  4*b],
    [0,  b, -2*b,  0,   -B],
    [0,2*b, -4*b,  B,    0]
])

assert C + C.T == S.zeros(5)

print("CONSTRAINT ORDER:", names)
S.pprint(C)
print("ANTISYMMETRY: PASS")
print("STAGE 1 PASS", flush=True)

print("STAGE 2: symbolic rank", flush=True)

rank = C.rank()
nullity = 5-rank

print("SYMBOLIC RANK =", rank)
print("SYMBOLIC NULLITY =", nullity)

assert rank == 2
assert nullity == 3

# Explicit nonzero 2x2 minor:
minor = C.extract([1,3],[1,3]).det()
assert S.expand(minor-b**2) == 0

print("NONZERO MINOR =", minor)
print("STAGE 2 PASS", flush=True)

print("STAGE 3: exact null vectors", flush=True)

# Each vector supplies coefficients of a
# weakly first-class constraint combination.
#
# The first is phi1.
#
# The second is 2*phi2+chi.
#
# The third is a B-dependent combination.

vectors = [
    S.Matrix([1,0,0,0,0]),
    S.Matrix([0,2,1,0,0]),
    S.Matrix([0,-B/b,0,2, -1])
]

# Verify the first two explicitly.
for i,vec in enumerate(vectors[:2],1):
    assert C*vec == S.zeros(5)
    print("NULL VECTOR", i, "PASS")

# Obtain the third independently from SymPy.
kernel = C.nullspace()

assert len(kernel) == 3

for vec in kernel:
    assert C*vec == S.zeros(5)

print("SYMBOLIC KERNEL DIMENSION =",
      len(kernel))
print("STAGE 3 PASS", flush=True)

print("STAGE 4: numerical witness", flush=True)

witness = {
    b:S.Rational(387,512),
    B:S.Rational(1204251,81920)
}

Cw = C.subs(witness)

assert Cw.rank() == 2
assert Cw + Cw.T == S.zeros(5)

print("WITNESS b =", witness[b])
print("WITNESS B =", witness[B])
print("WITNESS RANK =", Cw.rank())
print("STAGE 4 PASS", flush=True)

print("STAGE 5: conditional dimension count",
      flush=True)

# Canonical variables:
# y,p0: 8
# u,p1: 8
# q,P: 8
# I,pI: 6
# Total phase-space dimension = 30.

phase_dim = 30
constraints = 5

first_class = nullity
second_class = rank

assert first_class + second_class == constraints
assert second_class % 2 == 0

physical_phase_dim = (
    phase_dim
    -2*first_class
    -second_class
)

assert physical_phase_dim == 22

configuration_dof = physical_phase_dim//2
assert configuration_dof == 11

print("PHASE-SPACE DIMENSION =", phase_dim)
print("FIRST-CLASS COMBINATIONS =", first_class)
print("SECOND-CLASS CONSTRAINTS =", second_class)
print("REDUCED PHASE-SPACE DIMENSION =",
      physical_phase_dim)
print("CANONICAL CONFIGURATION DOF =",
      configuration_dof)

print("STAGE 5 PASS", flush=True)

print()
print("R26L COMPLETE")
print("WEAK CONSTRAINT MATRIX: CONSTRUCTED")
print("SYMBOLIC RANK: VERIFIED")
print("EXACT WITNESS RANK: VERIFIED")
print("FIRST/SECOND CLASS COUNT: CONDITIONAL")
print("A-NONZERO REGULAR BRANCH ONLY")
print("A-ZERO BRANCH: OPEN")
print("GAUGE GENERATORS: NOT CONSTRUCTED")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-start,3))
