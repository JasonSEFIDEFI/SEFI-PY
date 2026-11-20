"""CP098-T89B: exact rank and nullspace audit.

Input: exact on-surface bracket matrix from T89A.
No general-sector or gauge classification claimed.
"""
import sympy as sp

print("=== CP098-T89B ===", flush=True)

Q = sp.Rational

M = sp.Matrix([
    [0,       0,        0,        0,           0],
    [0,       0,        0, -Q(3,64),     Q(3,32)],
    [0,       0,        0,  Q(3,16),    -Q(3,8)],
    [0, Q(3,64), -Q(3,16),        0, -Q(64333,640)],
    [0,-Q(3,32),  Q(3,8), Q(64333,640),        0],
])

assert M + M.T == sp.zeros(5)

names = ("phiT", "phiN", "chiT", "chiN", "psiT")

print("MATRIX_ANTISYMMETRY_RESIDUAL =",
      M + M.T, flush=True)

rank = M.rank()
nullity = len(M.nullspace())

print("FULL_MATRIX_RANK =", rank, flush=True)
print("FULL_MATRIX_NULLITY =", nullity, flush=True)

assert rank + nullity == 5
assert rank % 2 == 0

# The submatrix for (phiN, chiT, chiN, psiT).
K = M[1:5, 1:5]

# Pfaffian for a 4x4 antisymmetric matrix.
pf = sp.factor(
    K[0,1]*K[2,3]
    - K[0,2]*K[1,3]
    + K[0,3]*K[1,2]
)

print("FOUR_BLOCK_PFAFFIAN =", pf, flush=True)
print("FOUR_BLOCK_DETERMINANT =", K.det(),
      flush=True)

assert sp.simplify(K.det() - pf**2) == 0

# Verify a nonzero 2x2 minor.
minor = sp.simplify(
    M.extract([1,3],[1,3]).det()
)

print("NONZERO_2X2_MINOR =", minor, flush=True)
assert minor != 0

print("EXACT_NULLSPACE_BASIS:", flush=True)

vectors = M.nullspace()

for index, vector in enumerate(vectors, 1):
    residual = M*vector

    print(f"NULL_VECTOR_{index} =", list(vector),
          flush=True)
    print(f"NULL_RESIDUAL_{index} =",
          list(residual), flush=True)

    assert residual == sp.zeros(5,1)

# Confirm basis vectors are independent.
V = sp.Matrix.hstack(*vectors)

assert V.rank() == nullity
assert M*V == sp.zeros(5,nullity)

# A single 5x5 antisymmetric matrix cannot
# have rank 5. Check the actual rank rather
# than relying only on that general fact.
if rank == 2:
    print("PASS: rank-two matrix at tested point")
elif rank == 4:
    print("PASS: rank-four matrix at tested point")
else:
    raise AssertionError(
        "STOP: unexpected constraint-matrix rank"
    )

print("PASS: exact Pfaffian and nullspace audit")
print("POINTWISE RANK: VERIFIED")
print("REGULAR-SECTOR CONSTANT RANK: OPEN")
print("FIRST-CLASS CLASSIFICATION: OPEN")
print("GAUGE GENERATORS: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T89B COMPLETE")
