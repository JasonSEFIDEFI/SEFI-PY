"""CP098-T89H: symbolic rank and nullspace certificate.

Uses the bracket structure established in T89F.
D = {chi_N, psi_T} remains unrestricted.

This is an algebraic consequence audit, not an
independent recalculation of canonical brackets.
"""
import sympy as sp

print("=== CP098-T89H ===", flush=True)

v, k = sp.symbols(
    "v k", positive=True, finite=True
)
A = sp.symbols("A", real=True, nonzero=True)
D = sp.symbols("D", real=True)

# Ordering:
# phi_T, phi_N, chi_T, chi_N, psi_T

a = -2*A/v**3
b = 2*A*k/v**2
c = 4*A*k/v**2
d = -4*A*k**2/v

M = sp.Matrix([
    [0, 0, 0,  0,  0],
    [0, 0, 0,  a,  b],
    [0, 0, 0,  c,  d],
    [0,-a,-c,  0,  D],
    [0,-b,-d, -D,  0],
])

assert M + M.T == sp.zeros(5)

# Check the two proportionality relations
# responsible for the rank deficiency.
ratio1 = sp.simplify(c + 2*k*v*a)
ratio2 = sp.simplify(d + 2*k*v*b)

print("ROW_DEPENDENCY_1_RESIDUAL =",ratio1,
      flush=True)
print("ROW_DEPENDENCY_2_RESIDUAL =",ratio2,
      flush=True)

assert ratio1 == 0
assert ratio2 == 0

# Exact Pfaffian of the 4x4 submatrix.
K = M[1:5,1:5]

pf = sp.factor(
    K[0,1]*K[2,3]
    -K[0,2]*K[1,3]
    +K[0,3]*K[1,2]
)

minor = sp.factor(
    M.extract([1,3],[1,3]).det()
)

print("GENERAL_PFAFFIAN =",pf,flush=True)
print("GENERAL_MINOR =",minor,flush=True)

assert pf == 0
assert sp.simplify(
    minor-4*A**2/v**6
) == 0

# Three candidate null vectors, with no
# numerical values assigned to D.
n1 = sp.Matrix([1,0,0,0,0])

n2 = sp.Matrix([
    0,2*k*v,1,0,0
])

n3 = sp.Matrix([
    0,D/a,0,k*v,1
])

V = sp.Matrix.hstack(n1,n2,n3)

for index,n in enumerate((n1,n2,n3),1):
    residual = (M*n).applyfunc(sp.simplify)

    print(
        f"NULL_VECTOR_{index} =",
        list(n),flush=True
    )
    print(
        f"NULL_RESIDUAL_{index} =",
        list(residual),flush=True
    )

    assert residual == sp.zeros(5,1)

assert V.rank() == 3

# Rank upper bound follows from the
# Pfaffian and zero first row/column.
# Rank lower bound follows from minor != 0.
assert sp.simplify(minor).is_zero is False

rank = M.rank()

print("SYMBOLIC_RANK =",rank,flush=True)
print("SYMBOLIC_NULLITY =",5-rank,flush=True)

assert rank == 2
assert 5-rank == 3

# Reproduce both independently computed
# numerical matrices.
def check_case(name, values, expected):
    actual = M.subs(values)
    residual = (actual-expected).applyfunc(
        sp.simplify
    )

    print(
        f"{name}_REPRODUCTION_RESIDUAL =",
        residual,flush=True
    )

    assert residual == sp.zeros(5)

Q = sp.Rational

M89A = sp.Matrix([
    [0,0,0,0,0],
    [0,0,0,-Q(3,64),Q(3,32)],
    [0,0,0,Q(3,16),-Q(3,8)],
    [0,Q(3,64),-Q(3,16),0,-Q(64333,640)],
    [0,-Q(3,32),Q(3,8),Q(64333,640),0],
])

M89C = sp.Matrix([
    [0,0,0,0,0],
    [0,0,0,-Q(1,9),Q(2,3)],
    [0,0,0,Q(4,3),-8],
    [0,Q(1,9),-Q(4,3),0,-Q(24227,45)],
    [0,-Q(2,3),8,Q(24227,45),0],
])

check_case(
    "T89A",
    {v:4,k:Q(1,2),A:Q(3,2),
     D:-Q(64333,640)},
    M89A
)

check_case(
    "T89C",
    {v:3,k:2,A:Q(3,2),
     D:-Q(24227,45)},
    M89C
)

print("PASS: exact symbolic rank-two certificate")
print("PASS: three symbolic null directions")
print("PASS: both numerical matrices reproduced")
print("ALGEBRAIC RANK STRUCTURE: VERIFIED")
print("GLOBAL REGULAR-SECTOR ARGUMENT: TO REVIEW")
print("DIRAC CONSISTENCY: OPEN")
print("FIRST-CLASS CLASSIFICATION: OPEN")
print("GAUGE INTERPRETATION: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T89H COMPLETE")
