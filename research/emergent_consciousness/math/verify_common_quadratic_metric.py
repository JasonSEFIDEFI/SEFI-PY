import sympy as sp
import numpy as np

print("Checkpoint 024: Common quadratic metric")

S = sp.Matrix([[1,4],[0,1]])

D1 = sp.diag(
    sp.Rational(4,5),
    sp.Rational(3,10)
)
D2 = sp.diag(
    sp.Rational(2,5),
    sp.Rational(7,10)
)

B = S*D1*S.inv()
C = S*D2*S.inv()

P = S.inv().T*S.inv()

assert B == sp.Matrix([
    [sp.Rational(4,5),-2],
    [0,sp.Rational(3,10)]
])

assert C == sp.Matrix([
    [sp.Rational(2,5),sp.Rational(6,5)],
    [0,sp.Rational(7,10)]
])

assert P == sp.Matrix([[1,-4],[-4,17]])
assert P[0,0] > 0
assert P.det() > 0

print("PASS: exact matrices and positive metric")

q = sp.Rational(16,25)

for M in (B,C):
    difference = q*P-M.T*P*M

    # Positive semidefiniteness, exactly.
    assert difference == difference.T
    assert difference[0,0] >= 0
    assert difference.det() >= 0

    # Strict contraction relative to P.
    strict = P-M.T*P*M
    assert strict[0,0] > 0
    assert strict.det() > 0

print("PASS: exact common contraction certificate")

Bn = np.array(B, dtype=float)
Cn = np.array(C, dtype=float)
Pn = np.array(P, dtype=float)

norm_B = np.linalg.norm(Bn,2)
norm_C = np.linalg.norm(Cn,2)

assert norm_B > 1
assert norm_C > 1

print(f"euclidean_norm_B={norm_B:.9f}")
print(f"euclidean_norm_C={norm_C:.9f}")
print("PASS: Euclidean non-contraction")

# Arbitrary nonperiodic switching example.
schedule = [
    Bn,Cn,Cn,Bn,Bn,Cn,
    Bn,Cn,Bn,Cn,Cn,Bn
]

x = np.array([1.0,1.0])

def energy(x):
    return float(x @ Pn @ x)

initial = energy(x)

for n,M in enumerate(schedule,1):
    previous = energy(x)
    x = M @ x
    current = energy(x)

    assert current <= (
        float(q)*previous + 1e-12
    )

    assert current <= (
        float(q)**n*initial + 1e-12
    )

    print(
        f"step={n:2d} "
        f"P_energy={current:.9e} "
        f"bound={float(q)**n*initial:.9e}"
    )

print("PASS: nonperiodic switching verification")
print("PASS: Checkpoint 024 verified")
