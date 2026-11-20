"""CP098-T89G: exact non-adapted Lorentz audit.

Boost the T89A regular constraint-surface point.
Recompute canonical derivatives at the boosted point.
Check covariance, matrix equality, rank, and nullity.
"""
import runpy
from pathlib import Path
import sympy as sp

print("=== CP098-T89G ===", flush=True)

root = Path("research/foundational_audit/CP098/tests")
baseline = root / "CP098_T89A.py"

assert baseline.is_file(), "STOP: T89A missing"

# Load the previously verified constraint definitions.
# This also reruns T89A's original assertions.
data = runpy.run_path(str(baseline))

eta = data["eta"]
F = data["constraints"]
J = data["J"]
z0 = data["z0"]
M0 = data["M"]

e = sp.symbols("boost_eps", real=True)

Q = sp.Rational

# Exact Lorentz boost in the 0-2 plane.
Lambda = sp.Matrix([
    [Q(5,4), 0, Q(3,4), 0],
    [0, 1, 0, 0],
    [Q(3,4), 0, Q(5,4), 0],
    [0, 0, 0, 1],
])

metric_residual = sp.simplify(
    Lambda.T*eta*Lambda-eta
)

print("LORENTZ_METRIC_RESIDUAL =",
      metric_residual, flush=True)

assert metric_residual == sp.zeros(4)

# Contravariant vectors transform by Lambda.
# Covectors transform by inverse transpose.
cov_transform = Lambda.inv().T

assert sp.simplify(
    cov_transform.T*Lambda-sp.eye(4)
) == sp.zeros(4)

def transform_vector(z, start):
    return Lambda*z[start:start+4,0]

def transform_covector(z, start):
    return cov_transform*z[start:start+4,0]

u1 = transform_vector(z0,0)
q1 = transform_vector(z0,4)
P1 = transform_covector(z0,8)
j1 = transform_covector(z0,12)
p01 = transform_covector(z0,16)

z1 = sp.Matrix(
    list(u1)+list(q1)+list(P1)
    +list(j1)+list(p01)
    +[z0[20],z0[21]]
)

assert u1[2] != 0
assert q1[2] != 0

print("BOOSTED_U =",list(u1),flush=True)
print("BOOSTED_Q =",list(q1),flush=True)

# Verify the transformed point remains on
# the complete five-constraint surface.
for i,value in enumerate(F(z1),1):
    residual = sp.simplify(value)

    print(f"CONSTRAINT_{i}_RESIDUAL =",
          residual,flush=True)

    assert residual == 0, (
        f"STOP: constraint {i} not preserved"
    )

print("Computing boosted canonical gradients...",
      flush=True)

gradients = []

for i in range(5):
    entries = []

    for k in range(22):
        direction = sp.eye(22)[:,k]

        expression = F(z1+e*direction)[i]

        derivative = sp.diff(
            expression,e
        ).subs(e,0)

        entries.append(sp.simplify(derivative))

    gradients.append(sp.Matrix(entries))

    print(f"GRADIENT_{i+1}_COMPLETE",
          flush=True)

M1 = sp.zeros(5)

for i in range(5):
    for j in range(i+1,5):
        entry = sp.simplify(
            (gradients[i].T*J*gradients[j])[0]
        )

        M1[i,j] = entry
        M1[j,i] = -entry

assert M1+M1.T == sp.zeros(5)

print("BOOSTED_BRACKET_MATRIX =",flush=True)
sp.pprint(M1)

difference = sp.simplify(M1-M0)

print("LORENTZ_MATRIX_RESIDUAL =",
      difference,flush=True)

assert difference == sp.zeros(5), (
    "STOP: bracket matrix not preserved"
)

rank = M1.rank()
nullity = len(M1.nullspace())

print("BOOSTED_MATRIX_RANK =",rank,flush=True)
print("BOOSTED_MATRIX_NULLITY =",nullity,flush=True)

assert rank == 2
assert nullity == 3

K = M1[1:5,1:5]

pf = sp.factor(
    K[0,1]*K[2,3]
    -K[0,2]*K[1,3]
    +K[0,3]*K[1,2]
)

print("BOOSTED_PFAFFIAN =",pf,flush=True)

assert pf == 0

for i,n in enumerate(M1.nullspace(),1):
    assert M1*n == sp.zeros(5,1)
    print(f"NULL_VECTOR_{i} =",list(n),flush=True)

print("PASS: exact Lorentz transformation")
print("PASS: all five constraints preserved")
print("PASS: ten boosted canonical brackets")
print("PASS: bracket matrix covariance")
print("PASS: boosted rank two and nullity three")
print("NON-ADAPTED BOOST ORBIT: VERIFIED")
print("GENERAL REGULAR-SECTOR PROOF: OPEN")
print("FIRST-CLASS CLASSIFICATION: OPEN")
print("GAUGE GENERATORS: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T89G COMPLETE")
