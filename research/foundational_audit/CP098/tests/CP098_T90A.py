"""CP098-T90A: total-Hamiltonian tangency audit.

Five constraints, two primary multipliers.
Exact canonical brackets at the T89A regular point.

No gauge or first-class classification assumed.
"""
import runpy
from pathlib import Path
import sympy as sp

print("=== CP098-T90A ===", flush=True)

root = Path("research/foundational_audit/CP098/tests")
baseline = root / "CP098_T89A.py"

assert baseline.is_file(), "STOP: T89A missing"

data = runpy.run_path(str(baseline))

eta = data["eta"]
constraints = data["constraints"]
J = data["J"]
z0 = data["z0"]
M0 = data["M"]

eps = sp.symbols("eps", real=True)

lamT, lamN = sp.symbols(
    "lambda_T lambda_N", real=True
)

def dot(x,y):
    return (x.T*eta*y)[0]

def pair(p,x):
    return (p.T*x)[0]

def canonical_H(z):
    u = z[0:4,0]
    q = z[4:8,0]
    P = z[8:12,0]
    j = z[12:16,0]
    p0 = z[16:20,0]

    I = z[20]
    pI = z[21]

    A = (
        sp.Rational(3,2)
        + I + I**2/3
    )
    C = (
        sp.Rational(5,3)
        + I/2 + I**2/5
    )
    B = (
        sp.Rational(7,4)
        + sp.Rational(2,5)*pI
        + I*pI/3
        + pI**2/7
        + I**2/11
    )

    v = sp.sqrt(-dot(u,u))
    T = u/v

    w = q+T*dot(T,q)
    L = sp.sqrt(dot(w,w))
    k = L/v**2

    return (
        pair(p0,u)
        +pair(j,q)
        +v*(B+A*k**2)
        -v**5*k**2*dot(P,P)/(4*C)
    )

def gradient(fn):
    entries = []

    for i in range(22):
        direction = sp.eye(22)[:,i]

        expr = fn(z0+eps*direction)

        derivative = sp.diff(
            expr,eps
        ).subs(eps,0)

        entries.append(sp.simplify(derivative))

    return sp.Matrix(entries)

print("Computing constraint gradients...",
      flush=True)

G = []

for i in range(5):
    G.append(
        gradient(lambda z, i=i: constraints(z)[i])
    )
    print(f"CONSTRAINT_GRADIENT_{i+1}_COMPLETE",
          flush=True)

print("Computing canonical Hamiltonian gradient...",
      flush=True)

GH = gradient(canonical_H)

# Independent reconstruction of the
# five-by-five bracket matrix.
M = sp.zeros(5)

for i in range(5):
    for j in range(i+1,5):
        value = sp.simplify(
            (G[i].T*J*G[j])[0]
        )

        M[i,j] = value
        M[j,i] = -value

matrix_residual = M-M0

print("BRACKET_MATRIX_RESIDUAL =",
      matrix_residual,flush=True)

assert matrix_residual == sp.zeros(5)

# Hamiltonian tangency source vector:
# h_i = {Phi_i,H_c}.
h = sp.Matrix([
    sp.simplify((g.T*J*GH)[0])
    for g in G
])

print("HAMILTONIAN_BRACKET_VECTOR =",
      list(h),flush=True)

# Total-Hamiltonian tangency:
# {Phi_i,H_T} = h_i
#   +lambda_T M[i,0]
#   +lambda_N M[i,1].
tangency = sp.Matrix([
    sp.factor(
        h[i]
        +lamT*M[i,0]
        +lamN*M[i,1]
    )
    for i in range(5)
])

print("TANGENCY_EQUATIONS:",flush=True)

for i,expr in enumerate(tangency,1):
    print(f"DOT_PHI_{i} =",expr,flush=True)

# Analyze the multiplier system without
# assuming consistency or uniqueness.
coefficient = tangency.jacobian(
    [lamT,lamN]
)

constant = tangency.subs({
    lamT:0,
    lamN:0
})

augmented = coefficient.row_join(-constant)

rank_coefficient = coefficient.rank()
rank_augmented = augmented.rank()

print("MULTIPLIER_COEFFICIENT_RANK =",
      rank_coefficient,flush=True)
print("MULTIPLIER_AUGMENTED_RANK =",
      rank_augmented,flush=True)

solutions = sp.linsolve(
    (coefficient,-constant),
    (lamT,lamN)
)

print("MULTIPLIER_SOLUTION_SET =",
      solutions,flush=True)

if rank_augmented > rank_coefficient:
    print(
        "STOP-GATE: tangency equations incompatible "
        "at this phase-space point",
        flush=True
    )
    print(
        "POSSIBLE ADDITIONAL CONSTRAINTS: "
        "REQUIRE INDEPENDENT ANALYSIS",
        flush=True
    )
elif rank_coefficient < 2:
    print(
        "NOTICE: at least one multiplier remains "
        "undetermined at this point",
        flush=True
    )
else:
    print(
        "NOTICE: multipliers determined at this point",
        flush=True
    )

# Cross-check direct total-Hamiltonian gradient
# at one explicit multiplier assignment.
def total_H(z):
    f = constraints(z)
    return (
        canonical_H(z)
        +lamT*f[0]
        +lamN*f[1]
    )

GHT = gradient(total_H)

direct = sp.Matrix([
    sp.simplify((g.T*J*GHT)[0])
    for g in G
])

crosscheck = (direct-tangency).applyfunc(
    sp.simplify
)

print("TOTAL_HAMILTONIAN_CROSSCHECK =",
      list(crosscheck),flush=True)

assert crosscheck == sp.zeros(5,1)

print("PASS: exact total-Hamiltonian tangency audit")
print("POINTWISE CONSISTENCY EQUATIONS: COMPUTED")
print("TERTIARY-CONSTRAINT ANALYSIS: OPEN")
print("GLOBAL DIRAC CONSISTENCY: OPEN")
print("FIRST-CLASS CLASSIFICATION: OPEN")
print("GAUGE INTERPRETATION: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T90A COMPLETE")
