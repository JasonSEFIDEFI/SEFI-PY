import sympy as S
from time import monotonic

start = monotonic()

print("=== CP098-T90C-R26P ===", flush=True)

# Minkowski metric: (-,+,+,+)
eta = S.diag(-1, 1, 1, 1)

u = S.Matrix(S.symbols("u0:4"))
q = S.Matrix(S.symbols("q0:4"))

def dot(a, b):
    return (a.T * eta * b)[0]

def directional(expr, variables, vector):
    return sum(
        vector[i] * S.diff(expr, variables[i])
        for i in range(4)
    )

def exact_zero(expr, label):
    residual = S.cancel(S.expand(expr))
    if residual != 0:
        raise AssertionError(
            f"{label}: residual = {residual}"
        )
    print(label, "PASS", flush=True)

# Geometric invariants
s = -dot(u, u)
d = dot(u, q)
Q = dot(q, q)
W = Q + d**2/s

print()
print("STAGE 1: exact geometric identities", flush=True)

identities = {
    "u.du(s)=2s":
        directional(s, u, u) - 2*s,

    "q.du(s)=-2d":
        directional(s, u, q) + 2*d,

    "u.du(d)=d":
        directional(d, u, u) - d,

    "q.du(d)=Q":
        directional(d, u, q) - Q,

    "u.dq(d)=-s":
        directional(d, q, u) + s,

    "q.dq(d)=d":
        directional(d, q, q) - d,

    "u.du(W)=0":
        directional(W, u, u),

    "q.du(W)=2dW/s":
        directional(W, u, q) - 2*d*W/s,

    "u.dq(W)=0":
        directional(W, q, u),

    "q.dq(W)=2W":
        directional(W, q, q) - 2*W
}

for name, residual in identities.items():
    exact_zero(residual, name)

print("STAGE 1 PASS", flush=True)

print()
print("STAGE 2: exact energy derivative identities",
      flush=True)

ss, dd, ww, aa = S.symbols(
    "ss dd ww aa", real=True
)

# On s>0, v=sqrt(s).
# E=2*A*W/v^3.
# A depends on identity variables only.

E = 2*aa*ww/ss**S.Rational(3, 2)
b = E

def energy_derivative(ds, dw):
    return (
        S.diff(E, ss)*ds
        + S.diff(E, ww)*dw
    )

checks = {
    "u.du(E)=-3E":
        energy_derivative(2*ss, 0) + 3*E,

    "q.dq(E)=2E":
        energy_derivative(0, 2*ww) - 2*E,

    "u.dq(E)=0":
        energy_derivative(0, 0),

    "q.du(E)=10AdW/v5":
        energy_derivative(
            -2*dd,
            2*dd*ww/ss
        ) - 10*aa*dd*ww/ss**S.Rational(5, 2)
}

for name, residual in checks.items():
    exact_zero(residual, name)

print("STAGE 2 PASS", flush=True)

print()
print("STAGE 3: canonical Poisson bracket derivation",
      flush=True)

# Canonical pairs:
# (u,p1), (q,P).
#
# phi1 = u.P
# phi2 = q.P
# chi  = u.p1
# psi  = q.p1 + E
#
# Here J represents q.p1.
# phi1, phi2, chi, J are formal independent
# symbols used after applying the product rule.

phi1, phi2, chi, J = S.symbols(
    "phi1 phi2 chi J"
)

psi = J + E

# Product-rule results:
#
# {phi1,psi} = phi2 - chi - u.dq(E)
# {phi2,psi} = -J - q.dq(E)
# {chi,psi}  = J - u.du(E)

bracket_1 = phi2 - chi - energy_derivative(0, 0)

bracket_2 = -J - energy_derivative(0, 2*ww)

bracket_3 = J - energy_derivative(2*ss, 0)

exact_zero(
    bracket_1 - (phi2-chi),
    "{phi1,psi}=phi2-chi"
)

exact_zero(
    bracket_2 - (-psi-b),
    "{phi2,psi}=-psi-b"
)

exact_zero(
    bracket_3 - (psi+2*b),
    "{chi,psi}=psi+2b"
)

print("STAGE 3 PASS", flush=True)

print()
print("STAGE 4: regular-branch algebra",
      flush=True)

# b=2*A*W/v^3.
# On W>0, v>0, b=0 iff A=0.

A0, W0, v0 = S.symbols(
    "A0 W0 v0", real=True, nonzero=True
)

b0 = 2*A0*W0/v0**3

assert S.simplify(b0/A0 - 2*W0/v0**3) == 0

print("b/A=2W/v^3: PASS")
print("A!=0, W>0, v>0 IMPLIES b!=0: PASS")
print("A=0 BRANCH: NOT ANALYZED")

print("STAGE 4 PASS", flush=True)

print()
print("R26P COMPLETE")
print("GENERAL GEOMETRIC IDENTITIES: PASS")
print("ENERGY DERIVATIVE IDENTITIES: PASS")
print("THREE CANONICAL BRACKETS: DERIVED")
print("NO NUMERICAL WITNESSES REQUIRED")
print("REMAINING BRACKETS: OPEN")
print("GENERAL FULL-MATRIX PROOF: OPEN")
print("GAUGE GENERATORS: OPEN")
print("A-ZERO BRANCH: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start, 3)
)
