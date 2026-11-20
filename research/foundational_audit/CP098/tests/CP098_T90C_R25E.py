import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R25E ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

# First exact worldline jet.
u = S.Matrix([5, 3, 0, 0])
q = S.Matrix([1, 2, 2, 1])
j = S.Matrix([4, 1, 3, 6])

ss = -dot(u, u)
v = S.sqrt(ss)
vp = -dot(u, q)/v

d = dot(u, q)
Q = dot(q, q)
z = dot(u, j)
h = dot(q, j)
J = dot(j, j)

W = S.factor(Q + d**2/ss)
Z = S.factor(
    J + z**2/ss - (h + d*z/ss)**2/W
)

assert ss > 0 and W > 0

ia, is_, iw = S.symbols("ia is iw", real=True)
coords = (ia, is_, iw)

A = (3 + 2*ia)*is_
C = (2 + ia + S.Rational(3,2)*ia**2)*iw

Lgeom = -A*W/v**3 - C*Z/(v*W)

# R25B initial identity data.
positions = (S.Integer(0), S.Integer(1), S.Integer(1))
velocities = (S.Integer(1), S.Integer(0), S.Integer(0))

# Illustrative, explicitly positive kinetic and
# potential coefficients.
hs = (S.Integer(1),)*3
ks = (S.Integer(1),)*3

point = dict(zip(coords, positions))

print("STAGE 1: geometric identity sources", flush=True)

sources = [
    S.factor(S.diff(Lgeom, x).subs(point))
    for x in coords
]

print("SOURCES =", sources, flush=True)
print("v =", v, "vprime =", vp, flush=True)

print("STAGE 2: solve identity accelerations", flush=True)

acc = [
    S.factor(
        velocities[i]*vp/v
        + v*(sources[i] - v*ks[i]*positions[i])/hs[i]
    )
    for i in range(3)
]

print("IA'' =", acc[0])
print("IS'' =", acc[1])
print("IW'' =", acc[2])

# Verify the original Euler-Lagrange equations.
residuals = [
    S.factor(
        hs[i]*(acc[i]/v - velocities[i]*vp/v**2)
        + v*ks[i]*positions[i]
        - sources[i]
    )
    for i in range(3)
]

print("EL RESIDUALS =", residuals, flush=True)
assert residuals == [0, 0, 0]

print("STAGE 2 PASS", flush=True)

print("STAGE 3: induced coupling jets", flush=True)

lam = S.symbols("lam", real=True)

paths = {
    coords[i]: (
        positions[i]
        + velocities[i]*lam
        + acc[i]*lam**2/2
    )
    for i in range(3)
}

At = A.subs(paths)
Ct = C.subs(paths)

A0 = S.factor(At.subs(lam, 0))
A1 = S.factor(S.diff(At, lam).subs(lam, 0))
A2 = S.factor(S.diff(At, lam, 2).subs(lam, 0))

C0 = S.factor(Ct.subs(lam, 0))
C1 = S.factor(S.diff(Ct, lam).subs(lam, 0))
C2 = S.factor(S.diff(Ct, lam, 2).subs(lam, 0))

print("A, A', A'' =", (A0, A1, A2))
print("C, C', C'' =", (C0, C1, C2))

assert (A0, A1, C0, C1) == (3, 2, 2, 1)

difference = S.factor(C2 - 3)

print("R22-R24 PRESCRIBED C'' =", 3)
print("DYNAMIC C'' MINUS PRESCRIBED =", difference)
print("C'' COMPATIBLE =", difference == 0)

print("STAGE 3 PASS", flush=True)

print()
print("R25E COMPLETE")
print("IDENTITY ACCELERATIONS: SOLVED AT ONE JET")
print("IDENTITY EL RESIDUALS: EXACTLY ZERO")
print("PRESCRIBED COUPLING-JET COMPATIBILITY:",
      "PASS" if difference == 0 else "FAIL")
print("FULL COUPLED TRAJECTORY: OPEN")
print("WORLDLINE EQUATIONS: OPEN")
print("HAMILTONIAN RECONSTRUCTION: OPEN")
print("DIRAC SYSTEM: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =",
      round(monotonic()-start, 3))
