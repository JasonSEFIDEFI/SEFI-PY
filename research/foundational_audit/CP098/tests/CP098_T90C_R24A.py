import sympy as s
import time

print("=== CP098-T90C-R24A ===", flush=True)

start = time.monotonic()

def stage(message):
    elapsed = time.monotonic() - start
    print(f"[{elapsed:.1f}s] {message}", flush=True)

eta = s.diag(-1, 1, 1, 1)
t, e = s.symbols("t e", real=True)

def dot(x, y):
    return (x.T * eta * y)[0]

def lagrangian(u, q, j, A, C):
    v2 = -dot(u, u)
    v = s.sqrt(v2)

    # Avoid square roots in the transverse projector.
    uq = dot(u, q)
    w = q + uq*u/v2
    w2 = dot(w, w)

    uj = dot(u, j)
    wj = dot(w, j)

    R = j + uj*u/v2 - wj*w/w2

    return (
        -A*w2/v**3
        -C*dot(R, R)/(v*w2)
    )

def highest_momentum_zero(u, q, j, C):
    v2 = -dot(u, u)
    v = s.sqrt(v2)

    uq = dot(u, q)
    w = q + uq*u/v2
    w2 = dot(w, w)

    uj = dot(u, j)
    wj = dot(w, j)

    R = j + uj*u/v2 - wj*w/w2

    return -2*C*(eta*R)[0]/(v*w2)

u0 = s.Matrix([5, 3, 0, 0])
q0 = s.Matrix([1, 2, 2, 1])
j0 = s.Matrix([4, 1, 3, 6])
a0 = s.Matrix([2, -1, 4, 3])
b0 = s.Matrix([1, 3, -2, 5])

A0, A1 = s.Integer(3), s.Integer(2)
C0, C1, C2 = (
    s.Integer(2),
    s.Integer(1),
    s.Integer(3)
)

basis = s.Matrix([1, 0, 0, 0])

# First contribution: partial L / partial u^0.
stage("STAGE 1: partial-u derivative")

gu = s.diff(
    lagrangian(
        u0 + e*basis, q0, j0, A0, C0
    ),
    e
).subs(e, 0)

gu = s.simplify(gu)

stage("STAGE 1 COMPLETE")

# Second contribution:
# d/dlambda (partial L / partial q^0).
#
# Use first-order trajectories because this derivative
# needs only u'=q, q'=j, j'=a.
stage("STAGE 2: mixed directional derivative")

ut1 = u0 + q0*t
qt1 = q0 + j0*t
jt1 = j0 + a0*t

At1 = A0 + A1*t
Ct1 = C0 + C1*t

expr_q = lagrangian(
    ut1, qt1 + e*basis, jt1, At1, Ct1
)

# Differentiate with respect to the perturbation first.
# Then evaluate that perturbation at zero before
# taking the worldline derivative.
partial_q = s.diff(expr_q, e).subs(e, 0)

gqdot = s.diff(partial_q, t).subs(t, 0)
gqdot = s.simplify(gqdot)

stage("STAGE 2 COMPLETE")

# Third contribution:
# d^2/dlambda^2 (partial L / partial j^0).
#
# Here we use the analytically derived highest momentum,
# independently of the R22B polynomial engine.
stage("STAGE 3: second derivative of P_0")

ut2 = u0 + q0*t + j0*t**2/2
qt2 = q0 + j0*t + a0*t**2/2
jt2 = j0 + a0*t + b0*t**2/2

Ct2 = C0 + C1*t + C2*t**2/2

Pt0 = highest_momentum_zero(
    ut2, qt2, jt2, Ct2
)

Pddot = s.diff(Pt0, t, 2).subs(t, 0)
Pddot = s.simplify(Pddot)

stage("STAGE 3 COMPLETE")

p0_component = s.simplify(gu - gqdot + Pddot)

expected = s.Rational(
    1797334885453,
    504122720256
)

residual = s.simplify(p0_component - expected)

print()
print("PARTIAL-u =", gu)
print("q-GRADIENT DERIVATIVE =", gqdot)
print("SECOND P DERIVATIVE =", Pddot)
print("INDEPENDENT p0[0] =", p0_component)
print("R22B p0[0] =", expected)
print("RESIDUAL =", residual)

assert residual == 0, (
    "STOP: independent momentum component disagrees"
)

print()
stage("R24A PASS")
print("FIRST COMPONENT: EXACT MATCH")
print("REMAINING SEVEN COMPONENT COMPARISONS: OPEN")
print("FULL ACTION MOMENTUM: OPEN")
print("DIRAC SYSTEM: OPEN")
