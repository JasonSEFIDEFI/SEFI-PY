"""CP098-T83B: Exact canonical curvature transport.

Independent of T74.
Adapted Lorentz frame, unrestricted momenta.
Canonical differentiation precedes constraints.
"""

import sympy as sp

print("=== CP098-T83B ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, LT, LN = sp.symbols("a LT LN", real=True)
C = sp.symbols("C", real=True, nonzero=True)

P0, P1, P2, P3 = sp.symbols(
    "P0 P1 P2 P3", real=True
)

eps = sp.symbols("eps", real=True)
eta = sp.diag(-1, 1, 1, 1)

def dot(X, Y):
    return (X.T * eta * Y)[0]

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([P0, P1, P2, P3])

T = sp.Matrix([1, 0, 0, 0])
N = sp.Matrix([0, 1, 0, 0])

# Canonical equations before imposing constraints:
#
# udot = dHT/dp1 = q
#
# qdot = dHT/dP
#      = -v^5*k^2*eta*P/(2*C)
#        +LT*T+LN*N
#
# The multiplier parameters are held fixed.

udot = q
qdot = (
    -v**5*k**2*eta*P/(2*C)
    +LT*T+LN*N
)

def curvature(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    tangent = uu/speed
    w = qq+tangent*dot(tangent, qq)
    return sp.sqrt(dot(w, w))/speed**2

# Directional derivative of the exact
# geometric curvature definition.
#
# The flow vectors are evaluated at the
# initial point, as required for a
# first directional derivative.

ue = u+eps*udot
qe = q+eps*qdot

expression = curvature(ue, qe)

actual = sp.diff(expression, eps).subs(eps, 0)

print("PASS: canonical directional derivative built",
      flush=True)

# Off-surface prediction.
phiN = P1

expected = (
    LN/v**2
    -3*a*k/v
    -v**3*k**2*phiN/(2*C)
)

residual = sp.cancel(
    sp.together(actual-expected)
)

if residual != 0:
    raise AssertionError(
        f"STOP: off-surface residual = {residual}"
    )

print("PASS: exact off-surface curvature transport",
      flush=True)

# Apply the normal momentum constraint
# only after the canonical derivative.
surface_rate = sp.simplify(
    actual.subs(P1, 0)
)

T74_rate = LN/v**2-3*a*k/v

if sp.simplify(surface_rate-T74_rate) != 0:
    raise AssertionError(
        "STOP: T74 surface curvature mismatch"
    )

print("PASS: T74 regular-surface curvature rate",
      flush=True)

# Verify the off-surface correction is
# proportional to the normal constraint.
correction = sp.simplify(
    actual-T74_rate
)

if sp.simplify(
    correction+v**3*k**2*phiN/(2*C)
) != 0:
    raise AssertionError(
        "STOP: normal-constraint correction mismatch"
    )

print("PASS: off-surface correction identified")
print("PASS: constraints applied after differentiation")
print("PASS: generic internal C retained")
print("CURVATURE TRANSPORT: EXACTLY VERIFIED")
print("REMAINING SCALAR TRANSPORT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T83B COMPLETE")
