"""CP098-T83C: Exact canonical transport of a=-T.q.

Independent of T74.
Adapted Lorentz frame, unrestricted momenta.
Differentiate before imposing constraints.
"""

import sympy as sp

print("=== CP098-T83C ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, LT, LN = sp.symbols("a LT LN", real=True)
C = sp.symbols("C", real=True, nonzero=True)

P0, P1, P2, P3 = sp.symbols(
    "P0 P1 P2 P3", real=True
)

eps = sp.symbols("eps", real=True)
eta = sp.diag(-1, 1, 1, 1)

def dot(X, Y):
    return (X.T*eta*Y)[0]

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([P0, P1, P2, P3])

T = u/v
N = sp.Matrix([0, 1, 0, 0])

# Canonical flow with multipliers held fixed.
udot = q

qdot = (
    -v**5*k**2*eta*P/(2*C)
    +LT*T+LN*N
)

def scalar_a(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    tangent = uu/speed
    return -dot(tangent, qq)

actual = sp.diff(
    scalar_a(u+eps*udot, q+eps*qdot),
    eps
).subs(eps, 0)

print("PASS: canonical directional derivative built",
      flush=True)

# At this adapted point, phi_T=P.T=P0.
phiT = P0

expected = (
    LT-v**3*k**2
    +v**5*k**2*phiT/(2*C)
)

residual = sp.cancel(
    sp.together(actual-expected)
)

if residual != 0:
    raise AssertionError(
        f"STOP: off-surface residual = {residual}"
    )

print("PASS: exact off-surface a transport",
      flush=True)

surface_rate = sp.simplify(
    actual.subs(P0, 0)
)

T74_rate = LT-v**3*k**2

if sp.simplify(surface_rate-T74_rate) != 0:
    raise AssertionError(
        "STOP: T74 surface rate mismatch"
    )

print("PASS: T74 regular-surface a transport",
      flush=True)

correction = sp.simplify(actual-T74_rate)

if sp.simplify(
    correction-v**5*k**2*phiT/(2*C)
) != 0:
    raise AssertionError(
        "STOP: off-surface correction mismatch"
    )

print("PASS: off-surface correction identified")
print("PASS: constraints applied after differentiation")
print("PASS: generic internal C retained")
print("A TRANSPORT: EXACTLY VERIFIED")
print("REMAINING SCALAR TRANSPORT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T83C COMPLETE")
