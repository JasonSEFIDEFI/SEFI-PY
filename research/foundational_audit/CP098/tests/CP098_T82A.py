"""CP098-T82A: Exact canonical Jacobi implementation audit.

Purpose:
    Verify canonical Poisson-bracket signs and Jacobi identity
    using exact polynomial test functions.

Scope:
    Full 22-dimensional canonical phase space.
    No constraints imposed.
    No numerical differentiation.
    No physical validation claim.
"""

import sympy as sp

print("=== CP098-T82A ===", flush=True)

y = sp.symbols("y0:4", real=True)
u = sp.symbols("u0:4", real=True)
q = sp.symbols("q0:4", real=True)

p0 = sp.symbols("p0_0:4", real=True)
p1 = sp.symbols("p1_0:4", real=True)
P = sp.symbols("P0:4", real=True)

I, pI = sp.symbols("I pI", real=True)

pairs = (
    list(zip(y, p0))
    + list(zip(u, p1))
    + list(zip(q, P))
    + [(I, pI)]
)

assert len(pairs) == 13
# 4 + 4 + 4 + 1 = 13 canonical pairs,
# hence 26-dimensional full phase space.

def bracket(f, g):
    result = sp.S.Zero

    for coordinate, momentum in pairs:
        result += (
            sp.diff(f, coordinate)
            * sp.diff(g, momentum)
            - sp.diff(f, momentum)
            * sp.diff(g, coordinate)
        )

    return sp.expand(result)

# Canonical sign and normalization checks.
for coordinate, momentum in pairs:
    assert bracket(coordinate, momentum) == 1
    assert bracket(momentum, coordinate) == -1
    assert bracket(coordinate, coordinate) == 0
    assert bracket(momentum, momentum) == 0

print("PASS: canonical pair signs", flush=True)

# Cross-sector commutativity.
assert bracket(y[0], p1[0]) == 0
assert bracket(u[0], P[0]) == 0
assert bracket(q[0], p0[0]) == 0
assert bracket(I, p0[0]) == 0

print("PASS: cross-sector brackets", flush=True)

# Three independent polynomial functions.
# Include mixed dependencies across all sectors.
f = (
    y[0]*p0[1]
    + u[1]*p1[2]
    + q[2]*P[3]
    + I*pI**2
    + y[2]*u[0]*P[1]
)

g = (
    p0[0]*u[2]
    + p1[1]*q[0]
    + P[2]*y[3]
    + I**2*pI
    + q[1]*p0[2]*p1[3]
)

h = (
    y[1]*p1[0]
    + u[3]*P[0]
    + q[3]*p0[3]
    + I**3
    + pI*p0[1]*P[2]
)

jacobi = (
    bracket(f, bracket(g, h))
    + bracket(g, bracket(h, f))
    + bracket(h, bracket(f, g))
)

residual = sp.expand(jacobi)

if residual != 0:
    print("FAIL: Jacobi residual:", residual)
    raise SystemExit(1)

print("PASS: exact polynomial Jacobi identity")
print("PASS: all canonical sectors included")
print("GENERAL CP098 CONSTRAINT JACOBI AUDIT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T82A COMPLETE")
