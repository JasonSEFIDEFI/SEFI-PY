"""CP098-T85C: Off-surface compact bracket audit.

Compare T85B's compact geometric E with T70's
independent adapted canonical bracket.

All canonical momenta remain unrestricted.
No constraint substitutions or multiplier fixing.
"""

import sympy as sp

print("=== CP098-T85C ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)

A, C = sp.symbols("A C", real=True, nonzero=True)
DA = sp.symbols("DA", real=True)

eta = sp.diag(-1, 1, 1, 1)

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def covdot(p, z):
    return (p.T*eta*z)[0]

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("j0:4", real=True))
p0 = sp.Matrix(sp.symbols("z0:4", real=True))

speed = sp.sqrt(-dot(u, u))
T = u/speed
w = q+T*dot(T, q)
L = sp.sqrt(dot(w, w))
N = w/L

assert sp.simplify(speed-v) == 0
assert sp.simplify(L-v**2*k) == 0

print("PASS: unrestricted adapted geometry",
      flush=True)

# Compact formula used by T85B.

E_compact = (
    -L/(2*C)*(
        covdot(p1, P)
        +pair(p1, T)*covdot(T, P)
        -pair(p1, N)*covdot(N, P)
    )
    -A*L**2*covdot(N, P)/(C*speed**3)
    -pair(p0, N)/speed
    +2*A*L*dot(T, q)/speed**5
    +L*dot(T, q)*dot(P, P)/(2*C*speed)
    +2*L*DA/speed**3
)

P0, P1, P2, P3 = P
j0, j1, j2, j3 = p1
z0, z1, z2, z3 = p0

P_squared = -P0**2+P1**2+P2**2+P3**2

# Independent T70 canonical bracket,
# prior to any momentum constraints.

E_T70 = -(
    4*A*C*a*k
    +2*A*P1*k**2*v**4
    -4*C*DA*k*v**2
    +2*C*z1*v**2
    +a*k*v**4*P_squared
    +k*v**5*(P2*j2+P3*j3)
)/(2*C*v**3)

residual = sp.cancel(
    sp.together(E_compact-E_T70)
)

print("OFF_SURFACE_RESIDUAL =", residual,
      flush=True)

assert residual == 0, (
    "STOP: compact bracket differs from T70 "
    "before imposing momentum constraints"
)

print("PASS: exact off-surface agreement",
      flush=True)

# Explicitly retain the unconstrained
# primary momentum components.

for variable in (P0, P1, P2, P3):
    assert variable in E_compact.free_symbols

print("PASS: all P components retained",
      flush=True)

# Confirm that this comparison has not
# silently imposed secondary constraints.

assert j1 in E_compact.free_symbols
assert j2 in E_compact.free_symbols
assert j3 in E_compact.free_symbols
assert z1 in E_compact.free_symbols

print("PASS: unconstrained p1 and p0 retained",
      flush=True)

print("OFF-SURFACE BRACKET CONSISTENCY: VERIFIED")
print("FULL GENERAL-POINT PB AUDIT: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T85C COMPLETE")
