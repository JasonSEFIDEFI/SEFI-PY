"""CP098-T69: Independent normal-secondary preservation.

Reconstruct chi_N_dot using verified canonical
momentum evolution, moving-frame transport,
and internal coupling evolution.

Compare against the earlier F_N + d*lambda_N
preservation equation.

Regular sector only: v>0, kappa>0, A!=0, C!=0.
"""

import sympy as sp

print("=== CP098-T69 ===", flush=True)

v, k = sp.symbols("v k", positive=True)

a, r, s, x, y, n, b, c = sp.symbols(
    "a r s x y n b c", real=True
)

A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

DA = sp.symbols("DA", real=True)

LT, LN = sp.symbols(
    "lambda_T lambda_N", real=True
)

def pair(p, X):
    return (p.T * X)[0]

def equal(lhs, rhs):
    return sp.cancel(sp.together(lhs-rhs)) == 0

T = sp.Matrix([1, 0, 0, 0])
N = sp.Matrix([0, 1, 0, 0])

R = r**2+s**2
S = r*x+s*y

P = sp.Matrix([0, 0, r, s])

p1 = sp.Matrix([
    0,
    -2*A*k/v,
    x,
    y
])

p0 = sp.Matrix([
    A*k**2+v**4*k**2*R/(4*C)-B,
    n,
    b,
    c
])

p1dot = sp.Matrix([
    2*A*k**2,
    2*A*a*k/v**2-n-a*k*v**2*R/(2*C),
    a*LN*r/(k*v**3)-LT*r/v-b,
    a*LN*s/(k*v**3)-LT*s/v-c
])

Ndot = sp.Matrix([
    v*k,
    0,
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

vdot = a
kdot = LN/v**2-3*a*k/v

# From T62/T67:
# {A,H_T} = v*DA.
Adot = v*DA

assert equal(pair(p1,N), -2*A*k/v)
assert equal(pair(P,N), 0)
assert equal(pair(p1,T), 0)

print("PASS: regular constraints reconstructed",
      flush=True)

# Independently account for each contribution
# to the total derivative of chi_N.

momentum_term = sp.factor(pair(p1dot,N))
frame_term = sp.factor(pair(p1,Ndot))

coupling_term = sp.factor(2*k*Adot/v)

curvature_term = sp.factor(2*A*kdot/v)

speed_term = sp.factor(-2*A*k*vdot/v**2)

chiNdot = sp.factor(
    momentum_term
    + frame_term
    + coupling_term
    + curvature_term
    + speed_term
)

print("p1_dot.N =", momentum_term, flush=True)
print("p1.N_dot =", frame_term, flush=True)
print("A contribution =", coupling_term, flush=True)
print("kappa contribution =", curvature_term, flush=True)
print("v contribution =", speed_term, flush=True)

# Independently reconstruct the previously
# derived normal preservation equation.

d = 2*A/v**3

FN = (
    -n
    -6*A*a*k/v**2
    -v**3*k*S/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

expected = sp.factor(FN+d*LN)

assert equal(chiNdot,expected)

print("PASS: independent F_N reconstruction",
      flush=True)

assert equal(
    sp.diff(chiNdot,LN),
    d
)

print("PASS: normal multiplier coefficient",
      flush=True)

assert equal(
    sp.diff(chiNdot,LT),
    0
)

print("PASS: tangential multiplier independence",
      flush=True)

# The normal multiplier is fixed by preservation.
LN_fixed = -FN/d

assert equal(
    chiNdot.subs(LN,LN_fixed),
    0
)

print("PASS: chi_N preserved with fixed lambda_N",
      flush=True)

# A remains nonzero in this regular sector,
# so the multiplier equation has a unique
# solution at this stage.

assert equal(
    sp.diff(expected,LN),
    2*A/v**3
)

print("PASS: regular-sector multiplier uniqueness",
      flush=True)

# Confirm that transverse p0 components b,c
# do not enter this particular preservation
# equation even though they are retained
# in the complete regular family.

assert b not in chiNdot.free_symbols
assert c not in chiNdot.free_symbols

print("PASS: transverse p0 cancellation/absence",
      flush=True)

print("chi_N_dot =", chiNdot, flush=True)
print("lambda_N_fixed =", sp.factor(LN_fixed),
      flush=True)

print("NORMAL SECONDARY PRESERVATION: VERIFIED")
print("T50 NORMAL MULTIPLIER: REPRODUCED")
print("T64-T68 TRANSPORT: CONSISTENT")
print("FULL E EVOLUTION: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
