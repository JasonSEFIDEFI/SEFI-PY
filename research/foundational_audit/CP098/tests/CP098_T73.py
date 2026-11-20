"""CP098-T73: Independent preservation of psi_T.

Verify the regular-surface Hamiltonian derivative
of psi_T from canonical external transport and
internal Poisson brackets.

Do not assume psi_T preservation in the derivation.
"""

import sympy as sp

print("=== CP098-T73 ===", flush=True)

v, k = sp.symbols("v k", positive=True)

a, r, s, x, y, n, b, c = sp.symbols(
    "a r s x y n b c", real=True
)

A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

DA, DC = sp.symbols("DA DC", real=True)
LT, LN = sp.symbols("LT LN", real=True)

def eq(lhs, rhs):
    return sp.cancel(sp.together(lhs-rhs)) == 0

R = r**2+s**2
S = r*x+s*y

# Regular adapted constraint surface.
p0_time = (
    A*k**2
    +v**4*k**2*R/(4*C)
    -B
)

psi = (
    p0_time+B-A*k**2
    -v**4*k**2*R/(4*C)
)

assert eq(psi,0)
print("PASS: psi_T surface reconstructed",
      flush=True)

# External transport from T64-T68.
vdot = a
kdot = LN/v**2-3*a*k/v

Rdot = (
    -2*S
    -2*LN*R/(v**2*k)
)

# p0 is conserved, but T moves:
# d(p0.T)/dt = p0.Tdot = v*k*n.
p0Tdot = v*k*n

# Internal evolution from canonical brackets:
#
# DA = {A,B}
# DC = {C,B}
#
# A_dot = v*DA
# C_dot = v*DC
#
# B_dot = v*k^2 {B,A}
#       -v^5*k^2*R/4 {B,1/C}
#
# {B,A}=-DA
# {B,1/C}=DC/C^2.

Adot = v*DA
Cdot = v*DC

Bdot = (
    -v*k**2*DA
    -v**5*k**2*R*DC/(4*C**2)
)

print("PASS: internal evolution constructed",
      flush=True)

# Independently differentiate each term of psi_T.

Akap_dot = (
    Adot*k**2
    +2*A*k*kdot
)

Qdot = (
    4*v**3*vdot*k**2*R/(4*C)
    +v**4*2*k*kdot*R/(4*C)
    +v**4*k**2*Rdot/(4*C)
    -v**4*k**2*R*Cdot/(4*C**2)
)

psi_dot = sp.factor(
    p0Tdot
    +Bdot
    -Akap_dot
    -Qdot
)

print("PASS: psi_T derivative constructed",
      flush=True)

# Reconstruct the T69 normal constraint
# preservation equation independently.

FN = (
    -n
    -6*A*a*k/v**2
    -v**3*k*S/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

chiN_dot = FN+2*A*LN/v**3

assert eq(
    psi_dot,
    -v*k*chiN_dot
)

print("PASS: secondary preservation identity",
      flush=True)

# Fix lambda_N using the T69 equation.
LN_fixed = -FN*v**3/(2*A)

assert eq(
    psi_dot.subs(LN,LN_fixed),
    0
)

print("PASS: psi_T preserved with fixed lambda_N",
      flush=True)

# Check the internal C-gradient cancellation.
assert eq(sp.diff(psi_dot,DC),0)

print("PASS: DC contributions cancel",
      flush=True)

# The tangential multiplier must not appear.
assert LT not in psi_dot.free_symbols

print("PASS: tangential multiplier independence",
      flush=True)

# Transverse p0 components are retained in
# the regular family but do not enter this
# particular preservation equation.
assert b not in psi_dot.free_symbols
assert c not in psi_dot.free_symbols

print("PASS: transverse p0 independence",
      flush=True)

assert eq(
    sp.diff(psi_dot,LN),
    -2*A*k/v**2
)

print("PASS: normal multiplier coefficient",
      flush=True)

print("psi_T_dot =",psi_dot,flush=True)
print("lambda_N_fixed =",sp.factor(LN_fixed),
      flush=True)

print("SECONDARY CONSTRAINT PRESERVATION: VERIFIED")
print("REGULAR SURFACE TANGENCY: CONSISTENT")
print("FULL E EVOLUTION: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
