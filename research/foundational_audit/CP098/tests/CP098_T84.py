"""CP098-T84: Simultaneous regular-surface tangency.

Assemble five constraint derivatives from canonical
transport components and internal evolution.

Regular sector only.
Does not establish the unrestricted E bracket.
"""

import sympy as sp

print("=== CP098-T84 ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, r, s, x, y, n, b, c = sp.symbols(
    "a r s x y n b c", real=True
)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
DA, DC, LT, LN = sp.symbols(
    "DA DC LT LN", real=True
)

eta = sp.diag(-1, 1, 1, 1)

def pair(p, X):
    return (p.T*X)[0]

def dot(X, Y):
    return (X.T*eta*Y)[0]

def eq(f, g=0):
    return sp.cancel(sp.together(f-g)) == 0

def rate_pair(p, pdot, X, Xdot):
    return pair(pdot, X)+pair(p, Xdot)

R = r**2+s**2

T = sp.Matrix([1, 0, 0, 0])
N = sp.Matrix([0, 1, 0, 0])

P = sp.Matrix([0, 0, r, s])
p1 = sp.Matrix([0, -2*A*k/v, x, y])

p0 = sp.Matrix([
    A*k**2+v**4*k**2*R/(4*C)-B,
    n, b, c
])

Tdot = sp.Matrix([0, v*k, 0, 0])

Ndot = sp.Matrix([
    v*k, 0,
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

Pdot = sp.Matrix([
    0,
    v**3*k*R/(2*C),
    -x-LN*r/(v**2*k),
    -y-LN*s/(v**2*k)
])

p1dot = sp.Matrix([
    2*A*k**2,
    2*A*a*k/v**2-n-a*k*v**2*R/(2*C),
    a*LN*r/(k*v**3)-LT*r/v-b,
    a*LN*s/(k*v**3)-LT*s/v-c
])

p0dot = sp.zeros(4, 1)

vdot = a
kdot = LN/v**2-3*a*k/v
Adot = v*DA
Cdot = v*DC

Bdot = (
    -v*k**2*DA
    -v**5*k**2*R*DC/(4*C**2)
)

Rdot = 2*pair(Pdot, eta*P)

# Primary constraints.
phiTdot = rate_pair(P, Pdot, T, Tdot)
phiNdot = rate_pair(P, Pdot, N, Ndot)

# Tangential secondary.
chiTdot = rate_pair(p1, p1dot, T, Tdot)

# Normal secondary.
chiNdot = (
    rate_pair(p1, p1dot, N, Ndot)
    +2*Adot*k/v
    +2*A*kdot/v
    -2*A*k*vdot/v**2
)

# Hamiltonian secondary.
psiTdot = (
    rate_pair(p0, p0dot, T, Tdot)
    +Bdot
    -Adot*k**2
    -2*A*k*kdot
    -v**3*vdot*k**2*R/C
    -v**4*k*kdot*R/(2*C)
    -v**4*k**2*Rdot/(4*C)
    +v**4*k**2*R*Cdot/(4*C**2)
)

rates = {
    "phi_T": phiTdot,
    "phi_N": phiNdot,
    "chi_T": chiTdot,
    "chi_N": chiNdot,
    "psi_T": psiTdot,
}

print("PASS: all five derivatives assembled",
      flush=True)

# The first three must vanish without
# fixing either primary multiplier.

for name in ("phi_T", "phi_N", "chi_T"):
    assert eq(rates[name]), (
        f"STOP: {name} not preserved"
    )
    print(f"PASS: {name} identically preserved",
          flush=True)

# Normal secondary determines LN.
d = sp.diff(chiNdot, LN)

assert eq(d, 2*A/v**3)
assert eq(sp.diff(chiNdot, LT))

FN = sp.factor(chiNdot.subs(LN, 0))
LN_fixed = -FN/d

assert eq(chiNdot.subs(LN, LN_fixed))

print("PASS: normal multiplier uniquely fixed",
      flush=True)

# Independent secondary compatibility.
assert eq(psiTdot, -v*k*chiNdot)

print("PASS: secondary compatibility identity",
      flush=True)

# Simultaneous preservation, not merely
# separate solvability of equations.

for name, expression in rates.items():
    residual = sp.factor(
        expression.subs(LN, LN_fixed)
    )
    if not eq(residual):
        print(name, "RESIDUAL =", residual,
              flush=True)
        raise AssertionError(
            f"STOP: simultaneous tangency: {name}"
        )
    print(f"PASS: simultaneous {name} preservation",
          flush=True)

# Verify LT is not fixed by any of these
# surface consistency equations.

for name, expression in rates.items():
    assert eq(sp.diff(expression, LT)), (
        f"STOP: LT dependence in {name}"
    )

print("PASS: lambda_T remains unrestricted "
      "by these five equations", flush=True)

print("F_N =", FN, flush=True)
print("LAMBDA_N =", sp.factor(LN_fixed),
      flush=True)

print("SIMULTANEOUS REGULAR TANGENCY: VERIFIED")
print("FULL PHASE-SPACE PB AUDIT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T84 COMPLETE")
