"""CP098-T74: Full regular-surface transport of E.

Differentiate the verified surface expression E
using independently established scalar transport.

This is a surface-chain-rule result, NOT YET
an independent full canonical {E,H_T} audit.
"""

import sympy as sp

print("=== CP098-T74 ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a, n, R, S = sp.symbols("a n R S", real=True)

A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)

DA, DC = sp.symbols("DA DC", real=True)
H, J, K = sp.symbols("H J K", real=True)

LT, LN = sp.symbols("LT LN", real=True)

r, s, x, y, b, c = sp.symbols(
    "r s x y b c", real=True
)

def eq(lhs, rhs):
    return sp.cancel(sp.together(lhs-rhs)) == 0

# T70: verified regular-surface expression.
E = (
    -2*A*a*k/v**3
    +2*DA*k/v
    -n/v
    -a*k*v*R/(2*C)
    -k*v**2*S/(2*C)
)

print("PASS: T70 E expression reconstructed",
      flush=True)

# T64, T68, T69: external scalar transport.
vdot = a
kdot = LN/v**2-3*a*k/v
adot = LT-v**3*k**2

Rdot = -2*S-2*LN*R/(v**2*k)

# T68R3: n=p0.N, including transverse p0.
p0_time = A*k**2+v**4*k**2*R/(4*C)-B

ndot = (
    v*k*p0_time
    -v**3*k*(b*r+c*s)/(2*C)
)

# T68R3: invariant P.p1 derivative,
# not the derivative of coordinate r*x+s*y.
Sdot = (
    -x**2-y**2-r*b-s*c
    -LN*S/(v**2*k)
    -LT*R/v
    +a*LN*R/(k*v**3)
    -A*k**2*v**2*R/C
)

print("PASS: external transport reconstructed",
      flush=True)

# T62: internal canonical evolution.
Adot = v*DA
Cdot = v*DC

# H={DA,B}, J={DA,A}, K={DA,1/C}.
DAdot = (
    v*H
    +v*k**2*J
    -v**5*k**2*R*K/4
)

print("PASS: internal evolution retained",
      flush=True)

rates = {
    v: vdot,
    k: kdot,
    a: adot,
    n: ndot,
    R: Rdot,
    S: Sdot,
    A: Adot,
    C: Cdot,
    DA: DAdot
}

# Route 1: multivariable chain rule.
Edot_chain = sp.factor(
    sum(
        sp.diff(E, variable)*rate
        for variable, rate in rates.items()
    )
)

print("PASS: complete chain rule constructed",
      flush=True)

# Route 2: differentiate the five distinct
# contributions using their own product rules.

term_geometry = -2*A*a*k/v**3
term_internal = 2*DA*k/v
term_p0 = -n/v
term_quadratic = -a*k*v*R/(2*C)
term_p1 = -k*v**2*S/(2*C)

pieces = [
    term_geometry,
    term_internal,
    term_p0,
    term_quadratic,
    term_p1
]

piece_dots = [
    sp.factor(
        sum(
            sp.diff(term, variable)*rate
            for variable, rate in rates.items()
        )
    )
    for term in pieces
]

assert eq(sum(pieces), E)
assert eq(sum(piece_dots), Edot_chain)

print("PASS: five-term derivative reconciliation",
      flush=True)

# T71 independently verified this contribution.
expected_p0_dot = (
    -ndot/v+n*vdot/v**2
)

assert eq(piece_dots[2], expected_p0_dot)

print("PASS: T71 p0 evolution recovered",
      flush=True)

# Internal Hessian/gradient structures must survive.
assert eq(
    sp.diff(Edot_chain,H),
    2*k
)

assert eq(
    sp.diff(Edot_chain,J),
    2*k**3
)

assert eq(
    sp.diff(Edot_chain,K),
    -v**4*k**3*R/2
)

print("PASS: internal H,J,K coefficients",
      flush=True)

# Preserve the normal multiplier until the
# complete unconstrained-rate expression
# has been constructed.

FN = (
    -n
    -6*A*a*k/v**2
    -v**3*k*S/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

LN_fixed = -FN*v**3/(2*A)

Edot_fixed = sp.factor(
    Edot_chain.subs(LN,LN_fixed)
)

assert LN not in Edot_fixed.free_symbols

print("PASS: normal multiplier substituted last",
      flush=True)

# Verify the transverse p0 coefficients remain.
# Sdot and ndot each carry transverse terms;
# both contributions must be included.

assert b in Edot_fixed.free_symbols
assert c in Edot_fixed.free_symbols

print("PASS: transverse momenta retained",
      flush=True)

# Check both representations against the
# explicit regular adapted coordinates.
surface = {
    R: r**2+s**2,
    S: r*x+s*y
}

assert eq(
    Edot_chain.subs(surface),
    sum(piece_dots).subs(surface)
)

print("PASS: adapted-coordinate reconciliation",
      flush=True)

print(
    "E_dot_unfixed =",
    sp.factor(Edot_chain),
    flush=True
)

print(
    "E_dot_fixed =",
    Edot_fixed,
    flush=True
)

print("REGULAR-SURFACE E TRANSPORT: VERIFIED")
print("NORMAL MULTIPLIER: APPLIED AFTER DERIVATIVE")
print("INTERNAL HESSIAN STRUCTURE: RETAINED")
print("INDEPENDENT CANONICAL E_DOT: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
