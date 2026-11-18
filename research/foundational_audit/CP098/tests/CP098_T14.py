"""CP098-T14: General tangent-secondary Poisson bracket.

Calculates {psi_T, chi_T} without imposing constraints.

The candidate is tested independently. A mismatch is
reported rather than suppressed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

def metric(x, y):
    return (x.T * eta * y)[0]

def contract(p, x):
    return (p.T * x)[0]

us = sp.symbols("u0:4", real=True)
qs = sp.symbols("q0:4", real=True)
Ps = sp.symbols("P0:4", real=True)
rs = sp.symbols("r0:4", real=True)
ss = sp.symbols("s0:4", real=True)

u = sp.Matrix(us)
q = sp.Matrix(qs)
P = sp.Matrix(Ps)
p1 = sp.Matrix(rs)
p0 = sp.Matrix(ss)

A, B = sp.symbols("A B", real=True)
C = sp.symbols("C", real=True, nonzero=True)

v = sp.sqrt(-metric(u, u))
T = u/v

Tq = metric(T, q)
w = q + T*Tq
wn2 = metric(w, w)

kappa2 = wn2/v**4
P2 = metric(P, P)

psiT = (
    contract(p0, T)
    + B
    - A*kappa2
    - v**4*kappa2*P2/(4*C)
)

chiT = contract(p1, T)

# The bracket simplifies because chiT is independent
# of q and P. Only the (u,p1) canonical pair contributes.

# d(chiT)/d(p1_mu) = T^mu
# d(chiT)/d(u^mu) = p1_nu d(T^nu)/d(u^mu)

# psiT is independent of p1, so:
# {psiT,chiT} = sum_mu T^mu d(psiT)/d(u^mu)

print("=== CP098-T14: GENERAL TANGENT BRACKET ===", flush=True)

# Radial derivative under u -> alpha*u.
# v scales as alpha, T is unchanged,
# kappa^2 scales as alpha^-4 when q is fixed,
# and v^4*kappa^2 = w^2 is unchanged.
#
# Thus psiT's only radial dependence is -A*kappa^2.
# T^mu d/d(u^mu) = (1/v) u^mu d/d(u^mu).

bracket_derived = 4*A*kappa2/v

print("Derived bracket:", bracket_derived, flush=True)

# Independent candidate from the previous script.
candidate = (
    -contract(p0, T)/v
    + 4*A*kappa2/v
    + v**3*kappa2*P2/C
)

difference = sp.simplify(
    bracket_derived - candidate
)

print("Candidate difference:", difference, flush=True)

if difference == 0:
    print("PASS: candidate identity", flush=True)
else:
    print("CORRECTION REQUIRED: candidate does not match", flush=True)

# Verify the radial derivative directly at an exact
# regular timelike point without expensive global simplify.
point = dict(zip(
    us + qs + Ps + rs + ss,
    (
        3,0,0,0,
        1,18,0,0,
        0,0,1,0,
        0,0,0,0,
        83,2,0,0
    )
))

point.update({A: sp.Integer(1), B: sp.Integer(3), C: sp.Integer(1)})

radial = sum(
    T[i]*sp.diff(psiT, us[i])
    for i in range(4)
)

actual = sp.simplify(radial.subs(point))
expected = sp.simplify(bracket_derived.subs(point))

print("Direct derivative at test point:", actual, flush=True)
print("Derived expression at test point:", expected, flush=True)

assert actual == expected

print("PASS: direct derivative agrees", flush=True)

print()
print("LIMITATIONS:")
print("- General identity uses homogeneity argument")
print("- Direct derivative checked at one regular point")
print("- Full Dirac closure remains unresolved")
print("- Physical stability remains unproved")
