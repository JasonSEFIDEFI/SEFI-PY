"""CP098-T21A: Directional derivative optimization.

Independent symbolic calculation of the (u,p1)
contribution to {psi_T,chi_N}.

Regular adapted family:
v>0, k>0, a,t,n real.

Does not establish complete Dirac closure.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
a, t, n = sp.symbols("a t n", real=True)
e = sp.symbols("e", real=True)

eta = sp.diag(-1, 1, 1, 1)

def mdot(x, y):
    return (x.T * eta * y)[0]

def contract(p, x):
    return (p.T * x)[0]

# Adapted geometry at epsilon=0.
u = sp.Matrix([v, e, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

P = sp.Matrix([0, 0, 1, 0])

p0 = sp.Matrix([
    k**2 + v**4*k**2/4 - t**2/2,
    n,
    0,
    0,
])

A = sp.Integer(1)
C = sp.Integer(1)
B = t**2/2

# Use squared speed to avoid premature expansion
# of nested square-root expressions.
V2 = -mdot(u, u)

T = u/sp.sqrt(V2)
Tq = mdot(T, q)

w = q + T*Tq
w2 = mdot(w, w)

kappa2 = w2/V2**2
P2 = mdot(P, P)

psi = (
    contract(p0, T)
    + B
    - A*kappa2
    - V2**2*kappa2*P2/(4*C)
)

print("=== CP098-T21A ===", flush=True)

# At epsilon=0, N=(0,1,0,0).
# Therefore d/de is precisely N^mu d/du^mu.
#
# Differentiate before setting epsilon=0.

result = sp.diff(psi, e).subs(e, 0)
result = sp.factor(result)

print("Directional derivative =", result, flush=True)

# Independently reconstruct the expected expression
# from the five-parameter family.
#
# The canonical u,p1 contribution is the derivative
# calculated above; no closure condition is assumed.

expected = (
    n/v
    + a*(k**2/v + v**2*k**2/4)
    + 4*k**2*a/v**3
)

difference = sp.factor(result - expected)

print("Diagnostic candidate =", sp.factor(expected), flush=True)
print("Candidate difference =", difference, flush=True)

if difference == 0:
    print("PASS: diagnostic candidate matches")
else:
    print("CANDIDATE REJECTED: retain direct derivative")

# Independently verify a known geometry from T19.
known = {v: 3, k: 2, a: 1}
known_result = sp.factor(result.subs(known))

print("T19 geometry derivative =", known_result)

assert sp.simplify(
    known_result - (n/3 + sp.Rational(85, 27))
) == 0, "Mismatch with T19 u,p1 contribution"

print("PASS: agrees with T19")
print("GENERAL CLOSURE: NOT ESTABLISHED")
