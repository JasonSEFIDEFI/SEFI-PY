"""CP098-T22: Full q,P contribution with transverse momenta.

Tests the canonical (q,P) contribution to {psi_T,chi_N}
at an adapted regular geometry, retaining all four
q and P components through differentiation.
"""

import sympy as sp

v, k, A, C = sp.symbols(
    "v k A C", positive=True
)
x, y, r, s = sp.symbols("x y r s", real=True)

q = sp.Matrix(sp.symbols("q0:4", real=True))
P = sp.Matrix(sp.symbols("P0:4", real=True))

eta = sp.diag(-1, 1, 1, 1)
T = sp.Matrix([1, 0, 0, 0])

def mdot(a, b):
    return (a.T * eta * b)[0]

w = q + T * mdot(T, q)
w2 = mdot(w, w)
kappa = sp.sqrt(w2)/v**2
N = w/sp.sqrt(w2)

p1 = sp.Matrix([0, -2*A*k/v, x, y])

# Only q,P-dependent terms are needed for this
# canonical contribution.
psi = -A*kappa**2 - w2*mdot(P, P)/(4*C)
chi = (p1.T*N)[0] + 2*A*kappa/v

point = {
    q[0]: 0,
    q[1]: v**2*k,
    q[2]: 0,
    q[3]: 0,
    P[0]: 0,
    P[1]: 0,
    P[2]: r,
    P[3]: s,
}

def derivative_at(expr, variable):
    reduced = expr.subs({
        key: value
        for key, value in point.items()
        if key != variable
    })
    return sp.diff(reduced, variable).subs(
        variable, point[variable]
    )

print("=== CP098-T22 ===", flush=True)

contributions = []

for i in range(4):
    term = (
        derivative_at(psi, q[i])
        * derivative_at(chi, P[i])
        -
        derivative_at(psi, P[i])
        * derivative_at(chi, q[i])
    )
    term = sp.factor(term)
    contributions.append(term)
    print(f"component {i} = {term}", flush=True)

total = sp.factor(sum(contributions))

print("Full q,P contribution =", total, flush=True)

# Independent reduced expression obtained by
# differentiating the normalized normal vector.
expected = (
    v*k*(x*r + y*s)/(2*C)
)

difference = sp.factor(total - expected)

print("Diagnostic expectation =", expected)
print("Difference =", difference)

if difference == 0:
    print("PASS: transverse contribution verified")
else:
    print("MISMATCH: investigate before proceeding")

print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
