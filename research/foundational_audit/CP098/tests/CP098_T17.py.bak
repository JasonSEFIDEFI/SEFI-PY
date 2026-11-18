"""CP098-T17: Direct tangent-normal secondary bracket.

Exact regular adapted-point calculation.
All canonical variables remain independent until
Poisson-bracket differentiation is completed.

A pointwise pass is not a general closure proof.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

u = sp.Matrix(sp.symbols("u0:4", real=True))
q = sp.Matrix(sp.symbols("q0:4", real=True))
P = sp.Matrix(sp.symbols("P0:4", real=True))
p1 = sp.Matrix(sp.symbols("r0:4", real=True))
p0 = sp.Matrix(sp.symbols("s0:4", real=True))

I, ps = sp.symbols("I ps", real=True)

def mdot(x, y):
    return (x.T * eta * y)[0]

def contract(p, x):
    return (p.T * x)[0]

v = sp.sqrt(-mdot(u, u))
T = u/v
Tq = mdot(T, q)
w = q + T*Tq
w2 = mdot(w, w)
k2 = w2/v**4
k = sp.sqrt(k2)
N = w/sp.sqrt(w2)

# Explicit internal specialization:
# A=I, C=1, B=ps^2/2.
A = I
C = sp.Integer(1)
B = ps**2/2

P2 = mdot(P, P)

psi = (
    contract(p0, T)
    + B
    - A*k2
    - v**4*k2*P2/(4*C)
)

chi = contract(p1, N) + 2*A*k/v

point = {}

def assign(vec, values):
    point.update(zip(list(vec), map(sp.sympify, values)))

assign(u, [3, 0, 0, 0])
assign(q, [1, 18, 0, 0])
assign(P, [0, 0, 1, 0])
assign(p1, [0, -sp.Rational(4, 3), 0, 0])
assign(p0, [sp.Rational(169, 2), 2, 0, 0])

point[I] = sp.Integer(1)
point[ps] = sp.Integer(1)

def at(expr):
    return sp.simplify(expr.subs(point))

def canonical_pair(f, g, xs, momenta):
    total = sp.Integer(0)
    for x, p in zip(xs, momenta):
        total += (
            at(sp.diff(f, x))*at(sp.diff(g, p))
            - at(sp.diff(f, p))*at(sp.diff(g, x))
        )
    return sp.simplify(total)

print("=== CP098-T17 ===", flush=True)

# Constraint values at the chosen point.
phiT = contract(P, T)
phiN = contract(P, N)
chiT = contract(p1, T)

print("phi_T =", at(phiT))
print("phi_N =", at(phiN))
print("chi_T =", at(chiT))
print("chi_N =", at(chi))
print("psi_T =", at(psi))

assert all(
    at(expr) == 0
    for expr in (phiT, phiN, chiT, chi, psi)
)

print("PASS: five established constraints satisfied")

# Full canonical bracket, including internal variables.
part_u = canonical_pair(psi, chi, u, p1)
part_q = canonical_pair(psi, chi, q, P)
part_I = canonical_pair(psi, chi, [I], [ps])

actual = sp.simplify(part_u + part_q + part_I)

# Previously proposed normal-preservation expression.
D_A = ps

F_N = (
    -contract(p0, N)
    + 2*k*D_A
    + 6*A*k*Tq/v**2
    - v**3*k*mdot(p1, P)/(2*C)
    + v**2*k*Tq*P2/(2*C)
)

required = (
    4*A*Tq*k/v**3
    - F_N/v
)

print("u,p1 contribution =", part_u)
print("q,P contribution =", part_q)
print("internal contribution =", part_I)
print("Direct bracket =", actual)
print("T16 required bracket =", at(required))

difference = sp.simplify(actual - at(required))

print("Difference =", difference)

if difference == 0:
    print("PASS: T16 condition at this point")
else:
    print("MISMATCH: T16 condition fails at this point")
    print("STOP: investigate missing terms or assumptions")

print("STATUS: general closure not established")
