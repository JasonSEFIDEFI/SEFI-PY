"""CP098-T12: Exact canonical Hamiltonian decomposition.

Domain: timelike v>0, kappa>0, C!=0.
No constraints imposed.

Expected: Hc - reconstructed Hc = 0.
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
wn = sp.sqrt(metric(w, w))
kappa = wn/v**2
N = w/wn

P2 = metric(P, P)

Hc = (
    contract(p0, u)
    + contract(p1, q)
    + v*(B + A*kappa**2)
    - v**5*kappa**2*P2/(4*C)
)

chiT = contract(p1, T)
chiN = contract(p1, N) + 2*A*kappa/v

psiT = (
    contract(p0, T)
    + B
    - A*kappa**2
    - v**4*kappa**2*P2/(4*C)
)

reconstructed = (
    v*psiT
    + v**2*kappa*chiN
    - Tq*chiT
)

geometric_residual = sp.Matrix([
    sp.simplify(
        q[i] - (v**2*kappa*N[i] - Tq*T[i])
    )
    for i in range(4)
])

assert geometric_residual == sp.zeros(4, 1)

residual = sp.simplify(
    sp.expand(Hc - reconstructed)
)

assert residual == 0

print("CP098-T12: PASS")
print("Geometric decomposition: exact")
print("Hamiltonian decomposition: exact")
print("Constraints imposed: none")
print("Residual:", residual)
print("SymPy:", sp.__version__)
