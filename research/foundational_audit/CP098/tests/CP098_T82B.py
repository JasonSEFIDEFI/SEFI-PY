"""CP098-T82B: Exact CP098 constraint construction audit.

No constraint substitution.
No numerical finite differences.
No claims of physical validation.
"""

import sympy as sp

print("=== CP098-T82B ===", flush=True)

# Canonical coordinates and momenta.
y = sp.symbols("y0:4", real=True)
u = sp.symbols("u0:4", real=True)
q = sp.symbols("q0:4", real=True)

p0 = sp.symbols("p0_0:4", real=True)
p1 = sp.symbols("p1_0:4", real=True)
P = sp.symbols("P0:4", real=True)

I, pI = sp.symbols("I pI", real=True)

# Generic internal functions.
A = sp.Function("A")(I)
C = sp.Function("C")(I)
B = sp.Function("B")(I, pI)

eta = (-1, 1, 1, 1)

def metric(x, z):
    return sum(
        eta[j] * x[j] * z[j]
        for j in range(4)
    )

def pairing(p, x):
    return sum(p[j] * x[j] for j in range(4))

v = sp.sqrt(-metric(u, u))
T = tuple(u[j]/v for j in range(4))

Tq = metric(T, q)
w = tuple(q[j] + T[j]*Tq for j in range(4))

w2 = metric(w, w)
k = sp.sqrt(w2)/v**2
N = tuple(w[j]/sp.sqrt(w2) for j in range(4))

P2 = metric(P, P)

phiT = pairing(P, T)
phiN = pairing(P, N)

chiT = pairing(p1, T)
chiN = pairing(p1, N) + 2*A*k/v

psiT = (
    pairing(p0, T)
    + B
    - A*k**2
    - v**4*k**2*P2/(4*C)
)

Hc = (
    pairing(p0, u)
    + pairing(p1, q)
    + v*(B+A*k**2)
    - v**5*k**2*P2/(4*C)
)

LT, LN = sp.symbols("LT LN", real=True)
HT = Hc + LT*phiT + LN*phiN

def zero(expr):
    return sp.cancel(sp.together(expr)) == 0

# Exact geometric identities.
assert zero(metric(T, T)+1)
print("PASS: T.T = -1", flush=True)

assert zero(metric(T, N))
print("PASS: T.N = 0", flush=True)

assert zero(metric(N, N)-1)
print("PASS: N.N = 1", flush=True)

assert zero(pairing(p1, q) -
            (v**2*k*pairing(p1, N)-Tq*chiT))
print("PASS: q decomposition", flush=True)

# Hamiltonian decomposition on the entire
# regular domain, not merely on constraints.
decomposition = (
    v*psiT
    + v**2*k*chiN
    - Tq*chiT
)

assert zero(Hc-decomposition)
print("PASS: exact Hamiltonian decomposition", flush=True)

# Check full canonical pair inventory.
pairs = (
    list(zip(y, p0))
    + list(zip(u, p1))
    + list(zip(q, P))
    + [(I, pI)]
)

assert len(pairs) == 13
print("PASS: 13 canonical pairs", flush=True)

# The model is translation invariant:
# y does not occur in its Hamiltonian
# or five constraints.
for expr in (Hc, phiT, phiN, chiT, chiN, psiT):
    assert not any(expr.has(component) for component in y)

print("PASS: translation independence", flush=True)

print("CONSTRAINT EXPRESSIONS: CONSTRUCTED")
print("ACTUAL CONSTRAINT JACOBI IDENTITY: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("GAUGE GENERATOR: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T82B COMPLETE")
