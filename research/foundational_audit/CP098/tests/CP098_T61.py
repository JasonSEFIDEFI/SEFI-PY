"""CP098-T61: Internal derivative dependency audit.

Determine what additional information is required
to evolve E, given its dependence on DA={A,B}.

No specific physical coupling A(I) is assumed.
"""

import sympy as sp

print("=== CP098-T61 ===", flush=True)

I1, I2, p1, p2 = sp.symbols(
    "I1 I2 p1 p2", real=True
)
h1, h2 = sp.symbols(
    "h1 h2", positive=True
)
K1, K2, m = sp.symbols(
    "K1 K2 m", real=True
)

A = sp.Function("A")(I1, I2)
C = sp.Function("C")(I1, I2)

B = (
    m
    + p1**2/(2*h1)
    + p2**2/(2*h2)
    + K1*I1**2/2
    + K2*I2**2/2
)

coords = (I1, I2)
momenta = (p1, p2)

def PB(F, G):
    return sp.expand(sum(
        sp.diff(F, I)*sp.diff(G, p)
        - sp.diff(F, p)*sp.diff(G, I)
        for I, p in zip(coords, momenta)
    ))

DA = PB(A, B)

expected_DA = (
    sp.diff(A, I1)*p1/h1
    + sp.diff(A, I2)*p2/h2
)

assert sp.simplify(DA-expected_DA) == 0
print("PASS: DA canonical definition", flush=True)

DDA = sp.expand(PB(DA, B))

expected_DDA = (
    sp.diff(A, I1, 2)*p1**2/h1**2
    + 2*sp.diff(A, I1, I2)*p1*p2/(h1*h2)
    + sp.diff(A, I2, 2)*p2**2/h2**2
    - K1*I1*sp.diff(A, I1)/h1
    - K2*I2*sp.diff(A, I2)/h2
)

assert sp.simplify(DDA-expected_DDA) == 0

print("PASS: {DA,B} derived", flush=True)
print("{DA,B} =", DDA, flush=True)

# Check the contribution from a potential-like
# geometric Hamiltonian term V(I).
V = sp.Function("V")(I1, I2)

DA_V = PB(DA, V)

expected_DA_V = -(
    sp.diff(A, I1)*sp.diff(V, I1)/h1
    + sp.diff(A, I2)*sp.diff(V, I2)/h2
)

assert sp.simplify(DA_V-expected_DA_V) == 0

print("PASS: {DA,V} derived", flush=True)

# Construct two smooth couplings that have the
# same A and DA at one phase-space point but
# different Hessians.

A0, alpha, beta, pstar = sp.symbols(
    "A0 alpha beta pstar",
    real=True, nonzero=True
)

A_linear = A0 + alpha*I1
A_curved = A0 + alpha*I1 + beta*I1**2/2

point = {
    I1: 0,
    I2: 0,
    p1: pstar,
    p2: 0
}

values = [
    sp.simplify(F.subs(point))
    for F in (A_linear, A_curved)
]
assert values[0] == values[1]

DA_values = [
    sp.simplify(PB(F, B).subs(point))
    for F in (A_linear, A_curved)
]
assert sp.simplify(DA_values[0]-DA_values[1]) == 0

DDA_values = [
    sp.simplify(PB(PB(F, B), B).subs(point))
    for F in (A_linear, A_curved)
]

difference = sp.factor(
    DDA_values[1]-DDA_values[0]
)

assert sp.simplify(
    difference-beta*pstar**2/h1**2
) == 0

print("PASS: equal A and DA at comparison point",
      flush=True)
print("Difference in {DA,B} =", difference,
      flush=True)
print("PASS: Hessian dependence demonstrated",
      flush=True)

# C(I) also has nontrivial canonical evolution.
DC = PB(C, B)

expected_DC = (
    sp.diff(C, I1)*p1/h1
    + sp.diff(C, I2)*p2/h2
)

assert sp.simplify(DC-expected_DC) == 0

print("PASS: C evolution depends on its gradient",
      flush=True)

print("INTERNAL DERIVATIVE DEPENDENCIES:")
print("A gradient: required for DA")
print("A Hessian: required for {DA,B}")
print("C gradient: required for {C,B}")
print("Additional geometric terms: must be included")
print("FULL E EVOLUTION: NOT YET DERIVED")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
