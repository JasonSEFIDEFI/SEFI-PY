"""CP098-T62: Complete internal evolution of DA.

Compare direct canonical differentiation against
independent decomposition of {DA,H_int}.

A(I) and C(I) remain arbitrary smooth functions.
Two internal pairs provide a symbolic test case.
"""

import sympy as sp

print("=== CP098-T62 ===", flush=True)

I1, I2, p1, p2 = sp.symbols(
    "I1 I2 p1 p2", real=True
)
h1, h2 = sp.symbols("h1 h2", positive=True)
K1, K2, m = sp.symbols(
    "K1 K2 m", real=True
)
v, k = sp.symbols("v k", positive=True)
P2 = sp.symbols("P_squared", real=True)

A = sp.Function("A")(I1, I2)
C = sp.Function("C")(I1, I2)

coords = (I1, I2)
momenta = (p1, p2)
hs = (h1, h2)

def PB(F, G):
    return sp.expand(sum(
        sp.diff(F, I)*sp.diff(G, p)
        - sp.diff(F, p)*sp.diff(G, I)
        for I, p in zip(coords, momenta)
    ))

def zero(expr):
    return sp.cancel(sp.together(expr)) == 0

B = (
    m
    + p1**2/(2*h1)
    + p2**2/(2*h2)
    + K1*I1**2/2
    + K2*I2**2/2
)

DA = PB(A, B)

H_int = (
    v*B
    + v*A*k**2
    - v**5*k**2*P2/(4*C)
)

# Route 1: direct canonical differentiation.
direct = PB(DA, H_int)

print("PASS: direct internal bracket constructed",
      flush=True)

# Route 2: decompose by Hamiltonian terms.
DDA = PB(DA, B)

DA_A = -sum(
    sp.diff(A, I)**2/h
    for I, h in zip(coords, hs)
)

DA_invC = sum(
    sp.diff(A, I)*sp.diff(C, I)/(h*C**2)
    for I, h in zip(coords, hs)
)

assert zero(PB(DA, A)-DA_A)
assert zero(PB(DA, 1/C)-DA_invC)

print("PASS: independent coupling brackets",
      flush=True)

decomposed = (
    v*DDA
    + v*k**2*DA_A
    - v**5*k**2*P2*DA_invC/4
)

assert zero(direct-decomposed)

print("PASS: direct/decomposed agreement",
      flush=True)

# Verify that the A and C contributions
# are not silently omitted.
assert zero(
    decomposed.subs(P2, 0)
    - (v*DDA+v*k**2*DA_A)
)

print("PASS: A coupling retained", flush=True)

P2_coefficient = sp.diff(decomposed, P2)

assert zero(
    P2_coefficient
    + v**5*k**2*DA_invC/4
)

print("PASS: C coupling retained", flush=True)

# Cross-check with concrete smooth polynomial
# couplings, independent of symbolic derivatives.
alpha, beta, gamma, c0, c1, c2 = sp.symbols(
    "alpha beta gamma c0 c1 c2",
    real=True
)

A_poly = (
    alpha*I1
    + beta*I1*I2
    + gamma*I2**2
    + 2
)

C_poly = c0+c1*I1+c2*I2**2

replacement = {
    A:A_poly,
    C:C_poly
}

# Substitute the functions first, then calculate
# their canonical derivatives afresh.
DA_poly = PB(A_poly, B)
H_poly = H_int.subs(replacement)

direct_poly = PB(DA_poly, H_poly)

decomp_poly = decomposed.subs(replacement).doit()

assert zero(direct_poly-decomp_poly)

print("PASS: polynomial-coupling cross-check",
      flush=True)

print("GENERAL INTERNAL RESULT:")
print(
    "{DA,H_int} = v*{DA,B}"
    " + v*k^2*{DA,A}"
    " - v^5*k^2*P^2*{DA,1/C}/4"
)

print("INTERNAL DA EVOLUTION: VERIFIED")
print("EXTERNAL CONTRIBUTIONS TO E_DOT: OPEN")
print("FULL E EVOLUTION: OPEN")
print("FULL GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
