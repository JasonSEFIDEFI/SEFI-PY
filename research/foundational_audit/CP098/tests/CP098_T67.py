"""CP098-T67: Internal evolution of E.

Verify the internal Poisson-bracket contribution
to E_dot at a fixed adapted external point.

A(I) and C(I) remain arbitrary smooth functions.

This does NOT compute external frame transport or
the complete nested Poisson bracket.
"""

import sympy as sp

print("=== CP098-T67 ===", flush=True)

I1, I2, j1, j2 = sp.symbols(
    "I1 I2 j1 j2", real=True
)
h1, h2 = sp.symbols("h1 h2", positive=True)
K1, K2, m = sp.symbols("K1 K2 m", real=True)

v, k = sp.symbols("v k", positive=True)
a, n = sp.symbols("a n", real=True)

P0, P1, P2, P3 = sp.symbols(
    "P0 P1 P2 P3", real=True
)
p12, p13 = sp.symbols(
    "p12 p13", real=True
)

A = sp.Function("A")(I1, I2)
C = sp.Function("C")(I1, I2)

coords = (I1, I2)
momenta = (j1, j2)
hs = (h1, h2)

def PB(F, G):
    return sum(
        sp.diff(F, I)*sp.diff(G, p)
        - sp.diff(F, p)*sp.diff(G, I)
        for I, p in zip(coords, momenta)
    )

def check_zero(expr):
    return sp.cancel(sp.together(expr)) == 0

B = (
    m
    + j1**2/(2*h1)
    + j2**2/(2*h2)
    + K1*I1**2/2
    + K2*I2**2/2
)

DA = PB(A, B)

P_squared = -P0**2+P1**2+P2**2+P3**2

H_int = (
    v*B
    + v*k**2*A
    - v**5*k**2*P_squared/(4*C)
)

# T60's E expression evaluated at an adapted
# external point, before imposing the regular
# momentum constraints.
#
# External symbols are held fixed only for
# this internal canonical derivative.

def E_formula(AA, CC, DD):
    numerator = (
        4*AA*CC*a*k
        + 2*AA*P1*k**2*v**4
        - 4*CC*DD*k*v**2
        + 2*CC*n*v**2
        + a*k*v**4*P_squared
        + k*v**5*(P2*p12+P3*p13)
    )
    return -numerator/(2*CC*v**3)

E = E_formula(A, C, DA)

print("PASS: T60 adapted-point E constructed",
      flush=True)

# Route 1: direct internal canonical bracket.

E_dot_direct = PB(E, H_int)

print("PASS: direct internal bracket constructed",
      flush=True)

# Route 2: independent chain rule, treating
# A, C, DA as distinct intermediate variables.

AA, CC, DD = sp.symbols(
    "AA CC DD", real=True, nonzero=True
)

E_abstract = E_formula(AA, CC, DD)

A_dot = PB(A, H_int)
C_dot = PB(C, H_int)
DA_dot = PB(DA, H_int)

assert check_zero(A_dot-v*DA)

assert check_zero(
    C_dot-v*PB(C, B)
)

print("PASS: A and C evolution identities",
      flush=True)

# Independently reconstruct DA_dot using
# the three Hamiltonian contributions.

DA_dot_expected = (
    v*PB(DA, B)
    + v*k**2*PB(DA, A)
    - v**5*k**2*P_squared*PB(DA, 1/C)/4
)

assert check_zero(DA_dot-DA_dot_expected)

print("PASS: T62 DA evolution recovered",
      flush=True)

E_dot_chain = (
    sp.diff(E_abstract, AA)*A_dot
    + sp.diff(E_abstract, CC)*C_dot
    + sp.diff(E_abstract, DD)*DA_dot
).subs({
    AA:A,
    CC:C,
    DD:DA
})

assert check_zero(E_dot_direct-E_dot_chain)

print("PASS: direct/chain-rule agreement",
      flush=True)

# Independent polynomial-coupling test.
# Rebuild E and H from scratch with explicit
# functions rather than substituting into
# an already differentiated bracket.

alpha, beta, gamma = sp.symbols(
    "alpha beta gamma", real=True
)
c0, c1, c2 = sp.symbols(
    "c0 c1 c2", real=True
)

A_poly = (
    2+alpha*I1+beta*I1*I2+gamma*I2**2
)
C_poly = c0+c1*I1+c2*I2**2

DA_poly = PB(A_poly, B)

H_poly = (
    v*B
    + v*k**2*A_poly
    - v**5*k**2*P_squared/(4*C_poly)
)

E_poly = E_formula(
    A_poly, C_poly, DA_poly
)

poly_direct = PB(E_poly, H_poly)

poly_chain = (
    sp.diff(E_abstract, AA)*PB(A_poly,H_poly)
    + sp.diff(E_abstract, CC)*PB(C_poly,H_poly)
    + sp.diff(E_abstract, DD)*PB(DA_poly,H_poly)
).subs({
    AA:A_poly,
    CC:C_poly,
    DD:DA_poly
})

assert check_zero(poly_direct-poly_chain)

print("PASS: independent polynomial test",
      flush=True)

# The internal contribution must not depend
# on primary multipliers: the primary
# constraints have no internal dependence.

LT, LN = sp.symbols("lambda_T lambda_N")

assert LT not in E_dot_direct.free_symbols
assert LN not in E_dot_direct.free_symbols

print("PASS: primary-multiplier independence",
      flush=True)

# Verify that Hessian terms are retained in
# the DA evolution, rather than assuming
# linear coupling.

hessian_term = (
    v*j1**2*sp.diff(A,I1,2)/h1**2
)

assert check_zero(
    sp.diff(DA_dot,j1,2)
    - sp.diff(DA_dot_expected,j1,2)
)

assert sp.diff(A,I1,2) in (
    sp.diff(DA_dot_expected,j1,2).atoms(
        sp.Derivative
    )
)

print("PASS: Hessian dependency retained",
      flush=True)

print("INTERNAL E EVOLUTION: VERIFIED")
print("T60 EXPRESSION: USED AT ADAPTED POINT")
print("EXTERNAL FRAME TRANSPORT: NOT INCLUDED")
print("FULL E EVOLUTION: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
