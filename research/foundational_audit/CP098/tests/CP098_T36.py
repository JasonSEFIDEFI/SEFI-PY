"""CP098-T36: Full regular tangential consistency.

Independently derives the internal {psi_T,H_c} bracket
for arbitrary differentiable A(I), C(I), quadratic B.

Combines with the external brackets established in T35.
Tests cancellation without imposing it.
"""

import sympy as sp

I1, I2, p1, p2 = sp.symbols(
    "I1 I2 p1 p2", real=True
)
h1, h2 = sp.symbols("h1 h2", positive=True)
k1, k2, m = sp.symbols("k1 k2 m", real=True)

v, kap = sp.symbols("v kap", positive=True)
P2 = sp.symbols("P2", real=True)

A = sp.Function("A")(I1, I2)
C = sp.Function("C")(I1, I2)

B = (
    m
    + p1**2/(2*h1)
    + p2**2/(2*h2)
    + k1*I1**2/2
    + k2*I2**2/2
)

K = v**4 * kap**2 * P2 / 4

psi = B - A*kap**2 - K/C
H = v*B + v*A*kap**2 - v*K/C

pairs = [(I1, p1), (I2, p2)]

def PB(f, g):
    return sp.expand(sum(
        sp.diff(f, I)*sp.diff(g, p)
        - sp.diff(f, p)*sp.diff(g, I)
        for I, p in pairs
    ))

def check(label, expression):
    result = sp.simplify(expression)
    print(label, "=", result, flush=True)
    assert result == 0, label

print("=== CP098-T36 ===", flush=True)

direct_internal = PB(psi, H)

DA = (
    sp.diff(A, I1)*p1/h1
    + sp.diff(A, I2)*p2/h2
)

expected_internal = -2*v*kap**2*DA

print(
    "Direct internal bracket =",
    sp.simplify(direct_internal),
    flush=True
)

check(
    "Internal bracket residual",
    direct_internal - expected_internal
)

# T35 external result:
# {psi,Hc}_ext = -v*kap*FN_ext.
# Here FN_ext is an independent placeholder.
FN_ext = sp.symbols("FN_ext", real=True)

bracket_external = -v*kap*FN_ext

# T25 and T33:
FN_internal = 2*kap*DA
FN_full = FN_ext + FN_internal

# Normal preservation fixes lambda_N when A != 0.
lambda_N = -v**3*FN_full/(2*A)

# T35:
bracket_phi_T = sp.Integer(0)
bracket_phi_N = -2*A*kap/v**2

lambda_T = sp.symbols("lambda_T", real=True)

full_residual = (
    bracket_external
    + direct_internal
    + lambda_T*bracket_phi_T
    + lambda_N*bracket_phi_N
)

check("Full tangential residual", full_residual)

print("PASS: internal canonical bracket")
print("PASS: combined tangential cancellation")
print("REGULAR TANGENTIAL CONSISTENCY: VERIFIED")
print("OTHER CONSTRAINT PRESERVATION: NOT YET COMPLETE")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
print("DYNAMICAL CONSISTENCY: OPEN")
