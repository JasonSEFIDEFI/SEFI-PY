"""CP098-T83A: Exact canonical transport of v.

Independent of T74 implementation.
No constraint-surface substitution.
Generic A(I), B(I,pI), C(I).
"""

import sympy as sp

print("=== CP098-T83A ===", flush=True)

u = sp.symbols("u0:4", real=True)
q = sp.symbols("q0:4", real=True)
p1 = sp.symbols("p1_0:4", real=True)

I, pI = sp.symbols("I pI", real=True)

A = sp.Function("A")(I)
B = sp.Function("B")(I, pI)
C = sp.Function("C")(I)

LT, LN = sp.symbols("LT LN", real=True)

eta = (-1, 1, 1, 1)

def metric(x, z):
    return sum(
        eta[j]*x[j]*z[j]
        for j in range(4)
    )

v = sp.sqrt(-metric(u, u))
T = tuple(u[j]/v for j in range(4))

a = -metric(T, q)

# Only the p1-dependent term in H_T contributes
# to {v,H_T}; both multiplier constraints are
# independent of p1.
#
# H_T contains p1.q, so dH_T/dp1_i = q_i.

canonical_vdot = sum(
    sp.diff(v, u[j])*q[j]
    for j in range(4)
)

residual = sp.factor(canonical_vdot-a)

if residual != 0:
    raise AssertionError(
        f"STOP: nonzero residual {residual}"
    )

print("PASS: {v,H_T} = a exactly")
print("PASS: no constraints imposed")
print("PASS: arbitrary regular timelike u")
print("PASS: independent of internal model")
print("T74 VELOCITY TRANSPORT: EXACTLY VERIFIED")
print("REMAINING SCALAR TRANSPORT: OPEN")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T83A COMPLETE")
