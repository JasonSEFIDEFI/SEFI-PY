"""CP098-T83G: Exact internal canonical transport.

Independent of T74 implementation.
Arbitrary smooth internal functions A(I), C(I), B(I,pI).
No numerical specialization or constraint substitution.
"""

import sympy as sp

print("=== CP098-T83G ===", flush=True)

I, pI = sp.symbols("I pI", real=True)
v, k, R = sp.symbols("v k R", real=True)

A = sp.Function("A")(I)
C = sp.Function("C")(I)
B = sp.Function("B")(I, pI)

def pb(f, g):
    return (
        sp.diff(f, I)*sp.diff(g, pI)
        -sp.diff(f, pI)*sp.diff(g, I)
    )

def equal(f, g):
    return sp.simplify(sp.expand(f-g)) == 0

# Internal-dependent part of the full Hamiltonian.
# External multipliers multiply constraints that
# have no internal dependence in the CP098 model.

H_internal = (
    v*B
    +v*A*k**2
    -v**5*k**2*R/(4*C)
)

DA = pb(A, B)
DC = pb(C, B)

H = pb(DA, B)
J = pb(DA, A)
K = pb(DA, 1/C)

print("PASS: internal Poisson structure constructed",
      flush=True)

# Independently calculate the canonical rates.

Adot = pb(A, H_internal)
Cdot = pb(C, H_internal)
DAdot = pb(DA, H_internal)

print("PASS: all three canonical rates constructed",
      flush=True)

expected_Adot = v*DA
expected_Cdot = v*DC

expected_DAdot = (
    v*H
    +v*k**2*J
    -v**5*k**2*R*K/4
)

checks = [
    ("A", Adot, expected_Adot),
    ("C", Cdot, expected_Cdot),
    ("DA", DAdot, expected_DAdot),
]

for name, actual, expected in checks:
    if not equal(actual, expected):
        print(
            f"STOP: {name} transport residual =",
            sp.simplify(actual-expected),
            flush=True
        )
        raise AssertionError(
            f"Internal transport mismatch: {name}"
        )

    print(
        f"PASS: exact canonical {name} transport",
        flush=True
    )

# Verify that A and C have no mutual internal
# Poisson bracket, as required by A=A(I), C=C(I).

if not equal(pb(A, C), 0):
    raise AssertionError(
        "STOP: unexpected internal A-C bracket"
    )

print("PASS: {A,C} = 0", flush=True)

# Confirm the coefficient of each independent
# Hamiltonian contribution in DAdot.

if not equal(
    pb(DA, H_internal),
    v*pb(DA, B)
    +v*k**2*pb(DA, A)
    -v**5*k**2*R*pb(DA, 1/C)/4
):
    raise AssertionError(
        "STOP: internal Hamiltonian decomposition"
    )

print("PASS: higher-bracket decomposition", flush=True)

print("PASS: generic internal functions retained")
print("PASS: no constraint surface imposed")
print("PASS: no numerical approximation")
print("INTERNAL TRANSPORT: EXACTLY VERIFIED")
print("GENERAL E_DOT CANONICAL PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T83G COMPLETE")
