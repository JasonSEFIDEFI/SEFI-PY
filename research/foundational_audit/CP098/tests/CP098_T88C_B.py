"""CP098-T88C-B: internal canonical derivative audit.

Exact polynomial realization of A(I), C(I), B(I,pI).
Reference values reproduce T88C-A.
No physical interpretation assumed.
"""
import sympy as sp

print("=== CP098-T88C-B ===", flush=True)

I, pI = sp.symbols("I pI", real=True)

# Reference point.
I0 = sp.S.Zero
pI0 = sp.S.Zero

# Smooth, nonconstant test functions.
A = sp.Rational(3,2) + I + I**2/3
C = sp.Rational(5,3) + I/2 + I**2/5
B = (
    sp.Rational(7,4)
    + sp.Rational(2,5)*pI
    + I*pI/3
    + pI**2/7
    + I**2/11
)

def PB(f,g):
    return sp.expand(
        sp.diff(f,I)*sp.diff(g,pI)
        -sp.diff(f,pI)*sp.diff(g,I)
    )

def at0(f):
    return sp.simplify(
        f.subs({I:I0,pI:pI0})
    )

DA = PB(A,B)
DC = PB(C,B)

H = PB(DA,B)
J = PB(DA,A)
K = PB(DA,1/C)

checks = {
    "A": (at0(A),sp.Rational(3,2)),
    "C": (at0(C),sp.Rational(5,3)),
    "B": (at0(B),sp.Rational(7,4)),
    "DA": (at0(DA),sp.Rational(2,5)),
}

for name,(actual,expected) in checks.items():
    residual = sp.simplify(actual-expected)

    print(
        f"{name}_REFERENCE_RESIDUAL =",
        residual,flush=True
    )

    assert residual == 0, (
        f"STOP: {name} reference mismatch"
    )

# Independently check antisymmetry and
# the derivation rule for the internal PB.
f = A*B
g = C

antisymmetry = sp.simplify(
    PB(A,B)+PB(B,A)
)

leibniz = sp.simplify(
    PB(f,g)-A*PB(B,g)-B*PB(A,g)
)

print(
    "INTERNAL_ANTISYMMETRY_RESIDUAL =",
    antisymmetry,flush=True
)
print(
    "INTERNAL_LEIBNIZ_RESIDUAL =",
    leibniz,flush=True
)

assert antisymmetry == 0
assert leibniz == 0

# Confirm higher internal brackets are
# calculated, not replaced with constants.
for name,expr in (
    ("DC",DC),
    ("H",H),
    ("J",J),
    ("K",K)
):
    value = at0(expr)
    assert value.is_finite

    print(
        f"{name}_REFERENCE_VALUE =",
        value,flush=True
    )

print("PASS: explicit internal canonical realization")
print("PASS: internal first and higher derivatives retained")
print("FULL NESTED JACOBI AUDIT: OPEN")
print("GENERAL INTERNAL-FUNCTION PROOF: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T88C-B COMPLETE")
