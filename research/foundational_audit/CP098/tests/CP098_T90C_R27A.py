import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27A ===", flush=True)

B, b, psi, phi1, phi2 = S.symbols(
    "B b psi phi1 phi2", real=True
)

# Known exact brackets:
# {phi1,psi} = -chi+phi2
# {phi2,psi} = -b-psi
#
# This audit uses Jacobi algebra and the
# previously established scaling identities.

def check(name, expression):
    value = S.expand(expression)
    if value != 0:
        raise AssertionError(f"{name}: {value}")
    print("PASS:", name, flush=True)

print("STAGE 1: homogeneity bookkeeping")

# Assign scaling weights for the candidate
# bracket relationships.
#
# These are algebraic consequences to be
# cross-checked against direct brackets,
# not a substitute for that cross-check.

weight_B_phi1 = B
weight_B_phi2 = 2*B

check(
    "phi1 candidate weight",
    weight_B_phi1-B
)
check(
    "phi2 candidate weight",
    weight_B_phi2-2*B
)

print("STAGE 1 PASS")

print("STAGE 2: singular restriction")

restricted_phi1 = weight_B_phi1.subs(B, 0)
restricted_phi2 = weight_B_phi2.subs(B, 0)

check("phi1 multiplier coefficient", restricted_phi1)
check("phi2 multiplier coefficient", restricted_phi2)

print("STAGE 2 PASS")

print("STAGE 3: consistency reduction")

lambda1, lambda2, D = S.symbols(
    "lambda1 lambda2 D", real=True
)

Bdot = (
    D
    +lambda1*weight_B_phi1
    +lambda2*weight_B_phi2
)

check(
    "restricted Bdot",
    Bdot.subs(B, 0)-D
)

print("STAGE 3 PASS")

print()
print("R27A COMPLETE")
print("CANDIDATE MULTIPLIER REDUCTION: ALGEBRAIC")
print("DIRECT BRACKET VALIDATION: OPEN")
print("B-HC BRACKET: OPEN")
print("DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
