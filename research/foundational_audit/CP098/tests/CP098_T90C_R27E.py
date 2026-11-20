import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27E ===", flush=True)

v, W, C = S.symbols("v W C", positive=True)
d, A, M, J = S.symbols("d A M J", real=True)
e, f, g, P2 = S.symbols("e f g P2", real=True)
phi1, phi2 = S.symbols("phi1 phi2", real=True)

s = v**2

B = (
    -e
    -v*W*f/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(phi2+d*phi1/s)
    +2*W*J/v**2
)

# Canonical generator transformations.
# These hold for the scalar contractions
# in any spacetime dimension.

D1 = {
    v:0,
    W:0,
    d:-s,
    e:g,
    f:-P2,
    g:0,
    P2:0,
    phi1:0,
    phi2:phi1,
}

D2 = {
    v:0,
    W:2*W,
    d:d,
    e:e,
    f:-f,
    g:0,
    P2:-2*P2,
    phi1:-phi1,
    phi2:0,
}

def action(expr, rules):
    return sum(
        S.diff(expr, x)*dx
        for x, dx in rules.items()
    )

def check(name, expr):
    result = S.factor(expr)
    if result != 0:
        raise AssertionError(f"{name}: {result}")
    print("PASS:", name, flush=True)

print("STAGE 1: derive exact coefficients", flush=True)

K1 = S.factor(action(B, D1))
K2 = S.factor(action(B, D2))

print("K1 =", K1, flush=True)
print("K2 =", K2, flush=True)
print("STAGE 1 PASS", flush=True)

print("STAGE 2: independent compact formulas", flush=True)

expected_K1 = (
    -g
    +v*W*P2/(2*C)
    -s*W*P2/(4*C*v)
    -5*A*s*W/v**5
    -s*M/v
)

# Under phi2, the A-dependent term involving
# phi2 and d*phi1/s scales with weight 2.
# Its complete contribution is retained.

expected_K2 = (
    -e
    -v*W*f/(2*C)
    +d*W*P2/(4*C*v)
    +15*A*d*W/v**5
    +d*M/v
    -4*A*W/(C*v**2)*(phi2+d*phi1/s)
    +4*W*J/v**2
)

check("K1 exact identity", K1-expected_K1)
check("K2 exact identity", K2-expected_K2)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: singular-stratum reduction", flush=True)

singular_K1 = S.factor(K1.subs(A, 0))
singular_K2 = S.factor(K2.subs(A, 0))

print("K1 at A=0 =", singular_K1, flush=True)
print("K2 at A=0 =", singular_K2, flush=True)

# B=0 is an additional condition.
# Do not set unrelated phase-space quantities
# to zero without justification.

print("STAGE 3 PASS", flush=True)

print()
print("R27E COMPLETE")
print("MULTIPLIER COEFFICIENTS: DERIVED")
print("SINGULAR A=0 REDUCTION: DERIVED")
print("B=0 CONSISTENCY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3)
)
