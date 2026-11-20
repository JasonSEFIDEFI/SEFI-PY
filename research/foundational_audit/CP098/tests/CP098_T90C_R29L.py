import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29L ===", flush=True)

rho, z, b = S.symbols("rho z b", real=True)

# Numerator of D from R29K.
ND = (
    3694084*b**2
    + 7626496*b*rho*z
    + 802816*rho**4
    + 41807344*rho**2
    + 3813248*z**2
    + 551342037
)

# Polynomial numerator of Ddot from R29K,
# omitting its nonzero overall rational factor.
NT = (
    25858588*b**3*rho**2
    + 133463680*b**2*rho**3*z
    + 2774257084*b**2*rho*z
    + 5619712*b*rho**6
    + 305997776*b*rho**4
    + 106770944*b*rho**2*z**2
    + 4202944071*b*rho**2
    + 1374199248*b*z**2
    + 10841600*rho**5*z
    + 490463648*rho**3*z
    + 26692736*rho*z**3
    + 5491255866*rho*z
)

print("STAGE 1: exact polynomial setup", flush=True)

assert S.Poly(ND, b).degree() == 2
assert S.Poly(NT, b).degree() == 3

D = ND / S.Integer(14776336)
Ddot = -2*NT / (
    S.Integer(923521)*(112*rho**2+2883)
)

assert S.simplify(D*14776336-ND) == 0
assert S.simplify(
    Ddot*923521*(112*rho**2+2883)+2*NT
) == 0

print("PASS: R29K polynomial transcription", flush=True)

print("STAGE 2: restrict rho=2", flush=True)

A = S.Poly(ND.subs(rho, 2), b, z)
C = S.Poly(NT.subs(rho, 2), b, z)

print("D polynomial degree:", A.total_degree(), flush=True)
print("Ddot polynomial degree:", C.total_degree(), flush=True)

print("STAGE 3: eliminate b", flush=True)

R = S.factor(
    S.resultant(A.as_expr(), C.as_expr(), b)
)

assert R != 0

print("PASS: resultant nonzero", flush=True)
print("RESULTANT FACTORIZATION =", R, flush=True)

print("STAGE 4: real candidate roots", flush=True)

poly = S.Poly(R, z)
print("RESULTANT DEGREE =", poly.degree(), flush=True)

real_count = poly.count_roots(-S.oo, S.oo)

print("REAL ROOT COUNT =", real_count, flush=True)

if real_count == 0:
    print("PASS: no simultaneous real solutions at rho=2")
    print("R29L: RHO=2 SLICE EXCLUDED")
else:
    print("REAL CANDIDATES EXIST")
    print("NEXT: verify candidates in original equations")
    print("NO SIMULTANEOUS SOLUTION CLAIMED")

print("R29L EVALUATION COMPLETE")
print("GENERAL rho FAMILY: OPEN")
print("GENERAL-w EXTENSION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
