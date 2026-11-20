import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29M ===", flush=True)

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


print("STAGE 2: eliminate b with rho variable", flush=True)

A = S.Poly(ND, b)
C = S.Poly(NT, b)

R = S.resultant(A.as_expr(), C.as_expr(), b)
assert R != 0

print("PASS: nonzero symbolic resultant", flush=True)

print("STAGE 3: polynomial structure", flush=True)

poly = S.Poly(R, rho, z)

print("TOTAL DEGREE =", poly.total_degree(), flush=True)
print("NUMBER OF TERMS =", len(poly.terms()), flush=True)

odd_rho = [
    monom for monom, coeff in poly.terms()
    if monom[0] % 2
]
odd_z = [
    monom for monom, coeff in poly.terms()
    if monom[1] % 2
]

print("ODD RHO POWERS =", len(odd_rho), flush=True)
print("ODD Z POWERS =", len(odd_z), flush=True)

print("STAGE 4: rho=2 regression", flush=True)

R2 = S.Poly(R.subs(rho, 2), z)
assert R2.degree() == 6
assert R2.count_roots(-S.oo, S.oo) == 0

print("PASS: R29L rho=2 regression", flush=True)

print("STAGE 5: coefficient analysis", flush=True)

if not odd_rho and not odd_z:
    x, y = S.symbols("x y", nonnegative=True)

    transformed = S.Add(*[
        coeff*x**(monom[0]//2)*y**(monom[1]//2)
        for monom, coeff in poly.terms()
    ])

    xy = S.Poly(transformed, x, y)

    print("POLYNOMIAL IN x=rho^2, y=z^2", flush=True)
    print("DEGREES =", xy.degree(x), xy.degree(y), flush=True)

    coefficients = xy.coeffs()

    positive = sum(
        1 for c in coefficients if c.is_positive is True
    )
    negative = sum(
        1 for c in coefficients if c.is_negative is True
    )
    zero = sum(
        1 for c in coefficients if c.is_zero is True
    )

    print("POSITIVE COEFFICIENTS =", positive, flush=True)
    print("NEGATIVE COEFFICIENTS =", negative, flush=True)
    print("ZERO COEFFICIENTS =", zero, flush=True)

    assert positive + negative + zero == len(coefficients)

    if positive == len(coefficients):
        print("CERTIFIED: resultant positive for x>0,y>=0",
              flush=True)
    else:
        print("OPEN: coefficient signs do not establish positivity",
              flush=True)
else:
    print("OPEN: resultant not even in both variables",
          flush=True)

print("R29M COMPUTATION COMPLETE")
print("GENERAL REAL SOLUTION CLASSIFICATION: OPEN")
print("GENERAL-w EXTENSION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-start, 3), flush=True)
