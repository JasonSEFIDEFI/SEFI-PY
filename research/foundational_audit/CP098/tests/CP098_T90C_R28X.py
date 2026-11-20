
import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28X ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R27U.py")

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

B = ns["B"]
phi1 = ns["phi1"]

w, rho, z, t, a, b = S.symbols(
    "w rho z t a b", real=True
)

IS_value = 12/w**2 + S.Rational(448,961)*rho**2

M_value = 1 + (
    S.Rational(9,4)
    + IS_value**2
    + 1 + a**2 + b**2
)/2

point = {}

for sym, val in zip(u, [2,0,0,0]):
    point[sym] = S.Integer(val)

for sym, val in zip(q, [0,w,0,0]):
    point[sym] = S.sympify(val)

for sym, val in zip(P, [0,0,rho,0]):
    point[sym] = S.sympify(val)

for sym, val in zip(p1, [0,0,z,t]):
    point[sym] = S.sympify(val)

p0_values = [
    2*w**2*rho**2/31-M_value,
    -8*w*rho*z/31,
    0,
    0
]

for sym, val in zip(p0, p0_values):
    point[sym] = S.sympify(val)

point.update({
    ns["IA"]: S.Rational(-3,2),
    ns["IS"]: IS_value,
    ns["IW"]: S.Integer(1),
    ns["pA"]: S.Integer(0),
    ns["pS"]: a,
    ns["pW"]: b
})

def at(expr):
    return S.factor(S.cancel(expr.subs(point)))

print("STAGE 1: full canonical derivatives", flush=True)

terms = []

for i in range(4):
    term = (
        u[i]*S.diff(B,q[i])
        - P[i]*S.diff(B,p1[i])
    )
    terms.append(term)
    print("Derivative", i+1, "of 4", flush=True)

print("STAGE 2: family substitution", flush=True)

result = S.Integer(0)

for i, term in enumerate(terms,1):
    result += at(term)
    result = S.factor(result)
    print("Substituted", i, "of 4", flush=True)

print("{B,phi1} on family =", result, flush=True)

if result == 0:
    print("PASS: family-wide {B,phi1}=0")
else:
    print("NONZERO FAMILY EXPRESSION")
    print("Factorization =", S.factor(result))

print("STAGE 3: primary regression", flush=True)

assert S.simplify(
    phi1.subs(point)
) == 0

print("PASS: phi1=0 throughout family")

print("R28X COMPLETE")
print("K_D1 FAMILY IDENTITY: STILL REQUIRES JACOBI ANALYSIS")
print("Ddot FAMILY CONDITION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3),
    flush=True
)
