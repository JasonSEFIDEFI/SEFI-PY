
import sympy as S
import runpy, contextlib, io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28Y ===", flush=True)

folder = Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(folder / "CP098_T90C_R27U.py"))

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]
B, H = ns["B"], ns["H"]
chi, phi2 = ns["chi"], ns["phi2"]

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [
        (ns["IA"], ns["pA"]),
        (ns["IS"], ns["pS"]),
        (ns["IW"], ns["pW"])
    ]
)

w, rho, z, t, a, b = S.symbols(
    "w rho z t a b", real=True
)

ISv = 12/w**2 + S.Rational(448,961)*rho**2
Mv = 1 + (
    S.Rational(9,4) + ISv**2
    + 1 + a**2 + b**2
)/2

point = {}

for vec, values in (
    (u, [2,0,0,0]),
    (q, [0,w,0,0]),
    (P, [0,0,rho,0]),
    (p1, [0,0,z,t]),
    (p0, [
        2*w**2*rho**2/31-Mv,
        -8*w*rho*z/31,
        0,0
    ])
):
    point.update({
        sym: S.sympify(value)
        for sym, value in zip(vec, values)
    })

point.update({
    ns["IA"]: S.Rational(-3,2),
    ns["IS"]: ISv,
    ns["IW"]: S.Integer(1),
    ns["pA"]: S.Integer(0),
    ns["pS"]: a,
    ns["pW"]: b
})

def at(expr):
    return S.factor(S.cancel(expr.subs(point)))

def bracket_at(f, g):
    result = S.Integer(0)
    for i, (x, p) in enumerate(pairs, 1):
        result += (
            at(S.diff(f,x))*at(S.diff(g,p))
            - at(S.diff(f,p))*at(S.diff(g,x))
        )
        print("pair", i, "of 11", flush=True)
    return S.factor(result)

print("STAGE 1: construct unrestricted inner brackets",
      flush=True)

Bphi1 = sum(
    u[i]*S.diff(B,q[i])
    - P[i]*S.diff(B,p1[i])
    for i in range(4)
)

G = chi - phi2

print("PASS: inner expressions constructed", flush=True)

print("STAGE 2: {{B,phi1},H} on family", flush=True)
term1 = bracket_at(Bphi1, H)
print("TERM 1 =", term1, flush=True)

print("STAGE 3: {B,chi-phi2} on family", flush=True)
term2 = bracket_at(B, G)
print("TERM 2 =", term2, flush=True)

K = S.factor(term1 + term2)
print("K_D1 FAMILY =", K, flush=True)

if K == 0:
    print("PASS: K_D1 vanishes on six-parameter family")
else:
    print("RESULT: K_D1 not identically zero")
    print("FACTORED:", K, flush=True)

print("R28Y COMPLETE")
print("Ddot FAMILY CONDITION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =", round(perf_counter()-start,3))
