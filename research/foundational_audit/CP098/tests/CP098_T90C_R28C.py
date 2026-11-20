import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28C ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R27U.py")
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

w, rho = S.symbols("w rho", real=True)
W, z = S.symbols("W z", positive=True)

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

point = dict(ns["point"])
point[q[1]] = w
point[P[2]] = rho
point[p1[3]] = S.Integer(0)

IS0 = 12/w**2 + S.Rational(448,961)*rho**2
point[ns["IS"]] = IS0

IA0 = -S.Rational(3,2)
C0 = S.Rational(31,8)
M0 = 1 + (IA0**2 + IS0**2 + 1)/2
point[p0[0]] = w**2*rho**2/(4*C0) - M0

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: joint symbolic constraints", flush=True)

constraints = [
    ("phi1",ns["phi1"]),
    ("phi2",ns["phi2"]),
    ("chi",ns["chi"]),
    ("psi",ns["psi"]),
    ("theta",ns["theta"]),
    ("B",ns["B"]),
    ("A",ns["A"]),
    ("J",ns["J"]),
    ("Q",ns["Q"]),
    ("pA",ns["pA"])
]

for name, expr in constraints:
    residual = at(expr)
    assert residual == 0, (name, residual)
    print("PASS:", name, flush=True)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: canonical D(w,rho)", flush=True)

pairs = (
    list(zip(u,p1))
    + list(zip(q,P))
    + [
        (ns["IA"],ns["pA"]),
        (ns["IS"],ns["pS"]),
        (ns["IW"],ns["pW"])
    ]
)

parts = []

for i,(x,p) in enumerate(pairs,1):
    bx = at(S.diff(ns["B"],x))
    bp = at(S.diff(ns["B"],p))
    hx = at(S.diff(ns["H"],x))
    hp = at(S.diff(ns["H"],p))

    term = S.cancel(bx*hp-bp*hx)
    parts.append(term)
    print(f"pair {i:02d} = {term}", flush=True)

D = S.cancel(sum(parts))
print("D(w,rho) =", S.factor(D), flush=True)

# Exact regression against R27Y (w=1).
R27Y = (
    S.Rational(50176,923521)*rho**4
    + S.Rational(2719,961)*rho**2
    + S.Rational(597,16)
)

assert S.cancel(D.subs(w,1)-R27Y) == 0
print("PASS: full R27Y polynomial recovered", flush=True)

# Exact regression against R28A (rho=1).
R28A = (
    476656*w**6
    + 20196757*w**4
    + 41330688*w**2
    + 531948096
)/(14776336*w**2)

assert S.cancel(D.subs(rho,1)-R28A) == 0
print("PASS: full R28A rational function recovered",
      flush=True)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: joint-domain sign certificate",
      flush=True)

num, den = S.fraction(S.cancel(D))
num = S.expand(num)
den = S.expand(den)

np = S.Poly(num,w,rho)
dp = S.Poly(den,w,rho)

def even_exponents(poly):
    return all(
        a % 2 == 0 and b % 2 == 0
        for (a,b),coeff in poly.terms()
    )

if even_exponents(np) and even_exponents(dp):
    nWz = S.Poly.from_dict(
        {
            (a//2,b//2):coeff
            for (a,b),coeff in np.terms()
        },
        (W,z)
    ).as_expr()

    dWz = S.Poly.from_dict(
        {
            (a//2,b//2):coeff
            for (a,b),coeff in dp.terms()
        },
        (W,z)
    ).as_expr()

    assert S.cancel(
        D-(nWz/dWz).subs({W:w**2,z:rho**2})
    ) == 0

    print("numerator(W,z) =",S.factor(nWz))
    print("denominator(W,z) =",S.factor(dWz))

    ncoeff = S.Poly(nWz,W,z).coeffs()
    dcoeff = S.Poly(dWz,W,z).coeffs()

    # Missing monomials have zero coefficients.
    # Only present monomials need positive coefficients.
    n_positive = all(c>0 for c in ncoeff)
    d_positive = all(c>0 for c in dcoeff)

    n_constant_W = S.cancel(nWz.subs(W,0))
    n_constant_z = S.cancel(nWz.subs(z,0))
    d_at_z0 = S.cancel(dWz.subs(z,0))

    if (
        n_positive
        and d_positive
        and n_constant_z != 0
        and d_at_z0 != 0
    ):
        print("PASS: D>0 for W>0 and z>=0")
        print("PASS: no D=0 point in joint family")
    else:
        print("SIGN NOT CERTIFIED")
        print("Further exact sign analysis required")
else:
    print("EVEN-EXPONENT TEST NOT APPLICABLE")
    print("Further exact sign analysis required")

print("STAGE 3 COMPLETE", flush=True)

print()
print("R28C COMPLETE")
print("JOINT SYMBOLIC CONSTRAINTS: CHECKED")
print("R27Y/R28A REGRESSIONS: CHECKED")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
