import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28G ===", flush=True)

source = Path(__file__).with_name(
    "CP098_T90C_R27U.py"
)
assert source.is_file(), "R27U source missing"

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

w = S.symbols("w", real=True, nonzero=True)
rho, z, t = S.symbols("rho z t", real=True)

u, q, P = ns["u"], ns["q"], ns["P"]
p0, p1 = ns["p0"], ns["p1"]

point = dict(ns["point"])

point[q[1]] = w
point[P[2]] = rho
point[p1[2]] = z
point[p1[3]] = t

W0 = w**2
P20 = rho**2

IS0 = 12/W0 + S.Rational(448,961)*P20
point[ns["IS"]] = IS0

IA0 = -S.Rational(3,2)
C0 = S.Rational(31,8)
M0 = 1 + (IA0**2 + IS0**2 + 1)/2

point[p0[0]] = W0*P20/(4*C0) - M0
point[p0[1]] = -S.Rational(8,31)*w*rho*z

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: exact constraint audit",
      flush=True)

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
    assert residual == 0, (name,residual)
    print("PASS:",name,flush=True)

print("STAGE 1 PASS",flush=True)

print("STAGE 2: canonical bracket",flush=True)

pairs = (
    list(zip(u,p1))
    + list(zip(q,P))
    + [
        (ns["IA"],ns["pA"]),
        (ns["IS"],ns["pS"]),
        (ns["IW"],ns["pW"])
    ]
)

assert len(pairs) == 11
parts = []

for i,(x,p) in enumerate(pairs,1):
    bx = at(S.diff(ns["B"],x))
    bp = at(S.diff(ns["B"],p))
    hx = at(S.diff(ns["H"],x))
    hp = at(S.diff(ns["H"],p))

    term = S.cancel(bx*hp-bp*hx)
    parts.append(term)
    print(f"pair {i:02d} = {term}",flush=True)

D = S.cancel(sum(parts))
print("D(w,rho,z,t) =",S.factor(D),flush=True)
print("STAGE 2 PASS",flush=True)

print("STAGE 3: exact regression identities",
      flush=True)

R28C = (
    476656*W0**3*P20
    +802816*W0**2*P20**2
    +19393941*W0**2
    +41330688*W0*P20
    +531948096
)/(14776336*W0)

R28E = R28C + S.Rational(8,31)*W0*t**2

R28F = (
    S.Rational(50176,923521)*rho**4
    +S.Rational(2719,961)*rho**2
    +S.Rational(8,31)*z**2
    +S.Rational(597,16)
)

assert S.cancel(D.subs(z,0)-R28E) == 0
print("PASS: R28E recovered at z=0",
      flush=True)

assert S.cancel(
    D.subs({w:1,t:0})-R28F
) == 0
print("PASS: R28F recovered at w=1,t=0",
      flush=True)

print("STAGE 3 PASS",flush=True)

print("STAGE 4: new correction structure",
      flush=True)

correction = S.factor(D-R28E)
print("D minus R28E =",correction)

poly_z = S.Poly(correction,z)
print("degree in z =",poly_z.degree())
print("z coefficients =",poly_z.all_coeffs())

# Do not assume the new correction is positive.
# Test for mixed terms before making any sign claim.
numerator,denominator = S.fraction(
    S.cancel(D)
)

print("numerator =",S.factor(numerator))
print("denominator =",S.factor(denominator))

poly = S.Poly(numerator,w,rho,z,t)
print("total polynomial degree =",poly.total_degree())
print("number of terms =",len(poly.terms()))

even = all(
    all(power % 2 == 0 for power in powers)
    for powers,coeff in poly.terms()
)

print("all numerator exponents even =",even)

if even:
    coeffs = poly.coeffs()
    print("all present coefficients positive =",
          all(c>0 for c in coeffs))
else:
    print("SIGN CERTIFICATE: further analysis needed")

print("STAGE 4 COMPLETE",flush=True)

print()
print("R28G COMPLETE")
print("FOUR-VARIABLE FAMILY: DERIVED")
print("R28E/R28F REGRESSIONS: VERIFIED")
print("GLOBAL SIGN: REQUIRES INTERPRETATION")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
