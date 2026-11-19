import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28K ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R27U.py")
assert source.is_file()

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

u,q,P = ns["u"],ns["q"],ns["P"]
p0,p1 = ns["p0"],ns["p1"]

pairs = (
    list(zip(u,p1))
    + list(zip(q,P))
    + [
        (ns["IA"],ns["pA"]),
        (ns["IS"],ns["pS"]),
        (ns["IW"],ns["pW"])
    ]
)
assert len(pairs)==11

z = S.sqrt(S.Rational(593952197,123008))
b = -S.Rational(32,31)*z

point = dict(ns["point"])
point[q[1]] = S.Integer(1)
point[P[2]] = S.Integer(1)
point[p1[2]] = z
point[p1[3]] = S.Integer(0)
point[ns["pA"]] = S.Integer(0)
point[ns["pS"]] = S.Integer(0)
point[ns["pW"]] = b

IS0 = 12 + S.Rational(448,961)
point[ns["IS"]] = IS0

M0 = 1 + (
    S.Rational(9,4) + IS0**2 + 1 + b**2
)/2

point[p0[0]] = S.Rational(2,31)-M0
point[p0[1]] = -S.Rational(8,31)*z

def at(expr):
    return S.simplify(expr.subs(point))

def bracket_at(f,g):
    terms=[]
    for x,p in pairs:
        fx=at(S.diff(f,x))
        fp=at(S.diff(f,p))
        gx=at(S.diff(g,x))
        gp=at(S.diff(g,p))
        terms.append(S.cancel(fx*gp-fp*gx))
    return S.simplify(sum(terms))

print("STAGE 1: exact witness",flush=True)

for name in (
    "phi1","phi2","chi","psi","theta",
    "B","A","J","Q","pA"
):
    residual=at(ns[name])
    assert residual==0,(name,residual)
    print("PASS:",name,flush=True)

D=bracket_at(ns["B"],ns["H"])
assert D==0,("D",D)
print("PASS: D=0",flush=True)
print("STAGE 1 PASS",flush=True)

print("STAGE 2: Q preservation",flush=True)

F=bracket_at(ns["Q"],ns["H"])
L1=bracket_at(ns["Q"],ns["phi1"])
L2=bracket_at(ns["Q"],ns["phi2"])

print("{Q,Hc} =",F,flush=True)
print("{Q,phi1} =",L1,flush=True)
print("{Q,phi2} =",L2,flush=True)

# R27Q/R27R derived that lambda1 drops out
# on this branch. Verify rather than assume.
assert L1==0,("lambda1 coupling survives",L1)

# Independent coefficient regression from R27R.
expected_L2 = -(
    2883*1*S.Rational(1,1)
    +112*1*1
)/961
assert S.simplify(L2-expected_L2)==0,(
    "R27R multiplier coefficient mismatch",
    L2,expected_L2
)
assert L2!=0

lambda2=S.simplify(-F/L2)
assert S.simplify(F+lambda2*L2)==0

print("PASS: Qdot=0 fixes lambda2",flush=True)
print("lambda2 =",lambda2,flush=True)
print("STAGE 2 PASS",flush=True)

print("STAGE 3: existing consistency checks",
      flush=True)

# On the exact-zero witness Bdot=D, because
# the previously derived K1 and K2 vanish.
K1=bracket_at(ns["B"],ns["phi1"])
K2=bracket_at(ns["B"],ns["phi2"])

print("{B,phi1} =",K1,flush=True)
print("{B,phi2} =",K2,flush=True)

assert K1==0
assert K2==0
assert D==0

print("PASS: Bdot=0 independently of multipliers",
      flush=True)
print("STAGE 3 PASS",flush=True)

print()
print("R28K COMPLETE")
print("EXACT D=0 WITNESS: VERIFIED")
print("Q PRESERVATION: CHECKED")
print("B PRESERVATION: CHECKED")
print("Ddot PRESERVATION: NOT YET CHECKED")
print("OTHER HIGHER TANGENCY CONDITIONS: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-start,3))
