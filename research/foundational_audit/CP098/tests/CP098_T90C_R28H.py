import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28H ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R27U.py")
assert source.is_file()

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(source))

u,q,P = ns["u"],ns["q"],ns["P"]
p0,p1 = ns["p0"],ns["p1"]

pairs = (
    list(zip(u,p1))
    + list(zip(q,P))
    + [(ns["IA"],ns["pA"]),
       (ns["IS"],ns["pS"]),
       (ns["IW"],ns["pW"])]
)

print("STAGE 1: canonical derivatives", flush=True)

derivs = [
    (S.diff(ns["B"],x), S.diff(ns["B"],p),
     S.diff(ns["H"],x), S.diff(ns["H"],p))
    for x,p in pairs
]
assert len(derivs)==11
print("PASS: 11 canonical pairs", flush=True)

constraints = [
    ("phi1",ns["phi1"]), ("phi2",ns["phi2"]),
    ("chi",ns["chi"]), ("psi",ns["psi"]),
    ("theta",ns["theta"]), ("B",ns["B"]),
    ("A",ns["A"]), ("J",ns["J"]),
    ("Q",ns["Q"]), ("pA",ns["pA"])
]

def witness(w,rho,z,t,a,b):
    w,rho,z,t,a,b = map(
        S.Rational,(w,rho,z,t,a,b)
    )
    assert w!=0

    point = dict(ns["point"])
    point[q[1]]=w
    point[P[2]]=rho
    point[p1[2]]=z
    point[p1[3]]=t
    point[ns["pS"]]=a
    point[ns["pW"]]=b

    W=w**2
    r2=rho**2
    IS=12/W+S.Rational(448,961)*r2
    point[ns["IS"]]=IS

    M=1+(
        S.Rational(9,4)+IS**2+1
        +a**2+b**2
    )/2

    point[p0[0]]=W*r2/S.Rational(31,2)-M
    point[p0[1]]=-S.Rational(8,31)*w*rho*z

    return point

def at(expr,point):
    return S.cancel(expr.subs(point))

def baseline(w,rho,z,t):
    W=w**2
    r2=rho**2
    N=(
        476656*W**3*r2
        +802816*W**2*r2**2
        +19393941*W**2
        +41330688*W*r2
        +531948096
    )
    return S.cancel(
        N/(14776336*W)
        +S.Rational(8,31)*W*(z**2+t**2)
    )

print("STAGE 2: independent identity witnesses",
      flush=True)

cases = [
    (1,1,0,0,0,0),
    (1,1,0,0,1,0),
    (1,1,0,0,0,2),
    (1,1,0,0,1,2),
    (2,1,1,1,1,1),
    (S.Rational(1,2),2,1,0,1,2)
]

for case in cases:
    w,rho,z,t,a,b=map(S.Rational,case)
    point=witness(w,rho,z,t,a,b)

    for name,expr in constraints:
        residual=at(expr,point)
        assert residual==0,(case,name,residual)

    contributions=[]
    for Bx,Bp,Hx,Hp in derivs:
        contribution=S.cancel(
            at(Bx,point)*at(Hp,point)
            -at(Bp,point)*at(Hx,point)
        )
        contributions.append(contribution)

    D=S.cancel(sum(contributions))
    D0=baseline(w,rho,z,t)

    print("case =",case,flush=True)
    print("PASS: ten constraints",flush=True)
    print("D =",D,flush=True)
    print("identity correction =",S.factor(D-D0),
          flush=True)
    print("sign =",S.sign(D),flush=True)

    if a==0 and b==0:
        assert D==D0

    if w==1 and rho==1 and z==0 and t==0:
        expected=D0+(a**2+b**2)/4
        assert D==expected,(
            "R27W2 regression failed",case,D,expected
        )

print("STAGE 2 PASS",flush=True)

print()
print("R28H COMPLETE")
print("IDENTITY-MOMENTUM WITNESSES: CHECKED")
print("GLOBAL IDENTITY-MOMENTUM SIGN: OPEN")
print("GENERAL IA BRANCH: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(perf_counter()-start,3))
