
import sympy as S
import runpy, contextlib, io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28Z ===", flush=True)

folder = Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(folder/"CP098_T90C_R27U.py"))

u,q,P,p0,p1 = (ns[k] for k in ("u","q","P","p0","p1"))
H,phi1,phi2 = (ns[k] for k in ("H","phi1","phi2"))
IA,IS,IW,pA,pS,pW = (ns[k] for k in
                      ("IA","IS","IW","pA","pS","pW"))

pairs = list(zip(u,p1))+list(zip(q,P))+[
    (IA,pA),(IS,pS),(IW,pW)
]

eta = [-1,1,1,1]
s = -sum(eta[i]*u[i]**2 for i in range(4))
v = S.sqrt(s)
d = sum(eta[i]*u[i]*q[i] for i in range(4))
W = sum(eta[i]*q[i]**2 for i in range(4))+d**2/s
PP = sum(eta[i]*P[i]**2 for i in range(4))

Q = S.Rational(3,2)*v + 56*v*W*PP/(961*IW) - 2*IS*W/v**3

w,rho,z,t,a,b = S.symbols("w rho z t a b",real=True)
ISv = 12/w**2 + S.Rational(448,961)*rho**2
Mv = 1+(S.Rational(9,4)+ISv**2+1+a**2+b**2)/2

point = {}
for vec,values in (
    (u,[2,0,0,0]),
    (q,[0,w,0,0]),
    (P,[0,0,rho,0]),
    (p1,[0,0,z,t]),
    (p0,[2*w**2*rho**2/31-Mv,-8*w*rho*z/31,0,0])
):
    point.update({sym:S.sympify(val) for sym,val in zip(vec,values)})
point.update({
    IA:S.Rational(-3,2), IS:ISv, IW:S.Integer(1),
    pA:S.Integer(0),pS:a,pW:b
})

def at(expr):
    return S.factor(S.cancel(expr.subs(point)))

print("STAGE 1: Q constraint",flush=True)
assert at(Q)==0
print("PASS: Q=0 throughout family",flush=True)

def bracket_at(f,g,label):
    result=S.Integer(0)
    for i,(x,p) in enumerate(pairs,1):
        result += (
            at(S.diff(f,x))*at(S.diff(g,p))
            - at(S.diff(f,p))*at(S.diff(g,x))
        )
        print(label,"pair",i,"of",len(pairs),flush=True)
    return S.factor(result)

print("STAGE 2: Q tangency coefficients",flush=True)
F=bracket_at(Q,H,"F_Q")
K1=bracket_at(Q,phi1,"K_Q1")
K2=bracket_at(Q,phi2,"K_Q2")

print("F_Q =",F,flush=True)
print("K_Q1 =",K1,flush=True)
print("K_Q2 =",K2,flush=True)

if K1==0 and K2!=0:
    lam2=S.factor(-F/K2)
    assert S.simplify(F+lam2*K2)==0
    print("lambda2 FAMILY =",lam2,flush=True)
    print("PASS: Qdot fixes lambda2 wherever K_Q2 != 0")
elif K1==0 and K2==0:
    print("DEGENERATE: inspect F_Q before proceeding")
else:
    print("COUPLED MULTIPLIERS: lambda1 remains in Qdot")

print("R28Z COMPLETE")
print("Ddot FAMILY: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",round(perf_counter()-start,3))
