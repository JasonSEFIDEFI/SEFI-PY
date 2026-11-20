
import sympy as S
import runpy, contextlib, io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29A ===", flush=True)

folder = Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(folder/"CP098_T90C_R27U.py"))

u,q,P,p0,p1 = (ns[k] for k in ("u","q","P","p0","p1"))
B,H,phi2 = (ns[k] for k in ("B","H","phi2"))
IA,IS,IW,pA,pS,pW = (
    ns[k] for k in ("IA","IS","IW","pA","pS","pW")
)

pairs = list(zip(u,p1))+list(zip(q,P))+[
    (IA,pA),(IS,pS),(IW,pW)
]

rho,z = S.symbols("rho z",real=True)
w = S.Integer(1)
a = S.Integer(0)
b = -32*rho*z/31

ISv = 12+S.Rational(448,961)*rho**2
Mv = 1+(S.Rational(9,4)+ISv**2+1+b**2)/2

point = {}
for vec,values in (
    (u,[2,0,0,0]),
    (q,[0,1,0,0]),
    (P,[0,0,rho,0]),
    (p1,[0,0,z,0]),
    (p0,[2*rho**2/31-Mv,-8*rho*z/31,0,0])
):
    point.update({
        sym:S.sympify(val)
        for sym,val in zip(vec,values)
    })

point.update({
    IA:S.Rational(-3,2),
    IS:ISv,
    IW:S.Integer(1),
    pA:S.Integer(0),
    pS:a,
    pW:b
})

def at(expr):
    return S.cancel(expr.subs(point))

def bracket_at(f,g,label):
    result=S.Integer(0)
    for i,(x,p) in enumerate(pairs,1):
        result += (
            at(S.diff(f,x))*at(S.diff(g,p))
            - at(S.diff(f,p))*at(S.diff(g,x))
        )
        print(label,"pair",i,"of 11",flush=True)
    return S.factor(result)

print("STAGE 1: unrestricted inner brackets",flush=True)

Bphi2 = sum(
    S.diff(B,q[i])*q[i]
    - S.diff(B,P[i])*P[i]
    for i in range(4)
)

# phi2 = q.P generates delta(q)=q,
# delta(P)=-P, with all other canonical
# coordinates unchanged.
Hphi2 = sum(
    S.diff(H,q[i])*q[i]
    - S.diff(H,P[i])*P[i]
    for i in range(4)
)

print("PASS: inner brackets constructed",flush=True)

print("STAGE 2: Jacobi term 1",flush=True)
T1=bracket_at(Bphi2,H,"T1")
print("T1 =",T1,flush=True)

print("STAGE 3: Jacobi term 2",flush=True)
T2=bracket_at(B,Hphi2,"T2")
print("T2 =",T2,flush=True)

K=S.factor(T1+T2)
print("K_D2 =",K,flush=True)

print("STAGE 4: exact D=0 substitution",flush=True)

Dbase=(
    802816*rho**4+476656*rho**2
    +41330688*rho**2+19393941+531948096
)/14776336

z2=S.factor(961*Dbase/(8*(32*rho**2-31)))

# Reduce rational polynomial dependence on z
# modulo z^2-z2, preserving both signs of z.
num,den=S.fraction(S.cancel(K))
num_reduced=S.rem(
    S.Poly(num,z),
    S.Poly(z**2-z2,z),
).as_expr()
den_reduced=S.rem(
    S.Poly(den,z),
    S.Poly(z**2-z2,z),
).as_expr()

print("K_D2 reduced numerator =",S.factor(num_reduced),flush=True)
print("K_D2 reduced denominator =",S.factor(den_reduced),flush=True)
print("D=0 relation: z^2 =",z2,flush=True)

print("R29A COMPLETE")
print("DOMAIN: rho^2 > 31/32, z^2=z2")
print("F_D: OPEN")
print("Ddot: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",round(perf_counter()-start,3))
