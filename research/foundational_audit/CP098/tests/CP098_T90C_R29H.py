import sympy as S
import runpy, contextlib, io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29H ===", flush=True)

folder = Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(folder / "CP098_T90C_R27U.py"))

u,q,P,p0,p1 = (ns[k] for k in ("u","q","P","p0","p1"))
B,H,phi2 = (ns[k] for k in ("B","H","phi2"))
IA,IS,IW,pA,pS,pW = (
    ns[k] for k in ("IA","IS","IW","pA","pS","pW")
)

pairs = list(zip(u,p1)) + list(zip(q,P)) + [
    (IA,pA),(IS,pS),(IW,pW)
]

z,b = S.symbols("z b", real=True)
rho = S.Integer(2)
w = S.Integer(1)

ISv = 12 + S.Rational(448,961)*rho**2
Mv = 1 + (S.Rational(9,4)+ISv**2+1+b**2)/2

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
    pS:S.Integer(0),
    pW:b
})

def at(expr):
    return S.cancel(expr.subs(point))

print("STAGE 1: canonical derivatives", flush=True)

first = [
    (
        S.diff(B,x), S.diff(B,p),
        S.diff(H,x), S.diff(H,p)
    )
    for x,p in pairs
]

directions = {}
for x,p in pairs:
    directions[x] = at(S.diff(H,p))
    directions[p] = -at(S.diff(H,x))

active = [
    (var,val) for var,val in directions.items()
    if val != 0
]

print("Active Hamiltonian directions:",len(active),flush=True)

def directional(expr):
    result = S.Integer(0)
    for var,val in active:
        result += val*at(S.diff(expr,var))
    return S.cancel(result)

print("STAGE 2: first bracket D", flush=True)

Dcalc = S.Integer(0)
for Bx,Bp,Hx,Hp in first:
    Dcalc += at(Bx)*at(Hp)-at(Bp)*at(Hx)

Dcalc = S.factor(Dcalc)

Dbase = (
    802816*rho**4
    +476656*rho**2
    +41330688*rho**2
    +19393941
    +531948096
)/S.Integer(14776336)

Dexpected = (
    Dbase + 8*z**2/31 + b**2/4
    + 16*b*rho*z/31
)

assert S.simplify(Dcalc-Dexpected) == 0
print("PASS: canonical D matches R28Q", flush=True)

print("STAGE 3: compute F_D", flush=True)

F = S.Integer(0)

for i,(Bx,Bp,Hx,Hp) in enumerate(first,1):
    contribution = (
        directional(Bx)*at(Hp)
        + at(Bx)*directional(Hp)
        - directional(Bp)*at(Hx)
        - at(Bp)*directional(Hx)
    )
    F = S.cancel(F+contribution)
    print("F_D pair",i,"of 11",flush=True)

F = S.factor(F)
print("PASS: F_D constructed",flush=True)

print("STAGE 4: compute K_D2 via Jacobi",flush=True)

Bphi2 = sum(
    S.diff(B,q[i])*q[i]
    - S.diff(B,P[i])*P[i]
    for i in range(4)
)

Hphi2 = sum(
    S.diff(H,q[i])*q[i]
    - S.diff(H,P[i])*P[i]
    for i in range(4)
)

def bracket_at(f,g,label):
    result = S.Integer(0)
    for i,(x,p) in enumerate(pairs,1):
        result += (
            at(S.diff(f,x))*at(S.diff(g,p))
            - at(S.diff(f,p))*at(S.diff(g,x))
        )
        print(label,"pair",i,"of 11",flush=True)
    return S.factor(result)

T1 = bracket_at(Bphi2,H,"T1")
T2 = bracket_at(B,Hphi2,"T2")
K = S.factor(T1+T2)

print("PASS: K_D2 constructed",flush=True)

print("STAGE 5: exact off-center witnesses",flush=True)

z2min = S.factor(
    961*Dbase/(8*(32*rho**2-31))
)
zval = S.sqrt(z2min+1)

disc = S.factor(
    -Dbase
    +8*(32*rho**2-31)*zval**2/961
)
assert S.simplify(disc-S.Rational(776,961)) == 0

center = -S.Rational(32,31)*rho*zval

branches = {
    "PLUS": center+2*S.sqrt(disc),
    "MINUS": center-2*S.sqrt(disc)
}

# Values independently established by R29G.
lam = {
    "PLUS": (
        -S.Rational(1792,103261)*S.sqrt(194)
        +S.Rational(14,3201091)*S.sqrt(144209559530)
    ),
    "MINUS": (
        S.Rational(1792,103261)*S.sqrt(194)
        +S.Rational(14,3201091)*S.sqrt(144209559530)
    )
}

for name,bval in branches.items():
    print("BRANCH:",name,flush=True)
    subs = {z:zval,b:bval}

    assert S.simplify(Dcalc.subs(subs)) == 0
    print("PASS: D=0",flush=True)

    Fv = S.simplify(F.subs(subs))
    Kv = S.simplify(K.subs(subs))

    dot = S.simplify(Fv+lam[name]*Kv)

    print("F_D =",Fv,flush=True)
    print("K_D2 =",Kv,flush=True)
    print("lambda2 =",lam[name],flush=True)
    print("Ddot =",dot,flush=True)

    if dot == 0:
        print("PASS: Ddot=0 at this witness",flush=True)
    else:
        numeric = S.N(dot,30)
        print("Ddot numerical check =",numeric,flush=True)
        if numeric.is_zero is False:
            print("OBSTRUCTION: Ddot != 0 at this witness",
                  flush=True)
        else:
            print("STOP: exact nonzero status unresolved",
                  flush=True)

print("R29H EVALUATION COMPLETE")
print("GENERAL FAMILY CLOSURE: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",round(perf_counter()-start,3))
