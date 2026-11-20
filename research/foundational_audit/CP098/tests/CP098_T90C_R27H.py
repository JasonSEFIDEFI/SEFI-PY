import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27H ===", flush=True)

u = S.symbols("u0:2")
q = S.symbols("q0:2")
P = S.symbols("P0:2")
p0 = S.symbols("r0:2")
p1 = S.symbols("t0:2")
I = S.symbols("IA IS IW")
pI = S.symbols("pIA pIS pIW")

IA, IS, IW = I
pIA, pIS, pIW = pI

def mdot(a,b):
    return -a[0]*b[0]+a[1]*b[1]

def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

s = -mdot(u,u)
v = S.sqrt(s)
d = mdot(u,q)
W = mdot(q,q)+d**2/s

A = (3+2*IA)*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW

# Fixed rational parameters for this diagnostic:
# m=1, hA=hS=hW=1, kA=kS=kW=1.
M = 1 + sum(
    (I[i]**2+pI[i]**2)/2
    for i in range(3)
)
J = sum(
    S.diff(A,I[i])*pI[i]
    for i in range(3)
)

phi1 = dot(u,P)
phi2 = dot(q,P)
chi = dot(u,p1)
psi = dot(q,p1)+2*A*W/v**3

theta = (
    -dot(p0,u)+dot(p1,q)
    +v*W*mdot(P,P)/(4*C)
    +3*A*W/v**3-v*M
)

Hc = (
    dot(p0,u)+dot(p1,q)
    -v*W*mdot(P,P)/(4*C)
    +A*W/v**3+v*M
)

B = (
    -dot(p0,q)
    -v*W*mdot(p1,P)/(2*C)
    +d*W*mdot(P,P)/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(
        phi2+d*phi1/s
    )
    +2*W*J/v**2
)

pairs = (
    list(zip(u,p1))
    +list(zip(q,P))
    +list(zip(I,pI))
)

print("STAGE 1: prepare exact witnesses", flush=True)

def witness(pIS_value):
    # On this chosen surface:
    # u=(2,0), q=(0,1), P=p1=0.
    # theta=0 gives p0[0]=-M.
    # B=0 gives p0[1]=J/2.
    data = {
        u[0]:2, u[1]:0,
        q[0]:0, q[1]:1,
        P[0]:0, P[1]:0,
        p1[0]:0, p1[1]:0,
        IA:0, IS:0, IW:1,
        pIA:0, pIS:S.Integer(pIS_value), pIW:0,
    }
    mass = S.simplify(M.subs(data))
    jvalue = S.simplify(J.subs(data))
    data[p0[0]] = -mass
    data[p0[1]] = jvalue/2
    return data

cases = [
    ("J_NONZERO", witness(1)),
    ("J_ZERO", witness(0))
]

def ev(expr, data):
    return S.simplify(expr.subs(data))

def check_zero(label, expr, data):
    value = ev(expr,data)
    if value != 0:
        raise AssertionError(f"{label}: {value}")

for label,data in cases:
    assert ev(s,data)>0
    assert ev(W,data)>0
    assert ev(C,data)!=0
    for name,expr in [
        ("A",A),("B",B),
        ("phi1",phi1),("phi2",phi2),
        ("chi",chi),("psi",psi),
        ("theta",theta)
    ]:
        check_zero(f"{label} {name}",expr,data)
    print(
        f"PASS: {label}; J={ev(J,data)}",
        flush=True
    )

print("STAGE 1 PASS",flush=True)

print("STAGE 2: differentiate canonical brackets",flush=True)

# Differentiate once and evaluate each contribution
# before summing. No global expansion or simplification.
D_terms = [
    (
        S.diff(B,x), S.diff(Hc,px),
        S.diff(B,px), S.diff(Hc,x)
    )
    for x,px in pairs
]

print("STAGE 2 PASS",flush=True)

print("STAGE 3: evaluate D on both surfaces",flush=True)

for label,data in cases:
    contributions = []
    for Bx,Hpx,Bpx,Hx in D_terms:
        contribution = (
            ev(Bx,data)*ev(Hpx,data)
            -ev(Bpx,data)*ev(Hx,data)
        )
        contributions.append(S.simplify(contribution))

    D = S.simplify(sum(contributions))
    Jval = ev(J,data)
    coefficient = ev(2*W*J/v**2,data)

    print(f"{label}:",flush=True)
    print("  D =",D,flush=True)
    print("  J =",Jval,flush=True)
    print("  lambda2 coefficient =",coefficient,flush=True)

    if coefficient != 0:
        multiplier = S.simplify(-D/coefficient)
        print(
            "  lambda2 required by Bdot=0 =",
            multiplier,
            flush=True
        )
    elif D != 0:
        print(
            "  Bdot=0 OBSTRUCTION: D is nonzero "
            "and lambda2 cannot cancel it",
            flush=True
        )
    else:
        print(
            "  Bdot=0 holds at this witness "
            "without fixing lambda2",
            flush=True
        )

print("STAGE 3 PASS",flush=True)

print()
print("R27H COMPLETE")
print("SINGULAR WITNESSES: VERIFIED")
print("D ON TWO SINGULAR WITNESSES: EVALUATED")
print("GENERAL 4D SINGULAR CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
