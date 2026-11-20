import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27I ===", flush=True)

u = S.symbols("u0:4")
q = S.symbols("q0:4")
P = S.symbols("P0:4")
p0 = S.symbols("r0:4")
p1 = S.symbols("t0:4")
I = S.symbols("IA IS IW")
pI = S.symbols("pIA pIS pIW")

IA, IS, IW = I
pIA, pIS, pIW = pI

def mdot(a, b):
    return -a[0]*b[0] + sum(
        a[i]*b[i] for i in range(1, 4)
    )

def dot(a, b):
    return sum(a[i]*b[i] for i in range(4))

s = -mdot(u, u)
v = S.sqrt(s)
d = mdot(u, q)
W = mdot(q, q) + d**2/s

A = (3 + 2*IA)*IS
C = (2 + IA + S.Rational(3,2)*IA**2)*IW

# Exact rational diagnostic parameters:
# m=hA=hS=hW=kA=kS=kW=1.
M = 1 + sum(
    (I[i]**2 + pI[i]**2)/2
    for i in range(3)
)

J = sum(
    S.diff(A, I[i])*pI[i]
    for i in range(3)
)

phi1 = dot(u, P)
phi2 = dot(q, P)
chi = dot(u, p1)
psi = dot(q, p1) + 2*A*W/v**3

theta = (
    -dot(p0, u)
    +dot(p1, q)
    +v*W*mdot(P, P)/(4*C)
    +3*A*W/v**3
    -v*M
)

Hc = (
    dot(p0, u)
    +dot(p1, q)
    -v*W*mdot(P, P)/(4*C)
    +A*W/v**3
    +v*M
)

B = (
    -dot(p0, q)
    -v*W*mdot(p1, P)/(2*C)
    +d*W*mdot(P, P)/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(
        phi2 + d*phi1/s
    )
    +2*W*J/v**2
)

def ev(expr, data):
    return S.simplify(expr.subs(data))

def witness(pIS_value):
    data = {}

    for i in range(4):
        data[u[i]] = S.Integer(2 if i == 0 else 0)
        data[q[i]] = S.Integer(1 if i == 1 else 0)
        data[P[i]] = S.Integer(1 if i == 2 else 0)
        data[p1[i]] = S.Integer(0)
        data[p0[i]] = S.Integer(0)

    data.update({
        IA:0, IS:0, IW:1,
        pIA:0,
        pIS:S.Integer(pIS_value),
        pIW:0,
    })

    mass = ev(M, data)
    jvalue = ev(J, data)
    cvalue = ev(C, data)

    # theta=0 fixes p0[0].
    # B=0 fixes p0[1].
    data[p0[0]] = (
        S.Rational(1, 2)/(4*cvalue)
        - mass
    )
    data[p0[1]] = jvalue/2

    return data

cases = [
    ("J_NONZERO", witness(1)),
    ("J_ZERO", witness(0)),
]

print("STAGE 1: verify 4D constraint surfaces", flush=True)

for label, data in cases:
    if not (
        ev(s,data)>0
        and ev(W,data)>0
        and ev(C,data)!=0
        and ev(mdot(P,P),data)>0
    ):
        raise AssertionError(f"{label}: invalid regular patch")

    for name, expr in [
        ("A",A), ("B",B),
        ("phi1",phi1), ("phi2",phi2),
        ("chi",chi), ("psi",psi),
        ("theta",theta)
    ]:
        value = ev(expr,data)
        if value != 0:
            raise AssertionError(
                f"{label}: {name}={value}"
            )

    print(
        f"PASS: {label}; "
        f"J={ev(J,data)}, "
        f"P²={ev(mdot(P,P),data)}",
        flush=True
    )

print("STAGE 1 PASS", flush=True)

print("STAGE 2: construct canonical derivatives", flush=True)

pairs = (
    list(zip(u,p1))
    +list(zip(q,P))
    +list(zip(I,pI))
)

derivatives = [
    (
        S.diff(B,x),
        S.diff(Hc,px),
        S.diff(B,px),
        S.diff(Hc,x)
    )
    for x,px in pairs
]

print("STAGE 2 PASS", flush=True)

print("STAGE 3: evaluate exact brackets", flush=True)

for label,data in cases:
    terms = []

    for Bx,Hpx,Bpx,Hx in derivatives:
        term = (
            ev(Bx,data)*ev(Hpx,data)
            -ev(Bpx,data)*ev(Hx,data)
        )
        terms.append(S.simplify(term))

    D = S.simplify(sum(terms))
    coeff = ev(2*W*J/v**2,data)

    print(label, flush=True)
    print("  D =",D,flush=True)
    print("  multiplier coefficient =",coeff,flush=True)

    if coeff != 0:
        required = S.simplify(-D/coeff)
        print(
            "  lambda2 required =",
            required,
            flush=True
        )
    elif D != 0:
        print(
            "  OBSTRUCTION: D != 0, "
            "multiplier coefficient = 0",
            flush=True
        )
    else:
        print(
            "  B preservation satisfied at witness",
            flush=True
        )

print("STAGE 3 PASS", flush=True)

print()
print("R27I COMPLETE")
print("FOUR-DIMENSIONAL WITNESSES: VERIFIED")
print("NONZERO TRANSVERSE MOMENTUM: VERIFIED")
print("D ON SINGULAR WITNESSES: EVALUATED")
print("GENERAL SINGULAR DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
