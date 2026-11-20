import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27L ===", flush=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pA, pS, pW = S.symbols("pA pS pW", real=True)
hA, hS, hW = S.symbols(
    "hA hS hW", positive=True
)
kA, kS, kW = S.symbols(
    "kA kS kW", positive=True
)
v, W = S.symbols("v W", positive=True)
P2 = S.symbols("P2", nonnegative=True)

I = (IA, IS, IW)
p = (pA, pS, pW)
h = (hA, hS, hW)
k = (kA, kS, kW)

alpha = 3 + 2*IA
A = alpha*IS
C = (2+IA+S.Rational(3,2)*IA**2)*IW

# Identity-dependent part of Hc.
# Other terms have zero identity derivatives.
H = (
    -v*W*P2/(4*C)
    +A*W/v**3
    +v*sum(
        (k[i]*I[i]**2+p[i]**2/h[i])/2
        for i in range(3)
    )
)

def PB_identity(F,G):
    return sum(
        S.diff(F,I[i])*S.diff(G,p[i])
        -S.diff(F,p[i])*S.diff(G,I[i])
        for i in range(3)
    )

def check(label, expr):
    residual = S.cancel(expr)
    if residual != 0:
        raise AssertionError(
            f"{label}: {residual}"
        )
    print("PASS:",label,flush=True)

print("STAGE 1: general identity dynamics",flush=True)

J = (
    2*IS*pA/hA
    +alpha*pS/hS
)

Adot = PB_identity(A,H)
Jdot = PB_identity(J,H)

check("Adot = vJ",Adot-v*J)

pAdot = -S.diff(H,IA)
pSdot = -S.diff(H,IS)

Jdot_expected = (
    4*v*pA*pS/(hA*hS)
    +2*IS*pAdot/hA
    +alpha*pSdot/hS
)

check("general Jdot identity",Jdot-Jdot_expected)
print("STAGE 1 PASS",flush=True)

print("STAGE 2: IS=0 component",flush=True)

# If alpha != 0, J=0 forces pS=0.
branch_S = {IS:0,pS:0}

result_S = S.factor(Jdot.subs(branch_S))
expected_S = -alpha**2*W/(hS*v**3)

check(
    "IS=0, J=0: Jdot",
    result_S-expected_S
)

print("Jdot =",result_S,flush=True)
print(
    "For alpha != 0 and W,v,hS > 0, "
    "Jdot is strictly negative.",
    flush=True
)
print("STAGE 2 PASS",flush=True)

print("STAGE 3: IA=-3/2 component",flush=True)

# If IS != 0, J=0 forces pA=0.
branch_A = {
    IA:-S.Rational(3,2),
    pA:0
}

result_A = S.factor(Jdot.subs(branch_A))
expected_A = S.factor(
    2*IS*pAdot.subs(branch_A)/hA
)

check(
    "IA=-3/2, J=0: Jdot",
    result_A-expected_A
)

print("Jdot =",result_A,flush=True)
print(
    "Necessary tangency condition: pAdot=0 "
    "when IS != 0.",
    flush=True
)
print("STAGE 3 PASS",flush=True)

print("STAGE 4: intersection",flush=True)

intersection = {
    IA:-S.Rational(3,2),
    IS:0
}

result_X = S.factor(Jdot.subs(intersection))

check(
    "intersection Jdot",
    result_X-4*v*pA*pS/(hA*hS)
)

print("Jdot =",result_X,flush=True)
print(
    "Necessary tangency condition: pA*pS=0.",
    flush=True
)
print("STAGE 4 PASS",flush=True)

print()
print("R27L COMPLETE")
print("GENERAL IDENTITY-SECTOR JDOT: DERIVED")
print("IS=0 NONINTERSECTION: SECOND-ORDER OBSTRUCTION")
print("IA=-3/2 COMPONENT: CONDITIONAL")
print("INTERSECTION: CONDITIONAL")
print("B=0 PRESERVATION ON THESE COMPONENTS: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
