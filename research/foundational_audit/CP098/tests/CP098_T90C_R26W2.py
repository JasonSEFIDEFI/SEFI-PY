import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26W2 ===", flush=True)

def check(label, residual):
    result = S.factor(residual)
    if result != 0:
        print("FAIL:", label, result, flush=True)
        raise AssertionError(label)
    print("PASS:", label, flush=True)

IA, IS, IW = S.symbols("IA IS IW", real=True)
pIA, pIS, pIW = S.symbols("pIA pIS pIW", real=True)

hA, hS, hW = S.symbols(
    "hA hS hW", positive=True
)

kA, kS, kW = S.symbols(
    "kA kS kW", real=True
)

v, W = S.symbols("v W", positive=True)
P2 = S.symbols("P2", real=True)

I = (IA, IS, IW)
pI = (pIA, pIS, pIW)
h = (hA, hS, hW)
k = (kA, kS, kW)

mu1 = 3 + 2*IA
mu2 = 2 + IA + S.Rational(3,2)*IA**2

A = mu1*IS
C = mu2*IW

V = sum(k[i]*I[i]**2/2 for i in range(3))
T = sum(pI[i]**2/(2*h[i]) for i in range(3))

# Identity-dependent Hamiltonian terms.
# All other Hc terms are independent of I.

H = (
    -v*W*P2/(4*C)
    +A*W/v**3
    +v*(V+T)
)

print("STAGE 1: exact Hamilton equations", flush=True)

Idot = [S.diff(H,p) for p in pI]
pdot = [-S.diff(H,x) for x in I]

for i in range(3):
    check(
        f"identity velocity {i}",
        Idot[i]-v*pI[i]/h[i]
    )

expected_pIS = -mu1*W/v**3-v*kS*IS

check(
    "pIS evolution",
    pdot[1]-expected_pIS
)

expected_pIA = (
    -2*IS*W/v**3
    -v*W*P2*S.diff(C,IA)/(4*C**2)
    -v*kA*IA
)

check(
    "pIA evolution",
    pdot[0]-expected_pIA
)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: IS=0 component", flush=True)

component_S = S.factor(pdot[1].subs(IS,0))

check(
    "pIS-dot on IS=0",
    component_S+mu1*W/v**3
)

print("pIS-dot =", component_S)
print("For IA != -3/2, pIS-dot != 0.")
print("Therefore IS=pIS=0 is not invariant")
print("away from the component intersection.")

print("STAGE 2 PASS", flush=True)

print("STAGE 3: IA=-3/2 component", flush=True)

root = -S.Rational(3,2)

component_A = S.factor(pdot[0].subs(IA,root))

Croot = S.factor(C.subs(IA,root))
Cderiv = S.factor(S.diff(C,IA).subs(IA,root))

expected_root = (
    -2*IS*W/v**3
    -v*W*P2*Cderiv/(4*Croot**2)
    +S.Rational(3,2)*v*kA
)

check(
    "pIA-dot on mu1=0",
    component_A-expected_root
)

print("C at root =", Croot)
print("C_IA at root =", Cderiv)
print("pIA-dot =", component_A)

print("Necessary condition for IA=pIA=constant:")
print("pIA-dot must also vanish.")

print("STAGE 3 PASS", flush=True)

print("STAGE 4: intersection", flush=True)

intersection_S = S.factor(
    pdot[1].subs({IA:root,IS:0})
)

intersection_A = S.factor(
    pdot[0].subs({IA:root,IS:0})
)

check("pIS-dot at intersection",intersection_S)

print("pIS-dot =",intersection_S)
print("pIA-dot =",intersection_A)

print("Intersection invariance requires further")
print("momentum and higher-order consistency tests.")

print("STAGE 4 PASS", flush=True)

print()
print("R26W COMPLETE")
print("EXACT IDENTITY MOMENTUM EQUATIONS: PASS")
print("IS=0 COMPONENT OBSTRUCTION: VERIFIED")
print("MU1=0 COMPONENT: CONDITIONAL")
print("INTERSECTION INVARIANCE: OPEN")
print("B=0 PRESERVATION: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start,3)
)
