import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26V ===", flush=True)

IA, IS = S.symbols("IA IS", real=True)
pIA, pIS = S.symbols("pIA pIS", real=True)
hA, hS, v = S.symbols(
    "hA hS v", positive=True
)

def check(label, expression):
    result = S.factor(expression)
    if result != 0:
        print("FAIL:", label, result, flush=True)
        raise AssertionError(label)
    print("PASS:", label, flush=True)

A = (3 + 2*IA)*IS

# Hamiltonian identity-sector evolution.
IA_dot = v*pIA/hA
IS_dot = v*pIS/hS

A_dot = (
    S.diff(A,IA)*IA_dot
    + S.diff(A,IS)*IS_dot
)

print("STAGE 1: exact evolution", flush=True)

expected = v*(
    2*IS*pIA/hA
    +(3+2*IA)*pIS/hS
)

check("general A-dot identity", A_dot-expected)
print("A-dot =", S.factor(A_dot))
print("STAGE 1 PASS", flush=True)

print("STAGE 2: I_S=0 component", flush=True)

on_IS_zero = S.factor(A_dot.subs(IS,0))

check(
    "I_S=0 tangency expression",
    on_IS_zero-v*(3+2*IA)*pIS/hS
)

# Away from the intersection IA=-3/2,
# tangency requires pIS=0.
check(
    "pIS=0 makes A-dot vanish",
    on_IS_zero.subs(pIS,0)
)

print("Away from IA=-3/2:")
print("  Necessary tangency: pIS=0")
print("STAGE 2 PASS", flush=True)

print("STAGE 3: mu1=0 component", flush=True)

on_mu_zero = S.factor(
    A_dot.subs(IA,-S.Rational(3,2))
)

check(
    "mu1=0 tangency expression",
    on_mu_zero-2*v*IS*pIA/hA
)

check(
    "pIA=0 makes A-dot vanish",
    on_mu_zero.subs(pIA,0)
)

print("Away from IS=0:")
print("  Necessary tangency: pIA=0")
print("STAGE 3 PASS", flush=True)

print("STAGE 4: intersection", flush=True)

intersection = A_dot.subs({
    IA:-S.Rational(3,2),
    IS:0
})

check("A-dot vanishes at intersection",
      intersection)

print("At IA=-3/2 and IS=0:")
print("  A-dot=0 automatically")
print("  Component tangency and higher derivatives OPEN")
print("STAGE 4 PASS", flush=True)

print()
print("R26V COMPLETE")
print("FIRST-ORDER A=0 TANGENCY: VERIFIED")
print("FULL INVARIANCE: NOT PROVEN")
print("B=0 PRESERVATION: OPEN")
print("HIGHER-ORDER CONSISTENCY: OPEN")
print("GAUGE GENERATORS: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(monotonic()-start,3)
)
