import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27F2 ===", flush=True)

v, W, C = S.symbols("v W C", positive=True)
A, d, M, J = S.symbols("A d M J", real=True)
e, f, g, P2, r = S.symbols(
    "e f g P2 r", real=True
)
phi1, phi2 = S.symbols("phi1 phi2", real=True)
lam1, lam2, D = S.symbols("lam1 lam2 D")

b = 2*A*W/v**3

psi = r+b

theta = (
    -g+r+v*W*P2/(4*C)
    +3*A*W/v**3-v*M
)

B = (
    -e
    -v*W*f/(2*C)
    +d*W*P2/(4*C*v)
    +5*A*d*W/v**5
    +d*M/v
    -2*A*W/(C*v**2)*(
        phi2+d*phi1/v**2
    )
    +2*W*J/v**2
)

K1 = (
    -g
    +v*W*P2/(4*C)
    -5*A*W/v**3
    -v*M
)

K2 = (
    -e
    -v*W*f/(2*C)
    +d*W*P2/(4*C*v)
    +15*A*d*W/v**5
    +d*M/v
    -4*A*W/(C*v**2)*(
        phi2+d*phi1/v**2
    )
    +4*W*J/v**2
)

def check(name, expression):
    residual = S.cancel(expression)
    if residual != 0:
        raise AssertionError(
            f"{name}: residual = {residual}"
        )
    print("PASS:", name, flush=True)

print("STAGE 1: corrected exact identities", flush=True)

check(
    "K1 = theta - psi - 3b",
    K1-(theta-psi-3*b)
)

check(
    "K2 relative to B",
    K2-B-(
        10*A*d*W/v**5
        -2*A*W/(C*v**2)*(
            phi2+d*phi1/v**2
        )
        +2*W*J/v**2
    )
)

print("STAGE 1 PASS", flush=True)

print("STAGE 2: A=0 reduction", flush=True)

check(
    "K1 = theta - psi at A=0",
    (K1-theta+psi).subs(A,0)
)

check(
    "K2 = B + 2WJ/v² at A=0",
    (K2-B-2*W*J/v**2).subs(A,0)
)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: constraint-surface solutions", flush=True)

# At A=0, psi=0 implies r=0.
g_solutions = S.solve(
    theta.subs({A:0,r:0}), g
)

if len(g_solutions) != 1:
    raise AssertionError("Unexpected theta solution")

e_solutions = S.solve(B.subs(A,0), e)

if len(e_solutions) != 1:
    raise AssertionError("Unexpected B solution")

surface = {
    A:0,
    r:0,
    g:g_solutions[0],
    e:e_solutions[0],
    phi1:0,
    phi2:0,
}

check("psi=0", psi.subs(surface))
check("theta=0", theta.subs(surface))
check("B=0", B.subs(surface))
check("phi1=0", phi1.subs(surface))
check("phi2=0", phi2.subs(surface))

print("STAGE 3 PASS", flush=True)

print("STAGE 4: multiplier dependence", flush=True)

check(
    "K1 vanishes",
    K1.subs(surface)
)

check(
    "K2 retains 2WJ/v²",
    K2.subs(surface)-2*W*J/v**2
)

Bdot = D+lam1*K1+lam2*K2

check(
    "restricted Bdot",
    Bdot.subs(surface)
    -(D+lam2*2*W*J/v**2)
)

print("STAGE 4 PASS", flush=True)

print()
print("R27F2 COMPLETE")
print("CORRECTED CONSTRAINT IDENTITIES: PASS")
print("SINGULAR-SURFACE REDUCTION: PASS")
print("MULTIPLIER DEPENDENCE: PASS")
print("B-HC BRACKET: OPEN")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
