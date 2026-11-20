import sympy as S
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R27F ===", flush=True)

v, W, C = S.symbols("v W C", positive=True)
A, d, M, J = S.symbols("A d M J", real=True)
e, f, g, P2, r = S.symbols(
    "e f g P2 r", real=True
)
phi1, phi2 = S.symbols("phi1 phi2", real=True)
lam1, lam2, D = S.symbols("lam1 lam2 D")

# r = p1.q
# g = p0.u
# e = p0.q
# f = p1.eta.P
# D = {B,Hc}, not yet evaluated.

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

# Exact multiplier coefficients from R27E.
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
    result = S.cancel(expression)
    if result != 0:
        raise AssertionError(
            f"{name}: residual = {result}"
        )
    print("PASS:", name, flush=True)

print("STAGE 1: exact constraint identities", flush=True)

check(
    "K1 relative to theta and psi",
    K1-(theta-psi-6*b)
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

print("STAGE 2: A=0 restriction", flush=True)

check(
    "K1 at A=0",
    (K1-theta+psi).subs(A,0)
)

check(
    "K2 at A=0",
    (K2-B-2*W*J/v**2).subs(A,0)
)

print("STAGE 2 PASS", flush=True)

print("STAGE 3: five-constraint reduction", flush=True)

# On A=0, psi=0 implies r=0.
# theta=0 then determines g.
# B=0 determines e.
# These substitutions are derived, not guessed.

g_solution = S.solve(
    theta.subs({A:0,r:0}),g
)
assert len(g_solution)==1

e_solution = S.solve(
    B.subs(A,0),e
)
assert len(e_solution)==1

surface = {
    A:0,
    r:0,
    g:g_solution[0],
    e:e_solution[0]
}

check(
    "K1 vanishes on constraint surface",
    K1.subs(surface)
)

check(
    "K2 retains identity-sector term",
    K2.subs(surface)-2*W*J/v**2
)

print("STAGE 3 PASS", flush=True)

print("STAGE 4: B preservation equation", flush=True)

Bdot = D+lam1*K1+lam2*K2

check(
    "restricted Bdot",
    Bdot.subs(surface)
    -(D+lam2*2*W*J/v**2)
)

print("STAGE 4 PASS", flush=True)

print()
print("R27F COMPLETE")
print("FULL CONSTRAINT REDUCTION: VERIFIED")
print("K1 ON SINGULAR SURFACE: ZERO")
print("K2 ON SINGULAR SURFACE: 2*W*J/v**2")
print("B-HC BRACKET: OPEN")
print("DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start,3)
)
