import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R25C ===", flush=True)

lam = S.symbols("lam", real=True)
ia, is_, iw = S.symbols("ia is iw", real=True)

hA, hS, hW = S.symbols(
    "hA hS hW", positive=True
)
kA, kS, kW = S.symbols(
    "kA kS kW", real=True
)

eta = S.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

# First previously verified worldline jet.
u = S.Matrix([5, 3, 0, 0])
q = S.Matrix([1, 2, 2, 1])
j = S.Matrix([4, 1, 3, 6])

ss = -dot(u, u)
v = S.sqrt(ss)
vprime = -dot(u, q)/v

d = dot(u, q)
Q = dot(q, q)
z = dot(u, j)
h = dot(q, j)
J = dot(j, j)

W = S.factor(Q + d**2/ss)
Z = S.factor(
    J + z**2/ss
    - (h + d*z/ss)**2/W
)

assert ss > 0 and W > 0

print("STAGE 1: geometric quantities", flush=True)
print("v =", v)
print("vprime =", vprime)
print("W =", W)
print("Z =", Z)

# R25B's explicitly constructed coupling functions.
Afield = (3 + 2*ia)*is_
Cfield = (2 + ia + S.Rational(3,2)*ia**2)*iw

Lgeom = (
    -Afield*W/v**3
    -Cfield*Z/(v*W)
)

coords = (ia, is_, iw)
hs = (hA, hS, hW)
ks = (kA, kS, kW)

# The R25B witness:
# IA=lambda, IS=IW=1.
trajectories = (
    lam,
    S.Integer(1),
    S.Integer(1)
)

print("STAGE 2: identity Euler-Lagrange residuals",
      flush=True)

sources = [
    S.diff(Lgeom, x)
    for x in coords
]

# At lambda=0, compute each EL residual.
residuals = []

for n in range(3):
    I = trajectories[n]
    I0 = I.subs(lam, 0)
    I1 = S.diff(I, lam).subs(lam, 0)
    I2 = S.diff(I, lam, 2).subs(lam, 0)

    source = sources[n].subs({
        ia: 0, is_: 1, iw: 1
    })

    residual = S.factor(
        hs[n]*(I2/v - I1*vprime/v**2)
        + v*ks[n]*I0
        - source
    )
    residuals.append(residual)

    print(
        f"EL RESIDUAL {coords[n]} =",
        residual,
        flush=True
    )

print("STAGE 2 COMPLETE", flush=True)

print("STAGE 3: illustrative unit coefficients",
      flush=True)

unit_values = {
    hA: 1, hS: 1, hW: 1,
    kA: 1, kS: 1, kW: 1
}

unit_residuals = [
    S.factor(r.subs(unit_values))
    for r in residuals
]

print("UNIT-COEFFICIENT RESIDUALS =",
      unit_residuals, flush=True)

unit_solution = all(r == 0 for r in unit_residuals)

print("UNIT WITNESS SATISFIES EL EQUATIONS =",
      unit_solution, flush=True)

print("STAGE 4: pointwise coefficient conditions",
      flush=True)

# Solve only the three equations at this jet.
# This is NOT a solution along the full worldline.
conditions = S.solve(
    residuals,
    (hA, kS, kW),
    dict=True
)

print("POINTWISE SOLUTIONS =", conditions, flush=True)

for solution in conditions:
    check = [
        S.simplify(r.subs(solution))
        for r in residuals
    ]
    assert check == [0, 0, 0], (
        "STOP: coefficient solution failed substitution"
    )

print("POINTWISE SUBSTITUTION: PASS", flush=True)

print()
print("R25C COMPLETE")
print("IDENTITY EL RESIDUALS: DERIVED AT JET 1")
print("UNIT-COEFFICIENT WITNESS:", 
      "PASS" if unit_solution else "NOT A SOLUTION")
print("POINTWISE COEFFICIENT CONDITIONS: COMPUTED")
print("FULL TRAJECTORY SOLUTION: OPEN")
print("WORLDLINE EQUATIONS OF MOTION: OPEN")
print("HAMILTONIAN RECONSTRUCTION: OPEN")
print("DIRAC SYSTEM: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("ELAPSED SECONDS =",
      round(monotonic()-start, 3))
