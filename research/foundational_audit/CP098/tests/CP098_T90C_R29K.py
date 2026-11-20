import sympy as S
import runpy
import contextlib
import io

from pathlib import Path
from time import perf_counter

start = perf_counter()

print("=== CP098-T90C-R29K ===", flush=True)

folder = Path(__file__).parent

# ============================================================
# STAGE 0: Load canonical source
# ============================================================

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(
        str(folder / "CP098_T90C_R27U.py")
    )

u, q, P, p0, p1 = (
    ns[k] for k in ("u", "q", "P", "p0", "p1")
)

B = ns["B"]
H = ns["H"]
phi2 = ns["phi2"]

IA, IS, IW, pA, pS, pW = (
    ns[k]
    for k in ("IA", "IS", "IW", "pA", "pS", "pW")
)

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [(IA, pA), (IS, pS), (IW, pW)]
)

assert len(pairs) == 11

# ============================================================
# STAGE 1: Define independent family parameters
# ============================================================

rho, z, b = S.symbols(
    "rho z b",
    real=True
)

w = S.Integer(1)

ISv = (
    12
    + S.Rational(448, 961) * rho**2
)

Mv = (
    1
    + (
        S.Rational(9, 4)
        + ISv**2
        + 1
        + b**2
    ) / 2
)

point = {}

for vec, values in (
    (u, [2, 0, 0, 0]),
    (q, [0, 1, 0, 0]),
    (P, [0, 0, rho, 0]),
    (p1, [0, 0, z, 0]),
    (
        p0,
        [
            2 * rho**2 / 31 - Mv,
            -8 * rho * z / 31,
            0,
            0
        ]
    )
):
    point.update({
        sym: S.sympify(value)
        for sym, value in zip(vec, values)
    })

point.update({
    IA: S.Rational(-3, 2),
    IS: ISv,
    IW: S.Integer(1),
    pA: S.Integer(0),
    pS: S.Integer(0),
    pW: b
})

def at(expr):
    return S.cancel(expr.subs(point))

print(
    "INDEPENDENT PARAMETERS: rho, z, b",
    flush=True
)

print(
    "FIXED: w=1, t=0, a=0",
    flush=True
)

# ============================================================
# STAGE 2: Construct canonical derivatives
# ============================================================

print(
    "STAGE 2: canonical first derivatives",
    flush=True
)

first = []

for x, p in pairs:
    first.append(
        (
            S.diff(B, x),
            S.diff(B, p),
            S.diff(H, x),
            S.diff(H, p)
        )
    )

print(
    "PASS: 11 canonical derivative groups",
    flush=True
)

# ============================================================
# STAGE 3: Construct Hamiltonian vector field
# ============================================================

print(
    "STAGE 3: Hamiltonian vector field",
    flush=True
)

directions = {}

for x, p in pairs:
    directions[x] = at(S.diff(H, p))
    directions[p] = -at(S.diff(H, x))

active = [
    (var, value)
    for var, value in directions.items()
    if value != 0
]

print(
    "Active Hamiltonian directions:",
    len(active),
    flush=True
)

def directional(expr):
    result = S.Integer(0)

    for var, value in active:
        result += value * at(S.diff(expr, var))

    return S.cancel(result)

# ============================================================
# STAGE 4: Compute D = {B,H}
# ============================================================

print(
    "STAGE 4: first bracket D",
    flush=True
)

Dcalc = S.Integer(0)

for Bx, Bp, Hx, Hp in first:
    Dcalc += (
        at(Bx) * at(Hp)
        - at(Bp) * at(Hx)
    )

Dcalc = S.factor(Dcalc)

Dbase = (
    802816 * rho**4
    + 476656 * rho**2
    + 41330688 * rho**2
    + 19393941
    + 531948096
) / S.Integer(14776336)

Dexpected = (
    Dbase
    + 8 * z**2 / 31
    + b**2 / 4
    + 16 * b * rho * z / 31
)

assert S.simplify(
    Dcalc - Dexpected
) == 0

print(
    "PASS: D matches R28Q general-b formula",
    flush=True
)

# ============================================================
# STAGE 5: Compute F_D = {{B,H},H}
# ============================================================

print(
    "STAGE 5: compute F_D",
    flush=True
)

F = S.Integer(0)

for i, (Bx, Bp, Hx, Hp) in enumerate(
    first, 1
):
    contribution = (
        directional(Bx) * at(Hp)
        + at(Bx) * directional(Hp)
        - directional(Bp) * at(Hx)
        - at(Bp) * directional(Hx)
    )

    F = S.cancel(F + contribution)

    print(
        "F_D pair",
        i,
        "of 11",
        flush=True
    )

F = S.factor(F)

print(
    "PASS: F_D constructed",
    flush=True
)

# ============================================================
# STAGE 6: Compute K_D2 via Jacobi identity
# ============================================================

print(
    "STAGE 6: compute K_D2 via Jacobi",
    flush=True
)

# phi2 generates:
# delta(q_i) = q_i
# delta(P_i) = -P_i
#
# All other canonical variables are unchanged.

Bphi2 = sum(
    S.diff(B, q[i]) * q[i]
    - S.diff(B, P[i]) * P[i]
    for i in range(4)
)

Hphi2 = sum(
    S.diff(H, q[i]) * q[i]
    - S.diff(H, P[i]) * P[i]
    for i in range(4)
)

def bracket_at(f, g, label):
    result = S.Integer(0)

    for i, (x, p) in enumerate(pairs, 1):
        result += (
            at(S.diff(f, x))
            * at(S.diff(g, p))
            - at(S.diff(f, p))
            * at(S.diff(g, x))
        )

        print(
            label,
            "pair",
            i,
            "of 11",
            flush=True
        )

    return S.factor(result)

T1 = bracket_at(
    Bphi2,
    H,
    "T1"
)

T2 = bracket_at(
    B,
    Hphi2,
    "T2"
)

K = S.factor(T1 + T2)

print(
    "PASS: K_D2 constructed",
    flush=True
)

# ============================================================
# STAGE 7: General reduced multiplier
# ============================================================

print(
    "STAGE 7: reduced multiplier",
    flush=True
)

lam2 = (
    -112 * rho * (b * rho + z)
    / (112 * rho**2 + 2883)
)

assert (
    112 * rho**2 + 2883
).is_positive is True

print(
    "PASS: multiplier denominator positive",
    flush=True
)

# ============================================================
# STAGE 8: Combine Ddot
# ============================================================

print(
    "STAGE 8: construct Ddot",
    flush=True
)

Ddot = S.factor(
    S.cancel(F + lam2 * K)
)

print(
    "PASS: Ddot constructed",
    flush=True
)

# ============================================================
# STAGE 9: Print exact results
# ============================================================

print(
    "\n===== EXACT GENERAL FAMILY RESULTS =====",
    flush=True
)

print(
    "D =",
    S.factor(Dcalc),
    flush=True
)

print(
    "F_D =",
    F,
    flush=True
)

print(
    "K_D2 =",
    K,
    flush=True
)

print(
    "lambda2 =",
    lam2,
    flush=True
)

print(
    "Ddot =",
    Ddot,
    flush=True
)

# ============================================================
# STAGE 10: R29H witness regressions
# ============================================================

print(
    "\nSTAGE 10: R29H regressions",
    flush=True
)

rv = S.Integer(2)

Dbase_v = S.factor(
    Dbase.subs(rho, rv)
)

z2min = S.factor(
    961 * Dbase_v
    / (8 * (32 * rv**2 - 31))
)

zv = S.sqrt(z2min + 1)

disc = S.Rational(776, 961)

center = (
    -S.Rational(32, 31)
    * rv * zv
)

branches = {
    "PLUS": center + 2*S.sqrt(disc),
    "MINUS": center - 2*S.sqrt(disc)
}

expected_dot = {
    "PLUS": (
        -S.Rational(
            2196291909748885,
            111262356920464396
        ) * S.sqrt(144209559530)
        + S.Rational(
            1152303896300159,
            9250279092157
        ) * S.sqrt(194)
    ),
    "MINUS": (
        -S.Rational(
            2196291909748885,
            111262356920464396
        ) * S.sqrt(144209559530)
        - S.Rational(
            1152303896300159,
            9250279092157
        ) * S.sqrt(194)
    )
}

for name, bv in branches.items():
    subs = {
        rho: rv,
        z: zv,
        b: bv
    }

    assert S.simplify(
        Dcalc.subs(subs)
    ) == 0

    actual = S.simplify(
        Ddot.subs(subs)
    )

    assert S.simplify(
        actual - expected_dot[name]
    ) == 0

    print(
        "PASS:",
        name,
        "R29H exact Ddot regression",
        flush=True
    )

print("\nR29K PASS", flush=True)

print(
    "GENERAL w=1 FAMILY Ddot: DERIVED",
    flush=True
)

print(
    "SIMULTANEOUS D=0 AND Ddot=0: OPEN",
    flush=True
)

print(
    "GENERAL-w EXTENSION: OPEN",
    flush=True
)

print(
    "FULL DIRAC CLOSURE: OPEN",
    flush=True
)

print(
    "PHYSICAL VALIDATION: NOT CLAIMED",
    flush=True
)

print(
    "ELAPSED SECONDS =",
    round(perf_counter() - start, 3),
    flush=True
)
