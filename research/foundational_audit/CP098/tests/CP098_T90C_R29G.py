import sympy as S
import runpy
import contextlib
import io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R29G ===", flush=True)

folder = Path(__file__).parent

with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path(str(folder / "CP098_T90C_R27U.py"))

u, q, P, p0, p1 = (
    ns[k] for k in ("u", "q", "P", "p0", "p1")
)
H, phi1, phi2 = (
    ns[k] for k in ("H", "phi1", "phi2")
)
IA, IS, IW, pA, pS, pW = (
    ns[k] for k in ("IA", "IS", "IW", "pA", "pS", "pW")
)

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [(IA, pA), (IS, pS), (IW, pW)]
)

eta = [-1, 1, 1, 1]
s = -sum(eta[i]*u[i]**2 for i in range(4))
v = S.sqrt(s)
d = sum(eta[i]*u[i]*q[i] for i in range(4))
W = sum(eta[i]*q[i]**2 for i in range(4)) + d**2/s
PP = sum(eta[i]*P[i]**2 for i in range(4))

Q = (
    S.Rational(3, 2)*v
    + 56*v*W*PP/(961*IW)
    - 2*IS*W/v**3
)

wval = S.Integer(1)
rval = S.Integer(2)

Dbase = (
    802816*rval**4*wval**4
    + 476656*rval**2*wval**6
    + 41330688*rval**2*wval**2
    + 19393941*wval**4
    + 531948096
)/(14776336*wval**2)

z2min = S.factor(
    961*Dbase/(8*wval**2*(32*rval**2-31))
)

zval = S.sqrt(z2min + 1)

disc = S.factor(
    -wval**2*Dbase
    + 8*wval**4*(32*rval**2-31)*zval**2/961
)

assert S.simplify(disc - S.Rational(776, 961)) == 0

center = -S.Rational(32, 31)*rval*zval

branches = {
    "PLUS": center + 2*S.sqrt(disc)/wval**2,
    "MINUS": center - 2*S.sqrt(disc)/wval**2
}

def make_point(bval):
    ISv = 12/wval**2 + S.Rational(448, 961)*rval**2
    Mv = 1 + (
        S.Rational(9, 4) + ISv**2 + 1 + bval**2
    )/2

    point = {}

    for vec, values in (
        (u, [2, 0, 0, 0]),
        (q, [0, wval, 0, 0]),
        (P, [0, 0, rval, 0]),
        (p1, [0, 0, zval, 0]),
        (p0, [
            2*wval**2*rval**2/31 - Mv,
            -8*wval*rval*zval/31,
            0, 0
        ])
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
        pW: bval
    })

    return point

print("STAGE 1: canonical derivatives", flush=True)

derivatives = {}

for label, expr in (
    ("Q", Q),
    ("H", H),
    ("phi1", phi1),
    ("phi2", phi2)
):
    derivatives[label] = [
        (
            S.diff(expr, x),
            S.diff(expr, p)
        )
        for x, p in pairs
    ]
    print("PASS:", label, "derivatives", flush=True)

def bracket_at(label1, label2, point):
    result = S.Integer(0)

    for (fx, fp), (gx, gp) in zip(
        derivatives[label1],
        derivatives[label2]
    ):
        result += (
            fx.subs(point)*gp.subs(point)
            - fp.subs(point)*gx.subs(point)
        )

    return S.simplify(result)

print("STAGE 2: exact branch evaluation", flush=True)

for name, bval in branches.items():
    print("BRANCH:", name, flush=True)

    point = make_point(bval)

    assert S.simplify(Q.subs(point)) == 0
    print("PASS: Q=0", flush=True)

    F = bracket_at("Q", "H", point)
    K1 = bracket_at("Q", "phi1", point)
    K2 = bracket_at("Q", "phi2", point)

    print("F_Q =", F, flush=True)
    print("K_Q1 =", K1, flush=True)
    print("K_Q2 =", K2, flush=True)

    if K1 == 0 and K2 != 0:
        lam2 = S.simplify(-F/K2)
        assert S.simplify(F + lam2*K2) == 0
        print("PASS: Qdot fixes lambda2", flush=True)
        print("lambda2 =", lam2, flush=True)
    elif K1 == 0 and K2 == 0:
        if F == 0:
            print("DEGENERATE: Qdot imposes no multiplier condition")
        else:
            print("STOP: Qdot cannot vanish on this branch")
    else:
        print("OPEN: coupled lambda1/lambda2 equation")

print("R29G EVALUATION COMPLETE")
print("Ddot: NOT YET TESTED")
print("FULL DIRAC CLOSURE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =", round(perf_counter()-start, 3))
