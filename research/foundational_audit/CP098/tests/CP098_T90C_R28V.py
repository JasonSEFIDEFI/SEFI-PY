
import sympy as S
import runpy, contextlib, io
from pathlib import Path
from time import perf_counter

start = perf_counter()
print("=== CP098-T90C-R28V ===", flush=True)

folder = Path(__file__).parent

with contextlib.redirect_stdout(io.StringIO()):
    old = runpy.run_path(str(folder / "CP098_T90C_R28K2.py"))
    new = runpy.run_path(str(folder / "CP098_T90C_R28Q.py"))

base = old["ns"]
B, H, phi1 = base["B"], base["H"], base["phi1"]

u, q, P = base["u"], base["q"], base["P"]
p0, p1 = base["p0"], base["p1"]

pairs = (
    list(zip(u, p1))
    + list(zip(q, P))
    + [
        (base["IA"], base["pA"]),
        (base["IS"], base["pS"]),
        (base["IW"], base["pW"])
    ]
)

assert len(pairs) == 11

def bracket(f, g):
    return sum(
        (
            S.diff(f, x)*S.diff(g, p)
            - S.diff(f, p)*S.diff(g, x)
        )
        for x, p in pairs
    )

print("STAGE 1: symbolic inner brackets", flush=True)

B_phi1 = bracket(B, phi1)
H_phi1 = bracket(H, phi1)

print("PASS: inner brackets constructed", flush=True)

# Jacobi RHS, evaluated by differentiating the
# unrestricted inner expressions before substitution.
for label, witness in (
    ("FIRST", old),
    ("SECOND", new)
):
    print("STAGE 2:", label, "witness", flush=True)

    at = witness["at"]

    left = S.Integer(0)
    right = S.Integer(0)

    for index, (x, p) in enumerate(pairs, 1):
        # Direct K_D1 via product rule, without
        # constructing the full expression for D.
        direction_x = at(S.diff(phi1, p))
        direction_p = -at(S.diff(phi1, x))

        # Direct evaluation of the Jacobi RHS.
        right += (
            at(S.diff(B_phi1, x))*at(S.diff(H, p))
            - at(S.diff(B_phi1, p))*at(S.diff(H, x))
            + at(S.diff(B, x))*at(S.diff(H_phi1, p))
            - at(S.diff(B, p))*at(S.diff(H_phi1, x))
        )

        print(label, "pair", index, "of 11", flush=True)

    right = S.simplify(right)

    # Earlier full second-derivative calculations:
    # R28M and R28S both established K_D1 = 0.
    expected = S.Integer(0)

    print(label, "Jacobi RHS =", right, flush=True)
    assert right == expected

    print("PASS:", label, "K_D1=0", flush=True)

print("R28V PASS: JACOBI RHS AT TWO WITNESSES")
print("FAMILY-WIDE IDENTITY: NOT YET ESTABLISHED")
print("FULL DIRAC CLOSURE: OPEN")
print(
    "ELAPSED SECONDS =",
    round(perf_counter()-start, 3),
    flush=True
)
