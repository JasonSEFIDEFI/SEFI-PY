import sympy as S
from time import monotonic

print("=== CP098-T90C-R24C ===", flush=True)
start = monotonic()

# Six independent Lorentz invariants.
s, d, Q, z, h, J, A, C = S.symbols(
    "s d Q z h J A C", real=True
)

W = Q + d**2/s
Z = J + z**2/s - (h + d*z/s)**2/W

L = -A*W/s**S.Rational(3, 2) - C*Z/(S.sqrt(s)*W)

variables = (s, d, Q, z, h, J)
grad = {x: S.diff(L, x) for x in variables}

eta = S.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

def calculate(u0, q0, j0, a0, b0):
    u, q, j, a, b = [
        S.Matrix(x) for x in (u0, q0, j0, a0, b0)
    ]

    values = {
        s: -dot(u, u),
        d: dot(u, q),
        Q: dot(q, q),
        z: dot(u, j),
        h: dot(q, j),
        J: dot(j, j),
        A: S.Integer(3),
        C: S.Integer(2)
    }

    # Exact first worldline derivatives of invariants.
    rates = {
        s: -2*dot(u, q),
        d: dot(q, q) + dot(u, j),
        Q: 2*dot(q, j),
        z: dot(q, j) + dot(u, a),
        h: dot(j, j) + dot(q, a),
        J: 2*dot(j, a),
        A: S.Integer(2),
        C: S.Integer(1)
    }

    if values[s] <= 0:
        raise ValueError("STOP: non-timelike u")
    if (W.subs(values)) <= 0:
        raise ValueError("STOP: nonpositive transverse norm")

    def at(expr):
        return S.simplify(expr.subs(values))

    def D(expr):
        return S.simplify(sum(
            S.diff(expr, x).subs(values)*rate
            for x, rate in rates.items()
        ))

    # Partial-u gradient, with all four covector components.
    gu = (
        -2*at(grad[s])*eta*u
        + at(grad[d])*eta*q
        + at(grad[z])*eta*j
    )

    # Total derivative of the partial-q gradient.
    gqdot = (
        D(grad[d])*eta*u
        + (at(grad[d]) + 2*D(grad[Q]))*eta*q
        + (2*at(grad[Q]) + D(grad[h]))*eta*j
        + at(grad[h])*eta*a
    )

    return gu.applyfunc(S.simplify), gqdot.applyfunc(S.simplify)

cases = [
    {
        "name": "JET 1",
        "u": [5, 3, 0, 0],
        "q": [1, 2, 2, 1],
        "j": [4, 1, 3, 6],
        "a": [2, -1, 4, 3],
        "b": [1, 3, -2, 5],
        "Pddot": [
            S.Rational(3250183453, 492307344),
            -S.Rational(18580716433, 1476922032),
            -S.Rational(248659083, 54700816),
            S.Rational(110869279, 492307344)
        ],
        "expected": [
            S.Rational(1797334885453, 504122720256),
            -S.Rational(7603294311121, 1512368160768),
            S.Rational(7357855225, 7876917504),
            S.Rational(20080631585, 15753835008)
        ],
        "expected_gu0": S.Rational(
            16903019833, 11723784192
        ),
        "expected_gqdot0": S.Rational(
            1128841411619, 252061360128
        )
    },
    {
        "name": "JET 2",
        "u": [13, 3, 4, 0],
        "q": [2, 3, 1, 4],
        "j": [1, 4, 2, 5],
        "a": [3, -2, 1, 4],
        "b": [2, 1, -3, 5],
        "Pddot": [
            S.Rational(
                3316234564794173, 17856129161341584
            ),
            -S.Rational(
                23174704243331285, 53568387484024752
            ),
            -S.Rational(
                1701604078088603, 5952043053780528
            ),
            S.Rational(
                3713100266877115, 17856129161341584
            )
        ],
        "expected": [
            S.Rational(
                170110180254221449591,
                1481058777158316343296
            ),
            -S.Rational(
                133849518421960573363,
                493686259052772114432
            ),
            -S.Rational(
                57948361657711656139,
                370264694289579085824
            ),
            S.Rational(
                84233369470633399,
                428547099872198016
            )
        ]
    }
]

for case in cases:
    print()
    print(case["name"], "START", flush=True)

    gu, gqdot = calculate(
        case["u"], case["q"], case["j"],
        case["a"], case["b"]
    )

    print("INVARIANT PARTIAL-u =", gu, flush=True)
    print("INVARIANT q-GRADIENT DERIVATIVE =",
          gqdot, flush=True)

    if case["name"] == "JET 1":
        assert gu[0] == case["expected_gu0"], (
            "STOP: R24A partial-u mismatch"
        )
        assert gqdot[0] == case["expected_gqdot0"], (
            "STOP: R24A q-gradient mismatch"
        )
        print("R24A GRADIENT CROSS-CHECK: PASS",
              flush=True)

    # P'' values are exact R24B results.
    # They are not independently recomputed here.
    Pddot = S.Matrix(case["Pddot"])

    p0 = (gu - gqdot + Pddot).applyfunc(S.simplify)
    expected = S.Matrix(case["expected"])
    residual = (p0 - expected).applyfunc(S.simplify)

    print("GEOMETRIC p0 =", p0, flush=True)
    print("FOUR-COMPONENT RESIDUAL =",
          residual, flush=True)

    assert residual == S.zeros(4, 1), (
        "STOP: geometric p0 disagreement"
    )

    print(case["name"], "PASS", flush=True)
    print("ELAPSED SECONDS =",
          round(monotonic() - start, 2), flush=True)

print()
print("R24C COMPLETE")
print("INVARIANT GRADIENT METHOD: VERIFIED AT TWO JETS")
print("GEOMETRIC p0: EIGHT EXACT COMPONENT MATCHES")
print("FULL SYMBOLIC p0 IDENTITY: NOT YET PROVEN")
print("FULL ACTION MOMENTUM: OPEN")
print("HAMILTONIAN RECONSTRUCTION: OPEN")
print("FULL DIRAC SYSTEM: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
