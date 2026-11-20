"""CP098-T87A: Independent geometric derivative audit.

Exact non-adapted points; all coordinate directions.
No momentum constraints and no compact-bracket assumptions.
"""
import sympy as sp

print("=== CP098-T87A ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)
eps = sp.symbols("eps", real=True)
basis = sp.eye(4)

def dot(x, y):
    return (x.T * eta * y)[0]

def geom(u, q):
    v = sp.sqrt(-dot(u, u))
    T = u / v
    tq = dot(T, q)
    w = q + T * tq
    L = sp.sqrt(dot(w, w))
    N = w / L
    k = L / v**2
    return v, T, L, N, k

def simplify(x):
    if isinstance(x, sp.MatrixBase):
        return x.applyfunc(sp.simplify)
    return sp.simplify(x)

def iszero(x):
    z = simplify(x)
    if isinstance(z, sp.MatrixBase):
        return all(t == 0 for t in z)
    return z == 0

def direct(u, q, h, variable):
    if variable == "u":
        ue, qe = u + eps*h, q
    else:
        ue, qe = u, q + eps*h

    return [
        simplify(sp.diff(obj, eps).subs(eps, 0))
        for obj in geom(ue, qe)
    ]

def predicted(u, q, h, variable):
    v, T, L, N, k = geom(u, q)
    tq = dot(T, q)

    if variable == "q":
        dv = sp.S.Zero
        dT = sp.zeros(4, 1)
        dL = dot(N, h)
        dN = (
            h + T*dot(T, h) - N*dot(N, h)
        ) / L
        dk = dL / v**2

    else:
        dv = -dot(T, h)
        dT = (h + T*dot(T, h)) / v

        # q is held fixed in this variation.
        dw = dT*tq + T*dot(dT, q)

        dL = dot(N, dw)
        dN = (dw - N*dL) / L
        dk = dL/v**2 - 2*L*dv/v**3

    return [dv, dT, dL, dN, dk]

cases = [
    (
        sp.Matrix([5, 3, 0, 0]),
        sp.Matrix([3, 5, 0, 0])
    ),
    (
        sp.Matrix([10, 2, 3, 1]),
        sp.Matrix([4, 1, 5, 2])
    ),
]

names = ("v", "T", "L", "N", "kappa")
total = 0

for case_id, (u, q) in enumerate(cases, 1):
    v, T, L, N, k = geom(u, q)

    assert sp.simplify(v**2 + dot(u, u)) == 0
    assert sp.simplify(dot(T, T) + 1) == 0
    assert sp.simplify(dot(T, N)) == 0
    assert sp.simplify(dot(N, N) - 1) == 0
    assert v.is_positive and L.is_positive

    print(
        f"CASE_{case_id}: regular geometry verified",
        flush=True
    )

    for variable in ("u", "q"):
        for i in range(4):
            h = basis[:, i]

            actual = direct(u, q, h, variable)
            expected = predicted(u, q, h, variable)

            for name, a, b in zip(names, actual, expected):
                residual = simplify(a - b)

                assert iszero(residual), (
                    f"STOP: case={case_id}, "
                    f"variable={variable}, "
                    f"direction={i}, "
                    f"quantity={name}, "
                    f"residual={residual}"
                )

                total += 1

        print(
            f"PASS: case {case_id}, all four "
            f"{variable}-directions",
            flush=True
        )

print("EXACT_COMPONENT_IDENTITIES =", total)
assert total == 80

print("PASS: independent geometric derivative checks")
print("ARBITRARY-POINT FORMAL IDENTITY: OPEN")
print("FULL CANONICAL BRACKET PROOF: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T87A COMPLETE")
