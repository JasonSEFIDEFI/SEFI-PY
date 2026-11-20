import sympy as S
from time import monotonic

start = monotonic()
print("=== CP098-T90C-R26B ===", flush=True)

eta = S.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

jets = [
    (
        "JET 1",
        [5, 3, 0, 0],
        [1, 2, 2, 1],
        [4, 1, 3, 6],
    ),
    (
        "JET 2",
        [13, 3, 4, 0],
        [2, 3, 1, 4],
        [1, 4, 2, 5],
    ),
]

C = S.Integer(2)

for name, uu, qq, jj in jets:
    print()
    print(name, flush=True)

    u = S.Matrix(uu)
    q = S.Matrix(qq)
    j = S.Matrix(jj)

    s = -dot(u, u)
    assert s > 0

    v = S.sqrt(s)
    d = dot(u, q)

    # Component of q orthogonal to u.
    n = S.simplify(q + d*u/s)

    W = S.factor(dot(n, n))
    assert W > 0

    # Component of j orthogonal to u and n.
    R = S.simplify(
        j
        + dot(u, j)*u/s
        - dot(n, j)*n/W
    )

    Z = S.factor(dot(R, R))

    print("STAGE 1: projection geometry", flush=True)

    assert S.simplify(dot(u, n)) == 0
    assert S.simplify(dot(u, R)) == 0
    assert S.simplify(dot(q, R)) == 0

    print("v =", v)
    print("W =", W)
    print("Z =", Z)
    print("PROJECTION ORTHOGONALITY: PASS")

    print("STAGE 2: direct momentum derivative",
          flush=True)

    x = S.Matrix(S.symbols("j0:4", real=True))

    Rx = S.simplify(
        x
        + dot(u, x)*u/s
        - dot(n, x)*n/W
    )

    Zx = S.expand(dot(Rx, Rx))

    Lx = -C*Zx/(v*W)

    P_direct = S.Matrix([
        S.diff(Lx, x[i])
        for i in range(4)
    ])

    P_expected = -2*C*eta*Rx/(v*W)

    assert all(
        S.simplify(P_direct[i]-P_expected[i]) == 0
        for i in range(4)
    )

    P = S.simplify(P_direct.subs(
        dict(zip(x, j))
    ))

    print("HIGHEST MOMENTUM =", P)
    print("MOMENTUM DERIVATIVE: PASS")

    print("STAGE 3: primary constraints",
          flush=True)

    phi1 = S.factor((u.T*P)[0])
    phi2 = S.factor((q.T*P)[0])

    print("u^T P =", phi1)
    print("q^T P =", phi2)

    assert phi1 == 0
    assert phi2 == 0

    # These are independent as momentum constraints
    # because u and q are linearly independent.
    constraint_jacobian = S.Matrix([
        list(u),
        list(q),
    ])

    assert constraint_jacobian.rank() == 2

    print("CONSTRAINT INDEPENDENCE: PASS")

    print("STAGE 4: highest-derivative Hessian",
          flush=True)

    H = S.hessian(Lx, list(x))

    # Direct projected Hessian, independently
    # assembled from u and n.
    eu = eta*u
    en = eta*n

    B = (
        eta
        + eu*eu.T/s
        - en*en.T/W
    )

    H_expected = -2*C*B/(v*W)

    assert all(
        S.simplify(
            H[i, k]-H_expected[i, k]
        ) == 0
        for i in range(4)
        for k in range(4)
    )

    assert H.rank() == 2
    assert H*u == S.zeros(4, 1)
    assert H*q == S.zeros(4, 1)

    print("HESSIAN RANK =", H.rank())
    print("HESSIAN NULLITY =", 4-H.rank())
    print("NULL VECTORS: u, q")
    print("HESSIAN CROSS-CHECK: PASS")

print()
print("R26B COMPLETE")
print("HIGHEST MOMENTUM MAP: VERIFIED AT TWO JETS")
print("TWO INDEPENDENT PRIMARY CONSTRAINT CANDIDATES")
print("HESSIAN RANK: TWO AT BOTH JETS")
print("HESSIAN NULLITY: TWO AT BOTH JETS")
print("GLOBAL CONSTRAINT COMPLETENESS: OPEN")
print("POISSON BRACKET ALGEBRA: OPEN")
print("DIRAC CONSISTENCY ALGORITHM: OPEN")
print("DEFI CORRESPONDENCE: OPEN")
print("PHYSICAL VALIDATION: NOT CLAIMED")
print("ELAPSED SECONDS =",
      round(monotonic()-start, 3))
