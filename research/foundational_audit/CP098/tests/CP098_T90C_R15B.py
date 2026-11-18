import sympy as s

print("=== CP098-T90C-R15B ===", flush=True)

eta = s.diag(-1, 1, 1, 1)

def dot(a, b):
    return (a.T * eta * b)[0]

A, B, C = s.symbols(
    "A B C", real=True, nonzero=True
)

jets = [
    ([5, 3, 0, 0], [1, 2, 2, 1], [4, 1, 3, 6]),
    ([13, 3, 4, 0], [2, 3, 1, 4], [1, 2, 5, 7])
]

for index, (uv, qv, jv) in enumerate(jets, 1):
    print("JET", index, flush=True)

    u = s.Matrix(uv)
    v = s.sqrt(-dot(u, u))
    assert v.is_positive

    T = u/v

    # Differentiate curvature with respect to q,
    # holding the timelike velocity u fixed.
    Q = s.symbols("q0:4")
    q = s.Matrix(Q)

    w = q + dot(T, q)*T
    k2 = s.factor(dot(w, w)/v**4)

    Lcurv = -v*A*k2

    gradient = s.Matrix([
        s.diff(Lcurv, x) for x in Q
    ])

    qsubs = dict(zip(Q, qv))
    w0 = w.subs(qsubs)

    k = s.sqrt(k2.subs(qsubs))
    assert k.is_positive

    N = w0/(v**2*k)

    expected_gradient = -2*A*k/v*eta*N

    curvature_residual = (
        gradient.subs(qsubs)-expected_gradient
    ).applyfunc(s.simplify)

    assert curvature_residual == s.zeros(4, 1)

    # Construct the exact transverse projector
    # at this fixed regular worldline jet.
    projector = (
        s.eye(4)
        + T*(T.T*eta)
        - N*(N.T*eta)
    )

    assert projector.rank() == 2

    # Differentiate torsion with respect to j,
    # holding u and q fixed.
    J = s.symbols("j0:4")
    j = s.Matrix(J)

    transverse_j = projector*j

    Ltors = (
        -C*dot(transverse_j, transverse_j)
        /(v**5*k**2)
    )

    P = s.Matrix([
        s.diff(Ltors, x) for x in J
    ])

    jsubs = dict(zip(J, jv))
    j0 = s.Matrix(jv)
    P0 = P.subs(jsubs)

    P2 = (P0.T*eta*P0)[0]

    # Reduced candidate: B already contains
    # the Legendre-transformed identity energy.
    Lreduced = -v*(B+A*k**2)+Ltors.subs(jsubs)

    H = (P0.T*j0)[0]-Lreduced

    expected_H = (
        v*(B+A*k**2)
        -v**5*k**2*P2/(4*C)
    )

    legendre_residual = s.simplify(H-expected_H)

    assert legendre_residual == 0

    print("v =", v)
    print("kappa =", k)
    print("PROJECTOR RANK =", projector.rank())
    print("CURVATURE RESIDUAL =", curvature_residual)
    print("LEGENDRE RESIDUAL =", legendre_residual)
    print("PASS", flush=True)

print()
print("R15B COMPLETE")
print("LIMIT: reduced candidate Legendre calculation")
print("FULL OSTROGRADSKY DERIVATION: OPEN")
print("DEFI TORSION IDENTIFICATION: OPEN")
