import sympy as s

print("=== CP098-T90C-R18B ===", flush=True)

eta = s.diag(-1, 1, 1, 1)
C = s.symbols("C", real=True, nonzero=True)

def dot(x, y):
    return (x.T * eta * y)[0]

def clean(M):
    return M.applyfunc(s.simplify)

jets = [
    ([5,3,0,0], [1,2,2,1]),
    ([13,3,4,0], [2,3,1,4])
]

for index, (uv, qv) in enumerate(jets, 1):
    u = s.Matrix(uv)
    q = s.Matrix(qv)

    v = s.sqrt(-dot(u,u))
    T = u/v
    w = q + dot(T,q)*T
    k = s.sqrt(dot(w,w))/v**2
    N = w/(v**2*k)

    assert v.is_positive
    assert k.is_positive

    projector = (
        s.eye(4)
        + T*(T.T*eta)
        - N*(N.T*eta)
    )

    # The Hessian is obtained by differentiating
    # the candidate torsion momentum.
    J = s.symbols("j0:4")
    j = s.Matrix(J)

    R = projector*j

    Ltors = -C*dot(R,R)/(v**5*k**2)

    P = s.Matrix([
        s.diff(Ltors,x) for x in J
    ])

    W = clean(P.jacobian(J))

    # Ostrogradsky differentiation contributes
    # minus the highest-derivative Hessian
    # to the dependence on a = y''''.
    snap_coefficient = -W

    expected = clean(
        2*C*eta*projector/(v**5*k**2)
    )

    assert clean(
        snap_coefficient-expected
    ) == s.zeros(4)

    assert W.rank() == 2
    assert snap_coefficient.rank() == 2

    assert clean(W*T) == s.zeros(4,1)
    assert clean(W*N) == s.zeros(4,1)

    # Explicitly test a nonzero transverse
    # fourth-derivative perturbation.
    a = s.Matrix([2,-1,4,3])

    response = clean(snap_coefficient*a)

    assert clean(
        projector*eta*response-eta*response
    ) == s.zeros(4,1)

    assert s.simplify(
        (response.T*T)[0]
    ) == 0

    assert s.simplify(
        (response.T*N)[0]
    ) == 0

    print()
    print("JET", index)
    print("v =", v)
    print("kappa =", k)
    print("HESSIAN RANK =", W.rank())
    print("SNAP RESPONSE =", response)
    print("TANGENT RESPONSE =", (response.T*T)[0])
    print("NORMAL RESPONSE =", (response.T*N)[0])
    print("PASS", flush=True)

print()
print("R18B COMPLETE")
print("FOURTH-DERIVATIVE COEFFICIENT: VERIFIED")
print("FULL LOWER MOMENTUM: NOT YET DERIVED")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
