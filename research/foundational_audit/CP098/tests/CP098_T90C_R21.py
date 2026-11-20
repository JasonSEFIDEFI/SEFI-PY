import sympy as s

print("=== CP098-T90C-R21 ===", flush=True)

eta = s.diag(-1, 1, 1, 1)
C = s.symbols("C", nonzero=True, real=True)

def dot(x, y):
    return (x.T * eta * y)[0]

def clean(x):
    return s.simplify(x)

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

    proj = (
        s.eye(4)
        + T*(T.T*eta)
        - N*(N.T*eta)
    )

    j_symbols = s.symbols("j0:4")
    j = s.Matrix(j_symbols)

    R = proj*j
    P = -2*C*eta*R/(v**5*k**2)

    W = P.jacobian(j)
    expected = -2*C*eta*proj/(v**5*k**2)

    assert (W-expected).applyfunc(clean) == s.zeros(4)
    assert W.rank() == 2

    assert (W*T).applyfunc(clean) == s.zeros(4,1)
    assert (W*N).applyfunc(clean) == s.zeros(4,1)

    # The fifth derivative enters P'' through
    # j'' = b. Its coefficient is W.
    b = s.Matrix([2,-1,4,3])
    response = (W*b).applyfunc(clean)

    # Response is covariant; raise its index
    # before applying the vector projector.
    raised = eta*response

    assert (
        proj*raised-raised
    ).applyfunc(clean) == s.zeros(4,1)

    assert clean((response.T*T)[0]) == 0
    assert clean((response.T*N)[0]) == 0

    print()
    print("JET", index)
    print("v =", v)
    print("kappa =", k)
    print("HESSIAN RANK =", W.rank())
    print("FIFTH-DERIVATIVE RESPONSE =", response)
    print("TANGENT RESPONSE =", clean((response.T*T)[0]))
    print("NORMAL RESPONSE =", clean((response.T*N)[0]))
    print("PASS", flush=True)

print()
print("R21 COMPLETE")
print("FIFTH-DERIVATIVE COEFFICIENT: VERIFIED")
print("FULL p0: NOT YET DERIVED")
print("FULL DIRAC SYSTEM: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
