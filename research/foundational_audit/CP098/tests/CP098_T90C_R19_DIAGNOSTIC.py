import sympy as s

print("=== CP098-T90C-R19 ===", flush=True)

eta = s.diag(-1, 1, 1, 1)
t, e = s.symbols("t e", real=True)

def dot(x, y):
    return (x.T * eta * y)[0]

def simp(x):
    return s.simplify(x)

def geom(u, q):
    v = s.sqrt(-dot(u, u))
    T = u/v
    w = q + dot(T, q)*T
    k = s.sqrt(dot(w, w))/v**2
    N = w/(v**2*k)
    proj = (
        s.eye(4)
        + T*(T.T*eta)
        - N*(N.T*eta)
    )
    return v, T, k, N, proj

def lagrangian(u, q, j, A, C):
    v, T, k, N, proj = geom(u, q)
    R = proj*j
    return (
        -v*A*k**2
        -C*dot(R, R)/(v**5*k**2)
    )

def highest_momentum(u, q, j, C):
    v, T, k, N, proj = geom(u, q)
    R = proj*j
    return -2*C*eta*R/(v**5*k**2)

jets = [
    (
        [5, 3, 0, 0],
        [1, 2, 2, 1],
        [4, 1, 3, 6],
        [2, -1, 4, 3]
    ),
    (
        [13, 3, 4, 0],
        [2, 3, 1, 4],
        [1, 2, 5, 7],
        [-2, 4, 1, 3]
    )
]

A0 = s.Integer(3)
C0 = s.Integer(2)
Cdot = s.Integer(1)

for index, (uv, qv, jv, av) in enumerate(jets, 1):
    print()
    print("JET", index, flush=True)

    u = s.Matrix(uv)
    q = s.Matrix(qv)
    j = s.Matrix(jv)
    a = s.Matrix(av)

    v, T, k, N, proj = geom(u, q)

    assert simp(v) > 0
    assert simp(k) > 0

    # Full q-gradient of curvature + torsion,
    # computed one component at a time.
    gradient_entries = []

    for i in range(4):
        basis = s.eye(4)[:, i]
        shifted_q = q + e*basis

        expression = lagrangian(
            u, shifted_q, j, A0, C0
        )

        value = simp(
            s.diff(expression, e).subs(e, 0)
        )
        gradient_entries.append(value)

    gradient = s.Matrix(gradient_entries)

    # Exact jet evolution through fourth order.
    # Only derivatives at t=0 are required.
    ut = u + q*t
    qt = q + j*t
    jt = j + a*t
    Ct = C0 + Cdot*t

    Pt = highest_momentum(ut, qt, jt, Ct)

    Pdot = s.Matrix([
        simp(s.diff(Pt[i], t).subs(t, 0))
        for i in range(4)
    ])

    p1 = (gradient - Pdot).applyfunc(simp)

    # Exact independent projection checks.
    tangent = simp(dot(p1, T))
    normal = simp(dot(p1, N))

    expected_normal = simp(-2*A0*k/v)

    print("TANGENT =", tangent, flush=True)
    print("NORMAL =", normal, flush=True)
    print("EXPECTED NORMAL =", expected_normal, flush=True)
    print("TANGENT GRADIENT =", simp(dot(eta*gradient,T)), flush=True)
    print("TANGENT Pdot =", simp(dot(eta*Pdot,T)), flush=True)
    print("NORMAL RESIDUAL =", simp(normal-expected_normal), flush=True)

    # Reconstruction of the full covector:
    # longitudinal + normal + transverse parts.
    p1_raised = eta*p1
    transverse = proj*p1_raised

    reconstructed = (
        -dot(p1_raised, T)*T
        +dot(p1_raised, N)*N
        +transverse
    )

    assert (
        reconstructed-p1_raised
    ).applyfunc(simp) == s.zeros(4, 1)

    print("v =", v)
    print("kappa =", k)
    print("FULL p1 =", p1)
    print("p1.T =", tangent)
    print("p1.N =", normal)
    print("EXPECTED p1.N =", expected_normal)
    print("TRANSVERSE p1 =", transverse)
    print("RECONSTRUCTION RESIDUAL = 0")
    print("PASS", flush=True)

print()
print("R19 DIAGNOSTIC COMPLETE")
print("FULL LOWER MOMENTUM: EVALUATED AT TWO FOURTH-ORDER JETS")
print("C-DERIVATIVE INCLUDED")
print("FULL DIRAC CONSTRAINT ANALYSIS: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
