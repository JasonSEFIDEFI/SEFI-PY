import sympy as s

print("=== CP098-T90C-R16 ===", flush=True)

eta = s.diag(-1, 1, 1, 1)
A, C = s.symbols("A C", real=True, nonzero=True)
e = s.symbols("e", real=True)

def dot(a, b):
    return (a.T * eta * b)[0]

def clean(x):
    return s.simplify(x)

jets = [
    ([5,3,0,0], [1,2,2,1], [4,1,3,6]),
    ([13,3,4,0], [2,3,1,4], [1,2,5,7])
]

for index, (uv,qv,jv) in enumerate(jets,1):
    u = s.Matrix(uv)
    q = s.Matrix(qv)
    j = s.Matrix(jv)

    v = s.sqrt(-dot(u,u))
    assert v.is_positive

    T = u/v
    w = q + dot(T,q)*T
    L = s.sqrt(dot(w,w))
    k = L/v**2
    N = w/L

    assert k.is_positive

    projector = (
        s.eye(4)
        + T*(T.T*eta)
        - N*(N.T*eta)
    )

    R = projector*j
    R2 = clean(dot(R,R))

    # Highest-derivative momentum, raised.
    P = -2*C*R/(v**5*k**2)

    # Transverse part of dN/dlambda.
    Nprime_transverse = R/(v**2*k)

    # Differentiate the torsion Lagrangian
    # in the fixed normal direction.
    qe = q + e*N
    we = qe + dot(T,qe)*T
    Le = s.sqrt(dot(we,we))
    ke = Le/v**2
    Ne = we/Le

    proje = (
        s.eye(4)
        + T*(T.T*eta)
        - Ne*(Ne.T*eta)
    )

    Re = proje*j

    Ltors_e = (
        -C*dot(Re,Re)/(v**5*ke**2)
    )

    normal_gradient = clean(
        s.diff(Ltors_e,e).subs(e,0)
    )

    # Since P.N = 0 identically,
    # -(dP/dlambda).N = P.(dN/dlambda).
    derivative_contribution = clean(
        dot(P,Nprime_transverse)
    )

    torsion_total = clean(
        normal_gradient + derivative_contribution
    )

    # Independent curvature directional derivative.
    Lcurv_e = -v*A*ke**2

    curvature_normal = clean(
        s.diff(Lcurv_e,e).subs(e,0)
    )

    curvature_expected = -2*A*k/v

    assert clean(dot(P,T)) == 0
    assert clean(dot(P,N)) == 0

    assert clean(
        normal_gradient
        -2*C*R2/(v**7*k**3)
    ) == 0

    assert clean(
        derivative_contribution
        +2*C*R2/(v**7*k**3)
    ) == 0

    assert torsion_total == 0
    assert clean(
        curvature_normal-curvature_expected
    ) == 0

    print()
    print("JET", index)
    print("v =", v)
    print("kappa =", k)
    print("TORSION NORMAL GRADIENT =",
          normal_gradient)
    print("TORSION MOMENTUM-DERIVATIVE TERM =",
          derivative_contribution)
    print("TORSION SUM =", torsion_total)
    print("CURVATURE RESIDUAL =",
          clean(curvature_normal-curvature_expected))
    print("PASS", flush=True)

print()
print("R16 COMPLETE")
print("NORMAL-MOMENTUM CANCELLATION: VERIFIED AT TWO JETS")
print("FULL OSTROGRADSKY SYSTEM: NOT YET VERIFIED")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
