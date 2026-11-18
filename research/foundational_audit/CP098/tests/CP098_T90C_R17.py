import sympy as s

print("=== CP098-T90C-R17 ===", flush=True)

eta = s.diag(-1, 1, 1, 1)
A, C, e = s.symbols("A C e", real=True)

def dot(a, b):
    return (a.T * eta * b)[0]

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
    k = s.sqrt(dot(w,w))/v**2
    N = w/(v**2*k)

    assert k.is_positive

    projector = (
        s.eye(4)
        + T*(T.T*eta)
        - N*(N.T*eta)
    )

    R = projector*j
    P = -2*C*R/(v**5*k**2)

    # Directional derivative with respect to q
    # along T, holding u and j fixed.
    qe = q + e*T
    we = qe + dot(T,qe)*T
    ke = s.sqrt(dot(we,we))/v**2
    Ne = we/(v**2*ke)

    proje = (
        s.eye(4)
        + T*(T.T*eta)
        - Ne*(Ne.T*eta)
    )

    Re = proje*j

    Le = (
        -v*A*ke**2
        -C*dot(Re,Re)/(v**5*ke**2)
    )

    tangent_gradient = s.simplify(
        s.diff(Le,e).subs(e,0)
    )

    # T' = w/v = v*k*N.
    Tprime = w/v

    derivative_contribution = s.simplify(
        dot(P,Tprime)
    )

    tangent_total = s.simplify(
        tangent_gradient+derivative_contribution
    )

    assert s.simplify(dot(P,T)) == 0
    assert s.simplify(dot(P,N)) == 0
    assert tangent_gradient == 0
    assert derivative_contribution == 0
    assert tangent_total == 0

    print()
    print("JET", index)
    print("v =", v)
    print("kappa =", k)
    print("TANGENT GRADIENT =", tangent_gradient)
    print("MOMENTUM-DERIVATIVE TERM =",
          derivative_contribution)
    print("TOTAL =", tangent_total)
    print("PASS", flush=True)

print()
print("R17 COMPLETE")
print("TANGENT CONSTRAINT: VERIFIED AT TWO JETS")
print("FULL OSTROGRADSKY SYSTEM: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
