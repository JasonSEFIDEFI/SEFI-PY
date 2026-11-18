import sympy as sp

print("=== CP098-T90C-R14 ===")

v, k, C = sp.symbols(
    "v k C", real=True, nonzero=True
)

j0, j1, j2, j3 = sp.symbols(
    "j0 j1 j2 j3", real=True
)

eta = sp.diag(-1, 1, 1, 1)

T = sp.Matrix([1, 0, 0, 0])
N = sp.Matrix([0, 1, 0, 0])
j = sp.Matrix([j0, j1, j2, j3])

def dot(a, b):
    return (a.T * eta * b)[0]

def simplify_matrix(M):
    return M.applyfunc(sp.simplify)

# Mixed-index projector orthogonal to T and N.
projector = (
    sp.eye(4)
    + T * (T.T * eta)
    - N * (N.T * eta)
)

jt = simplify_matrix(projector * j)

L = -C * dot(jt, jt) / (v**5 * k**2)

P = sp.Matrix([
    sp.diff(L, x)
    for x in j
])

Hessian = P.jacobian(j)

expected_P = (
    -2*C/(v**5*k**2) * eta * jt
)

expected_Hessian = (
    -2*C/(v**5*k**2) * eta * projector
)

assert simplify_matrix(P-expected_P) == sp.zeros(4,1)
assert simplify_matrix(Hessian-expected_Hessian) == sp.zeros(4)

assert projector.rank() == 2
assert Hessian.rank() == 2

assert simplify_matrix(projector*T) == sp.zeros(4,1)
assert simplify_matrix(projector*N) == sp.zeros(4,1)

assert sp.simplify((P.T*T)[0]) == 0
assert sp.simplify((P.T*N)[0]) == 0

P2 = sp.simplify((P.T*eta*P)[0])

legendre = sp.simplify((P.T*j)[0]-L)
expected_legendre = -v**5*k**2*P2/(4*C)

assert sp.simplify(
    legendre-expected_legendre
) == 0

print("PROJECTOR =", projector)
print("PROJECTOR RANK =", projector.rank())
print("HESSIAN RANK =", Hessian.rank())
print("MOMENTUM =", P)
print("P.T =", sp.simplify((P.T*T)[0]))
print("P.N =", sp.simplify((P.T*N)[0]))
print("LEGENDRE RESIDUAL =",
      sp.simplify(legendre-expected_legendre))

print()
print("PASS: transverse momentum")
print("PASS: Hessian rank two")
print("PASS: two Hessian null directions")
print("PASS: exact torsion Legendre identity")
print()
print("R14 COMPLETE")
print("LIMIT: conditional torsion identification")
print("FULL ACTION EQUIVALENCE: NOT ESTABLISHED")
