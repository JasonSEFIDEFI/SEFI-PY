"""CP098-T64: Moving-frame transport audit.

Calculate the canonical directional evolution of
v, kappa, T, and N on the regular constraint surface.

Independent checks:
  1. Direct differentiation of geometric definitions.
  2. Differentiation of w and its normalization.
  3. Orthonormality preservation.
  4. Consistency with T58 and T63.

No gauge generator or action invariance assumed.
"""

import sympy as sp

print("=== CP098-T64 ===", flush=True)

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a, r, s = sp.symbols("a r s", real=True)
C = sp.symbols("C", real=True, nonzero=True)
LT, LN = sp.symbols(
    "lambda_T lambda_N", real=True
)
e = sp.symbols("e", real=True)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])
P = sp.Matrix([0, 0, r, s])

def dot(X, Y):
    return (X.T*eta*Y)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq + T*dot(T, qq)
    wn = sp.sqrt(dot(w, w))
    N = w/wn
    kap = wn/speed**2
    return speed, T, w, wn, N, kap

def simplify_matrix(M):
    return M.applyfunc(sp.simplify)

def D(expr):
    return sp.simplify(
        sp.diff(expr, e).subs(e, 0)
    )

def DM(M):
    return M.applyfunc(D)

speed, T, w, wn, N, kap = geom(u, q)

assert sp.simplify(speed-v) == 0
assert sp.simplify(kap-k) == 0
assert sp.simplify(wn-v**2*k) == 0

print("PASS: adapted geometric definitions",
      flush=True)

# Hamilton equations evaluated on the regular
# five-constraint family.

udot = q

qdot = (
    -v**5*k**2*eta*P/(2*C)
    + LT*T
    + LN*N
)

expected_qdot = sp.Matrix([
    LT,
    LN,
    -v**5*k**2*r/(2*C),
    -v**5*k**2*s/(2*C)
])

assert simplify_matrix(
    qdot-expected_qdot
) == sp.zeros(4, 1)

print("PASS: canonical evolution vectors",
      flush=True)

# Route 1: differentiate the full geometric
# definitions while evolving u and q together.

ue = u+e*udot
qe = q+e*qdot

ve, Te, we, wne, Ne, ke = geom(ue, qe)

vdot = D(ve)
Tdot = DM(Te)
wdot = DM(we)
wndot = D(wne)
Ndot = DM(Ne)
kdot = D(ke)

expected_vdot = a

expected_Tdot = sp.Matrix([
    0, v*k, 0, 0
])

expected_wdot = sp.Matrix([
    v**3*k**2,
    LN-a*v*k,
    -v**5*k**2*r/(2*C),
    -v**5*k**2*s/(2*C)
])

expected_wndot = LN-a*v*k

expected_Ndot = sp.Matrix([
    v*k,
    0,
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

expected_kdot = (
    LN/v**2-3*a*k/v
)

assert sp.simplify(
    vdot-expected_vdot
) == 0

assert simplify_matrix(
    Tdot-expected_Tdot
) == sp.zeros(4, 1)

assert simplify_matrix(
    wdot-expected_wdot
) == sp.zeros(4, 1)

assert sp.simplify(
    wndot-expected_wndot
) == 0

assert simplify_matrix(
    Ndot-expected_Ndot
) == sp.zeros(4, 1)

assert sp.simplify(
    kdot-expected_kdot
) == 0

print("PASS: direct geometric differentiation",
      flush=True)

print("v_dot =", vdot, flush=True)
print("T_dot =", Tdot, flush=True)
print("w_dot =", wdot, flush=True)
print("|w|_dot =", wndot, flush=True)
print("N_dot =", Ndot, flush=True)
print("kappa_dot =", kdot, flush=True)

# Route 2: reconstruct N_dot from
# N=w/|w|, using independently differentiated
# w and |w|.

Ndot_reconstructed = simplify_matrix(
    wdot/wn-w*wndot/wn**2
)

assert simplify_matrix(
    Ndot_reconstructed-Ndot
) == sp.zeros(4, 1)

print("PASS: independent normal-vector reconstruction",
      flush=True)

# Check the defining frame relations and their
# time derivatives.

assert sp.simplify(dot(T,T)+1) == 0
assert sp.simplify(dot(N,N)-1) == 0
assert sp.simplify(dot(T,N)) == 0

assert sp.simplify(
    2*dot(T,Tdot)
) == 0

assert sp.simplify(
    2*dot(N,Ndot)
) == 0

assert sp.simplify(
    dot(Tdot,N)+dot(T,Ndot)
) == 0

print("PASS: orthonormality preserved",
      flush=True)

# T58: t=v*kappa.
# Its total derivative must agree with the
# independently established surface expression.

tdot = sp.factor(vdot*k+v*kdot)

expected_tdot = (
    -2*a*k+LN/v
)

assert sp.simplify(
    tdot-expected_tdot
) == 0

print("PASS: T58 geometric evolution recovered",
      flush=True)

# T63: tangent transport.
assert simplify_matrix(
    Tdot-v*k*N
) == sp.zeros(4, 1)

print("PASS: T63 tangent transport recovered",
      flush=True)

# Normal-vector transport has transverse terms
# depending on P. These must not be discarded.

transverse = sp.Matrix([
    Ndot[2], Ndot[3]
])

expected_transverse = sp.Matrix([
    -v**3*k*r/(2*C),
    -v**3*k*s/(2*C)
])

assert simplify_matrix(
    transverse-expected_transverse
) == sp.zeros(2, 1)

print("PASS: transverse normal transport retained",
      flush=True)

print("MOVING-FRAME TRANSPORT: VERIFIED")
print("v_dot = a")
print("T_dot = v*k*N")
print("kappa_dot = LN/v^2 - 3*a*k/v")
print(
    "N_dot = v*k*T "
    "- v^3*k*(r*e2+s*e3)/(2*C)"
)

print("FULL CANONICAL E EVOLUTION: OPEN")
print("FRAME-TRANSPORT DEPENDENCIES: VERIFIED")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
