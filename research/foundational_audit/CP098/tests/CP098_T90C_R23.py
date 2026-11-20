import sympy as s
import importlib.util
import contextlib
import io
from pathlib import Path

print("=== CP098-T90C-R23 ===", flush=True)

source = Path(__file__).with_name("CP098_T90C_R22B.py")
if not source.is_file():
    raise RuntimeError("STOP: R22B source missing")

# Load the original exact-jet implementation unchanged.
# Its first-jet self-tests run during loading.
spec = importlib.util.spec_from_file_location("cp098_r22b", source)
mod = importlib.util.module_from_spec(spec)

with contextlib.redirect_stdout(io.StringIO()) as captured:
    spec.loader.exec_module(mod)

if "R22B COMPLETE" not in captured.getvalue():
    raise RuntimeError("STOP: R22B self-tests did not complete")

print("R22B ENGINE: SELF-TESTS PASSED", flush=True)

# Independent regular fifth-order jet.
mod.u0 = [13, 3, 4, 0]
mod.q0 = [2, 3, 1, 4]
mod.j0 = [1, 4, 2, 5]
mod.a0 = [3, -2, 1, 4]
mod.b0 = [2, 1, -3, 5]

# Retain the same coupling Taylor coefficients.
mod.A0 = s.Integer(3)
mod.A1 = s.Integer(2)
mod.C0 = s.Integer(2)
mod.C1 = s.Integer(1)
mod.C2 = s.Integer(3)

eta = s.diag(-1, 1, 1, 1)

def dot(x, y):
    return (x.T * eta * y)[0]

def zero(x):
    return x.applyfunc(s.simplify) == s.zeros(*x.shape)

u = s.Matrix(mod.u0)
q = s.Matrix(mod.q0)
j = s.Matrix(mod.j0)
b = s.Matrix(mod.b0)

v = s.sqrt(-dot(u, u))
T = u/v
w = q + dot(T, q)*T
w2 = s.simplify(dot(w, w))

assert v.is_positive
assert w2.is_positive

k = s.sqrt(w2)/v**2
N = w/s.sqrt(w2)

projector = (
    s.eye(4)
    + T*(T.T*eta)
    - N*(N.T*eta)
)

print("v =", v)
print("kappa =", k, flush=True)

gu, gqdot, Pddot, p0 = mod.momentum_parts(
    mod.b0, progress=True
)

print("GEOMETRIC p0 =", p0, flush=True)

print("STAGE 4: fifth-derivative response", flush=True)

p0_zero_b = mod.momentum_parts([0,0,0,0])[3]

W = -2*mod.C0*eta*projector/(v*w2)

residual_b = (
    p0-p0_zero_b-W*b
).applyfunc(s.simplify)

assert zero(residual_b)
print("FIFTH-DERIVATIVE RESIDUAL =", residual_b)
print("PASS", flush=True)

print("STAGE 5: C-second-derivative response", flush=True)

p0_zero_c2 = mod.momentum_parts(
    mod.b0, c2=s.Integer(0)
)[3]

R = projector*j
P = -2*mod.C0*eta*R/(v*w2)

residual_c2 = (
    p0-p0_zero_c2-(mod.C2/mod.C0)*P
).applyfunc(s.simplify)

assert zero(residual_c2)
print("C-SECOND-DERIVATIVE RESIDUAL =", residual_c2)
print("PASS", flush=True)

print("STAGE 6: geometric Hamiltonian identity", flush=True)

Lgeom = (
    -mod.A0*w2/v**3
    -mod.C0*dot(R,R)/(v*w2)
)

# Use the previously established R20 projections.
# Since q belongs to span(T,N), these determine p1.q.
p1q = (-2*mod.A0*k/v)*dot(N,q)

p0u = s.simplify((p0.T*u)[0])
Pj = s.simplify((P.T*j)[0])

Hgeom = s.simplify(p0u+p1q+Pj-Lgeom)

print("p0.u =", p0u)
print("p1.q =", s.simplify(p1q))
print("P.j =", Pj)
print("Lgeom =", s.simplify(Lgeom))
print("Hgeom =", Hgeom)

assert Hgeom == 0
print("PASS", flush=True)

print()
print("R23 COMPLETE")
print("SECOND EXACT FIFTH-ORDER JET: VERIFIED")
print("GEOMETRIC p0: EVALUATED AT TWO JETS")
print("FULL p0 INDEPENDENT DERIVATION: OPEN")
print("FULL ACTION MOMENTUM: OPEN")
print("HAMILTONIAN RECONSTRUCTION: OPEN")
print("FULL DIRAC SYSTEM: OPEN")
print("DEFI TORSION IDENTIFICATION: UNRESOLVED")
