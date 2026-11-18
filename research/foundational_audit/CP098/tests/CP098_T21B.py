"""CP098-T21B: Five-parameter adapted-family closure.

v>0, k>0; a,t,n real.
Calculates all three canonical-pair contributions.
No general constraint-surface closure claim.
"""

import sympy as sp

v, k = sp.symbols("v k", positive=True)
a, t, n = sp.symbols("a t n", real=True)
z, Pz, I, ps = sp.symbols("z Pz I ps", real=True)

# Established independently in CP098-T21A.
part_u = (a*k*v**4 + 4*a*k + 2*n*v**2)/(2*v**3)

# Adapted geometry with a transverse q perturbation z.
# Keep z and Pz independent until differentiation.
w2 = v**4*k**2 + z**2
kappa = sp.sqrt(w2)/v**2

# Transverse canonical momentum P=(0,0,Pz,0).
# Other P components vanish only after the derivative.
psi_P = -w2*Pz**2/4

# p1=(0,-2k/v,0,0).
# The normal component changes as q changes.
chi_q = (-2*k/v)*(v**2*k)/sp.sqrt(w2) + 2*I*kappa/v

dpsi_dPz = sp.diff(psi_P, Pz)
dchi_dz = sp.diff(chi_q, z)

# {psi,chi}_{q,P} = -dpsi/dP * dchi/dq
part_q = sp.simplify(
    (-dpsi_dPz*dchi_dz).subs({
        z: 0, Pz: 1, I: 1
    })
)

# Internal canonical pair (I,ps):
# psi contains ps^2/2 - I*kappa^2.
# chi contains 2*I*kappa/v.
psi_internal = ps**2/2 - I*kappa**2
chi_internal = 2*I*kappa/v

part_I = (
    sp.diff(psi_internal, I)
    * sp.diff(chi_internal, ps)
    - sp.diff(psi_internal, ps)
    * sp.diff(chi_internal, I)
)

part_I = sp.simplify(part_I.subs({
    z: 0, I: 1, ps: t
}))

direct = sp.factor(part_u + part_q + part_I)

# Independent normal-preservation expression,
# evaluated for this adapted constrained family.
# T.q=-a, p0.N=n, p1.P=0, P^2=1, D_A=t.
F_N = (
    -n
    + 2*k*t
    - 6*k*a/v**2
    - v**2*k*a/2
)

required = sp.factor(
    -4*a*k/v**3 - F_N/v
)

difference = sp.factor(direct - required)

print("=== CP098-T21B ===", flush=True)
print("u,p1 contribution =", part_u, flush=True)
print("q,P contribution =", part_q, flush=True)
print("Internal contribution =", part_I, flush=True)
print("Direct bracket =", direct, flush=True)
print("T16 requirement =", required, flush=True)
print("Difference =", difference, flush=True)

assert part_q == 0
assert difference == 0

print("PASS: five-parameter adapted-family identity")
print("GENERAL DIRAC CLOSURE: NOT ESTABLISHED")
