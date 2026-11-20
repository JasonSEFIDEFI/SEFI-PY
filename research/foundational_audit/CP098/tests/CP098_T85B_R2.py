"""CP098-T85B: Canonical directional derivative of E.

Use a compact unrestricted geometric expression.
Differentiate along the canonical Hamiltonian vector
field BEFORE restricting to the regular surface.

Compare with the independently recorded T74 transport.
"""

import sympy as sp

print("=== CP098-T85B-R2 ===", flush=True)

v, k = sp.symbols("v k", positive=True)
a,n,r,s,x,y,b,c = sp.symbols(
    "a n r s x y b c", real=True
)
A,B,C = sp.symbols("A B C", real=True, nonzero=True)
DA,DC,H,J,K = sp.symbols("DA DC H J K", real=True)
LT,LN = sp.symbols("LT LN", real=True)
e = sp.symbols("e", real=True)

eta = sp.diag(-1,1,1,1)

def dot(X,Y):
    return (X.T*eta*Y)[0]

def pair(p,X):
    return (p.T*X)[0]

def covdot(p,z):
    return (p.T*eta*z)[0]

def eq(f,g=0):
    return sp.cancel(sp.together(f-g)) == 0

def geometry(u,q):
    speed = sp.sqrt(-dot(u,u))
    T = u/speed
    w = q+T*dot(T,q)
    L = sp.sqrt(dot(w,w))
    N = w/L
    return speed,T,N,L

def compact_E(u,q,P,p1,p0,AA,CC,DDA):
    speed,T,N,L = geometry(u,q)
    tq = dot(T,q)
    pp = dot(P,P)

    return (
        -L/(2*CC)*(
            covdot(p1,P)
            +pair(p1,T)*pair(P,T)
            -pair(p1,N)*pair(P,N)
        )
        -AA*L**2*pair(P,N)/(CC*speed**3)
        -pair(p0,N)/speed
        +2*AA*L*tq/speed**5
        +L*tq*pp/(2*CC*speed)
        +2*L*DDA/speed**3
    )

u = sp.Matrix([v,0,0,0])
q = sp.Matrix([a,v**2*k,0,0])
P = sp.Matrix([0,0,r,s])
p1 = sp.Matrix([0,-2*A*k/v,x,y])

R = r**2+s**2
S = r*x+s*y

p0time = A*k**2+v**4*k**2*R/(4*C)-B
p0 = sp.Matrix([p0time,n,b,c])

E_surface = (
    -2*A*a*k/v**3
    +2*DA*k/v
    -n/v
    -a*k*v*R/(2*C)
    -k*v**2*S/(2*C)
)

assert eq(
    sp.factor(compact_E(u,q,P,p1,p0,A,C,DA)),
    E_surface
)

# Check against T70 without imposing momentum constraints.
Pg = sp.Matrix(sp.symbols("Pg0:4", real=True))
jg = sp.Matrix(sp.symbols("jg0:4", real=True))
zg = sp.Matrix(sp.symbols("zg0:4", real=True))

Pg0,Pg1,Pg2,Pg3 = Pg
jg0,jg1,jg2,jg3 = jg
zg0,zg1,zg2,zg3 = zg

Pg_sq = -Pg0**2+Pg1**2+Pg2**2+Pg3**2

E_T70 = -(
    4*A*C*a*k
    +2*A*Pg1*k**2*v**4
    -4*C*DA*k*v**2
    +2*C*zg1*v**2
    +a*k*v**4*Pg_sq
    +k*v**5*(Pg2*jg2+Pg3*jg3)
)/(2*C*v**3)

off_surface_residual = sp.cancel(
    sp.together(
        compact_E(u,q,Pg,jg,zg,A,C,DA)-E_T70
    )
)

print(
    "OFF_SURFACE_T70_RESIDUAL =",
    off_surface_residual,
    flush=True
)

assert off_surface_residual == 0, (
    "STOP: corrected compact bracket disagrees with T70"
)

print("PASS: unrestricted adapted T70 agreement",
      flush=True)

print("PASS: compact E reproduces T70 surface",
      flush=True)

# Canonical vector field, evaluated at the
# regular surface only AFTER the general
# geometric bracket has been constructed.

udot = q

qdot = sp.Matrix([
    LT,
    LN,
    -v**5*k**2*r/(2*C),
    -v**5*k**2*s/(2*C)
])

Pdot = sp.Matrix([
    0,
    v**3*k*R/(2*C),
    -x-LN*r/(v**2*k),
    -y-LN*s/(v**2*k)
])

p1dot = sp.Matrix([
    2*A*k**2,
    2*A*a*k/v**2-n-a*k*v**2*R/(2*C),
    a*LN*r/(k*v**3)-LT*r/v-b,
    a*LN*s/(k*v**3)-LT*s/v-c
])

p0dot = sp.zeros(4,1)

Adot = v*DA
Cdot = v*DC
DAdot = (
    v*H+v*k**2*J-v**5*k**2*R*K/4
)

# Directional derivative via a first-order
# perturbation of every canonical variable.
# No fixed-frame differentiation of E_surface.

u_e = u+e*udot
q_e = q+e*qdot
P_e = P+e*Pdot
p1_e = p1+e*p1dot
p0_e = p0+e*p0dot

E_e = compact_E(
    u_e,q_e,P_e,p1_e,p0_e,
    A+e*Adot,
    C+e*Cdot,
    DA+e*DAdot
)

print("PASS: canonical perturbation constructed",
      flush=True)

canonical_dot = sp.diff(E_e,e).subs(e,0)

print("PASS: unrestricted E differentiated",
      flush=True)

# Independently reproduce T74's scalar
# transport formula for comparison only.

vdot = a
kdot = LN/v**2-3*a*k/v
adot = LT-v**3*k**2

Rdot = -2*S-2*LN*R/(v**2*k)

ndot = (
    v*k*p0time
    -v**3*k*(b*r+c*s)/(2*C)
)

Sdot = (
    -x**2-y**2-r*b-s*c
    -LN*S/(v**2*k)
    -LT*R/v
    +a*LN*R/(k*v**3)
    -A*k**2*v**2*R/C
)

vv,kk,aa,nn,RR,SS = sp.symbols(
    "vv kk aa nn RR SS", real=True
)
AA,CC,DD = sp.symbols("AA CC DD", real=True)

E_independent = (
    -2*AA*aa*kk/vv**3
    +2*DD*kk/vv
    -nn/vv
    -aa*kk*vv*RR/(2*CC)
    -kk*vv**2*SS/(2*CC)
)

rates = {
    vv:vdot, kk:kdot, aa:adot,
    nn:ndot, RR:Rdot, SS:Sdot,
    AA:Adot, CC:Cdot, DD:DAdot
}

t74_dot = sum(
    sp.diff(E_independent,z)*dz
    for z,dz in rates.items()
)

t74_dot = t74_dot.subs({
    vv:v,kk:k,aa:a,nn:n,RR:R,
    SS:S,AA:A,CC:C,DD:DA
})

print("PASS: T74 comparator constructed",
      flush=True)

residual = sp.cancel(
    sp.together(canonical_dot-t74_dot)
)

print("CANONICAL_T74_RESIDUAL =",residual,
      flush=True)

if residual != 0:
    raise AssertionError(
        "STOP: unrestricted canonical E_dot "
        "disagrees with T74"
    )

print("PASS: exact canonical/T74 agreement",
      flush=True)

# Verify the internal higher-bracket
# coefficients independently survive.

assert eq(sp.diff(canonical_dot,H),2*k)
assert eq(sp.diff(canonical_dot,J),2*k**3)
assert eq(
    sp.diff(canonical_dot,K),
    -v**4*k**3*R/2
)

print("PASS: internal Hessian coefficients",
      flush=True)

# Apply multiplier consistency only after
# establishing the unfixed identity.

FN = (
    -n-6*A*a*k/v**2
    -v**3*k*S/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

LN_fixed = -FN*v**3/(2*A)

fixed_residual = sp.cancel(
    sp.together(
        (canonical_dot-t74_dot).subs(
            LN,LN_fixed
        )
    )
)

assert fixed_residual == 0

print("PASS: fixed-multiplier agreement",
      flush=True)

print("CANONICAL E_DOT: VERIFIED IN REGULAR SECTOR")
print("UNFIXED MULTIPLIER COMPARISON: EXACT")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
print("T85B-R2 COMPLETE")
