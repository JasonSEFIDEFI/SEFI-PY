"""CP098-T58: Canonical evolution of t=v*kappa.

Compute {t,Hc}, {t,phi_T}, {t,phi_N}
by independent directional differentiation.

External momenta unrestricted until surface evaluation.
No gauge generator or coefficient evolution assumed.
"""

import sympy as sp

eta = sp.diag(-1, 1, 1, 1)

v, k = sp.symbols("v k", positive=True)
a = sp.symbols("a", real=True)
A, B, C = sp.symbols(
    "A B C", real=True, nonzero=True
)
LT, LN = sp.symbols(
    "lambda_T lambda_N", real=True
)
e = sp.symbols("e", real=True)

Pz = sp.symbols("P0:4", real=True)
p1z = sp.symbols("p10:14", real=True)
p0z = sp.symbols("p00:04", real=True)

P = sp.Matrix(Pz)
p1 = sp.Matrix(p1z)
p0 = sp.Matrix(p0z)

u = sp.Matrix([v, 0, 0, 0])
q = sp.Matrix([a, v**2*k, 0, 0])

def dot(X, Y):
    return (X.T*eta*Y)[0]

def pair(p, X):
    return (p.T*X)[0]

def geom(uu, qq):
    speed = sp.sqrt(-dot(uu, uu))
    T = uu/speed
    w = qq+T*dot(T,qq)
    w2 = dot(w,w)
    N = w/sp.sqrt(w2)
    kap = sp.sqrt(w2)/speed**2
    return speed,T,N,kap,w2

def tfun(uu,qq):
    speed,T,N,kap,w2 = geom(uu,qq)
    return speed*kap

def H(uu,qq):
    speed,T,N,kap,w2 = geom(uu,qq)
    return (
        pair(p0,uu)
        +pair(p1,qq)
        +speed*(B+A*kap**2)
        -speed*w2*dot(P,P)/(4*C)
    )

def phiT(uu,qq):
    speed,T,N,kap,w2 = geom(uu,qq)
    return pair(P,T)

def phiN(uu,qq):
    speed,T,N,kap,w2 = geom(uu,qq)
    return pair(P,N)

def D(expr):
    return sp.factor(sp.diff(expr,e).subs(e,0))

print("=== CP098-T58 ===", flush=True)

speed,T,N,kap,w2 = geom(u,q)

assert sp.simplify(tfun(u,q)-v*k)==0

print("PASS: adapted geometric definition",flush=True)

# For t=t(u,q), the canonical bracket is
#
# {t,F} = t_u.F_p1 + t_q.F_P.
#
# No derivatives of t with respect to momenta,
# internal variables, or y are present.

H_P = -speed*w2*eta*P/(2*C)

term_u = D(tfun(u+e*q,q))
term_q = D(tfun(u,q+e*H_P))

bracket_H = sp.factor(term_u+term_q)

print("t_u.H_p1 =",term_u,flush=True)
print("t_q.H_P =",term_q,flush=True)
print("{t,Hc} =",bracket_H,flush=True)

# Independently construct the momentum gradients
# of each primary constraint:
#
# (phi_T)_P=T, (phi_N)_P=N.

bracket_phiT = D(tfun(u,q+e*T))
bracket_phiN = D(tfun(u,q+e*N))

print("{t,phi_T} =",bracket_phiT,flush=True)
print("{t,phi_N} =",bracket_phiN,flush=True)

# Geometrically, shifting q in the T direction
# leaves w and t unchanged.
assert sp.factor(bracket_phiT)==0

# Shifting q along N increases sqrt(w^2)
# by one unit per parameter.
assert sp.factor(bracket_phiN-1/v)==0

print("PASS: independent primary brackets",flush=True)

total = sp.factor(
    bracket_H+LT*bracket_phiT+LN*bracket_phiN
)

print("{t,H_T} =",total,flush=True)

# Evaluate only now on the full regular family.
r,s,x,y,n,b,c = sp.symbols(
    "r s x y n b c",real=True
)

surface = {
    Pz[0]:0,
    Pz[1]:0,
    Pz[2]:r,
    Pz[3]:s,
    p1z[0]:0,
    p1z[1]:-2*A*k/v,
    p1z[2]:x,
    p1z[3]:y,
    p0z[0]:
        A*k**2+v**4*k**2*(r*r+s*s)/(4*C)-B,
    p0z[1]:n,
    p0z[2]:b,
    p0z[3]:c
}

reduced_H = sp.factor(bracket_H.subs(surface))
reduced_total = sp.factor(total.subs(surface))

print("Surface {t,Hc} =",reduced_H,flush=True)
print("Surface {t,H_T} =",reduced_total,flush=True)

assert sp.factor(
    reduced_total-reduced_H-LN/v
)==0

assert LT not in reduced_total.free_symbols

print("PASS: total evolution decomposition",flush=True)
print("PASS: tangential multiplier absent",flush=True)

# The normal multiplier is fixed by T50.
DA=sp.symbols("DA",real=True)

FN = (
    -n-6*A*a*k/v**2
    -v**3*k*(r*x+s*y)/(2*C)
    -v**2*k*a*(r*r+s*s)/(2*C)
    +2*k*DA
)

LN_fixed = -v**3*FN/(2*A)

fixed = sp.factor(
    reduced_total.subs(LN,LN_fixed)
)

print("Fixed-multiplier t evolution =",fixed,flush=True)

assert LN not in fixed.free_symbols

print("PASS: fixed-multiplier expression obtained",flush=True)

# Cross-check with a direct gradient construction
# at the adapted point for the u contribution:
# t = sqrt(w^2)/v, and u varies along q.
#
# This second route does not use the directional
# derivative of tfun.

u_grad = sp.Matrix([
    D(tfun(u+e*sp.eye(4)[:,i],q))
    for i in range(4)
])

q_grad = sp.Matrix([
    D(tfun(u,q+e*sp.eye(4)[:,i]))
    for i in range(4)
])

gradient_bracket = sp.factor(
    (u_grad.T*q)[0]+(q_grad.T*H_P)[0]
)

assert sp.factor(gradient_bracket-bracket_H)==0

print("PASS: independent gradient contraction",flush=True)

print("GEOMETRIC COEFFICIENT EVOLUTION: VERIFIED")
print("E/d EVOLUTION: NOT YET CALCULATED")
print("TIME-DEPENDENT GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
