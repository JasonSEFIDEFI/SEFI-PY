"""CP098-T38R: Efficient primary-preservation audit.

Independent directional canonical derivatives.

Canonical pairs:
    (u,p1)
    (q,P)

All directional derivatives are evaluated at the
regular five-constraint surface.

No closure assumption.
"""

import sympy as sp

eta = sp.diag(-1,1,1,1)

v,k = sp.symbols("v k",positive=True)
a,n,b,c,r,s,x,y = sp.symbols(
    "a n b c r s x y",real=True
)

A,B,C = sp.symbols(
    "A B C",real=True,nonzero=True
)

e = sp.symbols("e",real=True)

LT,LN = sp.symbols(
    "lambda_T lambda_N",real=True
)

DA = sp.symbols("DA",real=True)

def dot(X,Y):
    return (X.T*eta*Y)[0]

def pair(p,X):
    return (p.T*X)[0]

def geom(u,q):
    speed = sp.sqrt(-dot(u,u))
    T = u/speed
    tq = dot(T,q)
    w = q+T*tq
    w2 = dot(w,w)
    kap = sp.sqrt(w2)/speed**2
    N = w/sp.sqrt(w2)
    return speed,T,N,kap,w2

def phiT(u,q,P):
    return pair(P,geom(u,q)[1])

def phiN(u,q,P):
    return pair(P,geom(u,q)[2])

def H(u,q,P,p1,p0):
    speed,T,N,kap,w2 = geom(u,q)
    return (
        pair(p0,u)
        +pair(p1,q)
        +speed*(B+A*kap**2)
        -speed*w2*dot(P,P)/(4*C)
    )

def d(expr):
    return sp.factor(
        sp.diff(expr,e).subs(e,0)
    )

u = sp.Matrix([v,0,0,0])
q = sp.Matrix([a,v**2*k,0,0])
P = sp.Matrix([0,0,r,s])
p1 = sp.Matrix([0,-2*A*k/v,x,y])

R = r*r+s*s

p0 = sp.Matrix([
    A*k*k+v**4*k*k*R/(4*C)-B,
    n,b,c
])

speed,T,N,kap,w2 = geom(u,q)

print("=== CP098-T38R ===",flush=True)

assert sp.simplify(phiT(u,q,P)) == 0
assert sp.simplify(phiN(u,q,P)) == 0

print("PASS: primary constraints",flush=True)

# Canonical bracket {f,H}:
#
# f_u . H_p1
# + f_q . H_P
# - f_P . H_q
#
# H_p1=q
# H_P=-speed*w2*eta*P/(2C)
#
# f_P=T or N, respectively.

def bracket_H(f,normal):

    term_u = d(
        f(u+e*q,q,P)
    )

    term_q = d(
        f(
            u,
            q-e*speed*w2*(eta*P)/(2*C),
            P
        )
    )

    term_P = -d(
        H(
            u,
            q+e*normal,
            P,p1,p0
        )
    )

    return (
        sp.factor(term_u),
        sp.factor(term_q),
        sp.factor(term_P),
        sp.factor(term_u+term_q+term_P)
    )

print("Computing phi_T preservation...",flush=True)

TT = bracket_H(phiT,T)

for label,value in zip(
    ["u term","q term","P term","{phi_T,Hc}"],TT
):
    print(label,"=",value,flush=True)

print("Computing phi_N preservation...",flush=True)

NN = bracket_H(phiN,N)

for label,value in zip(
    ["u term","q term","P term","{phi_N,Hc}"],NN
):
    print(label,"=",value,flush=True)

# Primary-primary bracket:
#
# {f,g} = f_q.g_P - f_P.g_q
#
# Both f and g are independent of p1.
#
# Evaluate using independent directional derivatives.

print("Computing primary bracket...",flush=True)

BTN = sp.factor(
    d(phiT(u,q+e*N,P))
    -d(phiN(u,q+e*T,P))
)

BNT = sp.factor(
    d(phiN(u,q+e*T,P))
    -d(phiT(u,q+e*N,P))
)

assert sp.factor(BTN+BNT) == 0

print("{phi_T,phi_N} =",BTN,flush=True)
print("{phi_N,phi_T} =",BNT,flush=True)

eqT = sp.factor(TT[3]+LN*BTN)
eqN = sp.factor(NN[3]+LT*BNT)

print("dot(phi_T) =",eqT,flush=True)
print("dot(phi_N) =",eqN,flush=True)

# Normal multiplier from CP098-T34.
FN_chi = (
    -n
    -6*A*a*k/v**2
    -v**3*k*(r*x+s*y)/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)

LN_normal = -v**3*FN_chi/(2*A)

print(
    "dot(phi_T) after lambda_N =",
    sp.factor(eqT.subs(LN,LN_normal)),
    flush=True
)

print(
    "dot(phi_N) after lambda_N =",
    sp.factor(eqN.subs(LN,LN_normal)),
    flush=True
)

print("PRIMARY AUDIT CALCULATED")
print("REGULAR CLOSURE: PENDING INTERPRETATION")
