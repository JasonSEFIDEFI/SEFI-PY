"""CP098-T75: Numerical off-surface canonical E-dot audit.

Independent finite-difference evaluation of:
    E = {chi_N, psi_T}
    E_dot = directional derivative of E along X_H

Compare against T74 surface-chain-rule transport.

One internal DOF, nonconstant A and C.
Numerical evidence only, not a symbolic proof.
"""

import numpy as np
import sympy as sp

print("=== CP098-T75 ===", flush=True)

# State layout:
# u[0:4], q[4:8], P[8:12], p1[12:16],
# p0[16:20], I[20], pI[21]

eta = np.diag([-1., 1., 1., 1.])

def metric(x, y):
    return float(x @ eta @ y)

def geometry(z):
    u = z[0:4]
    q = z[4:8]
    vv = np.sqrt(-metric(u, u))
    T = u/vv
    w = q+T*metric(T,q)
    kk = np.sqrt(metric(w,w))/vv**2
    N = w/(vv**2*kk)
    return vv,kk,T,N

def internal(z):
    I, pI = z[20:22]
    A = 1.3+0.2*I+0.07*I**2
    C = 1.7+0.1*I+0.04*I**2
    B = 1.1+pI**2/(2*1.4)+0.3*I**2
    return A,B,C

def chi(z):
    v,k,T,N = geometry(z)
    A,B,C = internal(z)
    return z[12:16]@N+2*A*k/v

def psi(z):
    v,k,T,N = geometry(z)
    A,B,C = internal(z)
    return (
        z[16:20]@T+B-A*k*k
        -v**4*k*k*metric(z[8:12],z[8:12])/(4*C)
    )

def hc(z):
    v,k,T,N = geometry(z)
    A,B,C = internal(z)
    return (
        z[16:20]@z[0:4]
        +z[12:16]@z[4:8]
        +v*(B+A*k*k)
        -v**5*k*k*metric(z[8:12],z[8:12])/(4*C)
    )

def phiT(z):
    return z[8:12]@geometry(z)[2]

def phiN(z):
    return z[8:12]@geometry(z)[3]

def gradient(f,z,h):
    out = np.zeros(len(z))
    for i in range(len(z)):
        zp=z.copy()
        zm=z.copy()
        zp[i]+=h
        zm[i]-=h
        out[i]=(f(zp)-f(zm))/(2*h)
    return out

def poisson_grad(f,g):
    # Canonical ordering:
    # (u,p1), (q,P), (I,pI).
    return (
        np.dot(f[0:4],g[12:16])
        -np.dot(f[12:16],g[0:4])
        +np.dot(f[4:8],g[8:12])
        -np.dot(f[8:12],g[4:8])
        +f[20]*g[21]-f[21]*g[20]
    )

def canonical_E(z,h):
    return poisson_grad(
        gradient(chi,z,h),
        gradient(psi,z,h)
    )

def flow(z,LT,LN,h):
    def total(zz):
        return (
            hc(zz)+LT*phiT(zz)+LN*phiN(zz)
        )
    g=gradient(total,z,h)
    dz=np.zeros_like(z)
    dz[0:4]=g[12:16]
    dz[12:16]=-g[0:4]
    dz[4:8]=g[8:12]
    dz[8:12]=-g[4:8]
    dz[20]=g[21]
    dz[21]=-g[20]
    return dz

# Explicit regular adapted point.
v=1.2
k=0.7
a=0.15
r,s=0.23,-0.17
x,y=0.31,-0.19
n=0.27
b,c=0.11,-0.14
I,pI=0.4,0.35
LT=0.12

z=np.zeros(22)
z[0:4]=[v,0,0,0]
z[4:8]=[a,v*v*k,0,0]
z[8:12]=[0,0,r,s]
z[20:22]=[I,pI]

A,B,C=internal(z)
R=r*r+s*s
S=r*x+s*y

z[12:16]=[0,-2*A*k/v,x,y]
z[16:20]=[
    A*k*k+v**4*k*k*R/(4*C)-B,
    n,b,c
]

assert abs(chi(z))<1e-11
assert abs(psi(z))<1e-11
assert abs(phiT(z))<1e-11
assert abs(phiN(z))<1e-11

print("PASS: regular adapted point",flush=True)

# Exact internal derivatives for the chosen
# illustrative internal Hamiltonian.

ii,pp=sp.symbols("ii pp",real=True)
As=sp.Rational(13,10)+ii/5+7*ii**2/100
Cs=sp.Rational(17,10)+ii/10+ii**2/25
Bs=sp.Rational(11,10)+pp**2/sp.Rational(28,10)+3*ii**2/10

def ib(f,g):
    return sp.diff(f,ii)*sp.diff(g,pp)-sp.diff(f,pp)*sp.diff(g,ii)

DAs=ib(As,Bs)
DCs=ib(Cs,Bs)
Hs=ib(DAs,Bs)
Js=ib(DAs,As)
Ks=ib(DAs,1/Cs)

subs={ii:I,pp:pI}

DA=float(DAs.subs(subs))
DC=float(DCs.subs(subs))
H=float(Hs.subs(subs))
J=float(Js.subs(subs))
K=float(Ks.subs(subs))

FN=(
    -n-6*A*a*k/v**2
    -v**3*k*S/(2*C)
    -v**2*k*a*R/(2*C)
    +2*k*DA
)
LN=-FN*v**3/(2*A)

print("PASS: internal brackets and multiplier",flush=True)

# Independent T74 chain-rule prediction.
vv,kk,aa,nn,RR,SS=sp.symbols(
    "vv kk aa nn RR SS",real=True
)
AA,BB,CC,DD,DDC,HH,JJ,KK=sp.symbols(
    "AA BB CC DD DDC HH JJ KK",real=True
)
ltt,lnn=sp.symbols("ltt lnn",real=True)

Es=(
    -2*AA*aa*kk/vv**3
    +2*DD*kk/vv
    -nn/vv
    -aa*kk*vv*RR/(2*CC)
    -kk*vv**2*SS/(2*CC)
)

rates={
    vv:aa,
    kk:lnn/vv**2-3*aa*kk/vv,
    aa:ltt-vv**3*kk**2,
    nn:vv*kk*(
        AA*kk**2+vv**4*kk**2*RR/(4*CC)-BB
    )-vv**3*kk*(b*r+c*s)/(2*CC),
    RR:-2*SS-2*lnn*RR/(vv**2*kk),
    SS:(
        -x*x-y*y-r*b-s*c
        -lnn*SS/(vv**2*kk)
        -ltt*RR/vv
        +aa*lnn*RR/(kk*vv**3)
        -AA*kk**2*vv**2*RR/CC
    ),
    AA:vv*DD,
    CC:vv*DDC,
    DD:(
        vv*HH+vv*kk**2*JJ
        -vv**5*kk**2*RR*KK/4
    )
}

Edot_symbolic=sum(
    sp.diff(Es,t)*rate
    for t,rate in rates.items()
)

values={
    vv:v,kk:k,aa:a,nn:n,RR:R,SS:S,
    AA:A,BB:B,CC:C,DD:DA,DDC:DC,
    HH:H,JJ:J,KK:K,ltt:LT,lnn:LN
}

expected=float(Edot_symbolic.subs(values))
expected_E=float(Es.subs(values))

print("PASS: T74 reference independently assembled",
      flush=True)

# Numerical off-surface canonical derivative.
# Nested central differences: compare across
# multiple step sizes, not a single tolerance.

results=[]

for h in [0.003,0.002,0.001]:
    dt=4*h
    dz=flow(z,LT,LN,h/2)

    E0=canonical_E(z,h)

    Ep=canonical_E(z+dt*dz,h)
    Em=canonical_E(z-dt*dz,h)

    numerical=(Ep-Em)/(2*dt)

    err=abs(numerical-expected)
    results.append(err)

    print(
        f"h={h:.4g} "
        f"E={E0:.10g} "
        f"E_dot={numerical:.10g} "
        f"error={err:.3e}",
        flush=True
    )

    if abs(E0-expected_E)>0.005:
        raise AssertionError(
            "STOP: canonical E differs from T70"
        )

# The calculation must be reasonably accurate
# and improve as the step size decreases.
if not (
    results[-1]<0.01
    and results[-1]<results[0]
):
    raise AssertionError(
        "STOP: independent E_dot comparison failed"
    )

print("PASS: off-surface numerical convergence",
      flush=True)
print("T74 E_dot =",expected,flush=True)

print("CANONICAL E_DOT NUMERICAL CROSS-CHECK: PASS")
print("GENERAL SYMBOLIC CANONICAL PROOF: OPEN")
print("JACOBI CROSS-CHECK: OPEN")
print("GAUGE GENERATOR: OPEN")
print("ACTION INVARIANCE: OPEN")
print("SINGULAR SECTORS: OPEN")
print("PHYSICAL VALIDATION: NOT ESTABLISHED")
