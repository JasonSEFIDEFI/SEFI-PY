"""Explicit, dimensionally consistent engineering extension; no original equations replaced."""
import math

def rhs(phi,d,detuning,K0,alpha,beta):
    if K0<=0 or alpha<0 or beta<=0 or d<=-1:raise ValueError('invalid model domain')
    u=-K0*math.sin(phi)/(1+d)
    return detuning+u,alpha*u*u-beta*d

def critical_detuning(K0,alpha,beta):
    if K0<=0 or alpha<0 or beta<=0:raise ValueError('invalid parameters')
    lo=0.;hi=K0
    for _ in range(100):
        x=(lo+hi)/2
        if x*(1+alpha*x*x/beta)<K0:lo=x
        else:hi=x
    return (lo+hi)/2

def equilibrium(detuning,K0,alpha,beta):
    d=alpha*detuning*detuning/beta;s=detuning*(1+d)/K0
    if abs(s)>=1:return None
    return math.asin(s),d

def jacobian(phi,d,detuning,K0,alpha,beta):
    k=K0/(1+d);s=math.sin(phi);c=math.cos(phi)
    return ((-k*c,k*s/(1+d)),(2*alpha*k*k*s*c,-2*alpha*k*k*s*s/(1+d)-beta))

def rk4(phi,d,dt,detuning,K0,alpha,beta):
    def f(p,q):return rhs(p,q,detuning,K0,alpha,beta)
    a,b=f(phi,d);c,e=f(phi+dt*a/2,d+dt*b/2);g,h=f(phi+dt*c/2,d+dt*e/2);j,k=f(phi+dt*g,d+dt*h)
    return phi+dt*(a+2*c+2*g+j)/6,d+dt*(b+2*e+2*h+k)/6

def simulate(detuning,K0,alpha,beta,dt,duration=120.,late_window=30.):
    phi=d=0.;minimum=0.;late=[];samples=[]
    for i in range(round(duration/dt)+1):
        t=i*dt
        if t>=duration-late_window-1e-10:late.append((t,phi,d))
        if i%max(1,round(.1/dt))==0:samples.append((t,phi,d))
        if i==round(duration/dt):break
        phi,d=rk4(phi,d,dt,detuning,K0,alpha,beta);minimum=min(minimum,d)
    drift=(late[-1][1]-late[0][1])/(late[-1][0]-late[0][0]);spread=max(x[1] for x in late)-min(x[1] for x in late)
    return dict(drift=drift,late_phase_range=spread,final_phase=phi,final_degradation=d,min_degradation=minimum,numerically_locked=abs(drift)<.001 and spread<.01),samples
