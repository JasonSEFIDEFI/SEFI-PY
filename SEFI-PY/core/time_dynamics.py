"""Classical identity reduction from the Time manuscript, Jason D. Dutton.

Coordinates and the evolution parameter are dimensionless. Geometry is prescribed,
not varied. This module supplies model calculations, not empirical validation,
physical-time calibration, quantum dynamics, or an attracting synchronization law.
Only the Python standard library is required.
"""
from dataclasses import dataclass
from fractions import Fraction
import cmath
import math


def _finite(*values):
    if not all(math.isfinite(x) for x in values):
        raise ValueError("Inputs must be finite")


@dataclass(frozen=True)
class IdentityMetric:
    alpha: float = 1.0
    beta: float = 1.0
    gamma: float = 1.0
    delta: float = 1.0

    def __post_init__(self):
        _finite(*self.diagonal)
        if min(self.diagonal) <= 0:
            raise ValueError("Identity metric coefficients must be positive")

    @property
    def diagonal(self):
        return self.alpha, self.beta, self.gamma, self.delta


@dataclass(frozen=True)
class Coupling:
    """Supply a C2 coupling and its two analytic derivatives."""
    value: object
    first: object
    second: object

    def __post_init__(self):
        if not all(callable(f) for f in (self.value,self.first,self.second)):
            raise ValueError("Coupling and both derivatives must be callable")

    @classmethod
    def linear(cls, strength):
        _finite(strength)
        return cls(lambda x: strength*x, lambda x: strength, lambda x: 0.0)

    @classmethod
    def quadratic(cls, strength):
        _finite(strength)
        return cls(lambda x: strength*x*x, lambda x: 2*strength*x,
                   lambda x: 2*strength)


def _state(state):
    if len(state) != 8:
        raise ValueError("State order: IO, IA, IS, IW, vO, vA, vS, vW")
    _finite(*state)


@dataclass(frozen=True)
class IdentityModel:
    metric: IdentityMetric
    mu1: Coupling
    mu2: Coupling

    def rhs(self, state, kappa, sigma):
        """Base identity equations for the instantaneous prescribed geometry."""
        _state(state); _finite(kappa, sigma)
        o, a, s, w, vo, va, vs, vw = state
        m = self.metric
        return (vo, va, vs, vw, -o,
                -a-(self.mu1.first(a)*kappa**2*s+self.mu2.first(a)*sigma**2*w)/m.beta,
                -s-self.mu1.value(a)*kappa**2/m.gamma,
                -w-self.mu2.value(a)*sigma**2/m.delta)

    def energy(self, state, kappa, sigma):
        """Conserved for constant geometry; otherwise external work is possible."""
        _state(state); _finite(kappa, sigma)
        _, a, s, w = state[:4]
        return sum(m*(q*q+v*v)/2 for m,q,v in
                   zip(self.metric.diagonal,state[:4],state[4:])) + \
            self.mu1.value(a)*kappa**2*s+self.mu2.value(a)*sigma**2*w

    def prescribed_power(self, state, kappa, sigma, kappa_rate, sigma_rate):
        """Explicit partial-tau derivative of the potential along a solution."""
        _state(state); _finite(kappa,sigma,kappa_rate,sigma_rate)
        _, a, s, w = state[:4]
        return 2*self.mu1.value(a)*kappa*kappa_rate*s + \
            2*self.mu2.value(a)*sigma*sigma_rate*w

    def effective_potential(self, a, kappa, sigma):
        """Algebraic completion of squares, not dynamical elimination."""
        _finite(a,kappa,sigma); m=self.metric
        return m.beta*a*a/2-kappa**4*self.mu1.value(a)**2/(2*m.gamma) \
            -sigma**4*self.mu2.value(a)**2/(2*m.delta)

    def effective_derivatives(self, a, kappa, sigma):
        _finite(a,kappa,sigma); m=self.metric
        first=m.beta*a; second=m.beta
        for mu, factor in [(self.mu1,kappa**4/m.gamma),(self.mu2,sigma**4/m.delta)]:
            first-=factor*mu.value(a)*mu.first(a)
            second-=factor*(mu.first(a)**2+mu.value(a)*mu.second(a))
        return first, second

    def equilibrium(self, a, kappa, sigma, tolerance=1e-10):
        """Return a stationary state only when Veff'(a)=0 within tolerance."""
        if tolerance <= 0: raise ValueError("Tolerance must be positive")
        if abs(self.effective_derivatives(a,kappa,sigma)[0]) > tolerance:
            raise ValueError("Candidate is not an equilibrium")
        return (0.,a,-self.mu1.value(a)*kappa**2/self.metric.gamma,
                -self.mu2.value(a)*sigma**2/self.metric.delta,0.,0.,0.,0.)

    def equilibrium_spectrum(self, a, kappa, sigma):
        """Normalized three-sector Hessian; exclude the decoupled Origin mode."""
        state=self.equilibrium(a,kappa,sigma); m=self.metric
        d=1+(self.mu1.second(a)*kappa**2*state[2]+
             self.mu2.second(a)*sigma**2*state[3])/m.beta
        aa=self.mu1.first(a)*kappa**2/math.sqrt(m.beta*m.gamma)
        bb=self.mu2.first(a)*sigma**2/math.sqrt(m.beta*m.delta)
        rad=math.hypot(d-1,2*math.hypot(aa,bb))
        return tuple(sorted((1.,(d+1-rad)/2,(d+1+rad)/2)))


def origin_energy(q, velocity, alpha=1., omega=1.):
    _finite(q,velocity,alpha,omega)
    if alpha <= 0 or omega <= 0: raise ValueError("Positive alpha and omega required")
    return alpha*(velocity**2+omega**2*q**2)/2


def clock_phase(q, velocity, omega=1.):
    """Directed oscillator phase in [0,2pi); undefined at zero amplitude."""
    origin_energy(q,velocity,omega=omega)
    if q == 0 and velocity == 0: raise ValueError("Zero-amplitude clock has no phase")
    return math.atan2(-velocity/omega,q) % math.tau


def phase_history(times, coordinates, velocities, omega=1., initial_phase=None):
    """Lift sampled harmonic phase; reject gaps that can alias cycle counts.

    Requires the stated constant-frequency oscillator model and adequate sampling.
    Cycle history is additional memory, not an instantaneous-state observable.
    """
    if not (len(times)==len(coordinates)==len(velocities)) or not times:
        raise ValueError("Nonempty arrays of equal length required")
    _finite(*times)
    phases=[clock_phase(q,v,omega) for q,v in zip(coordinates,velocities)]
    initial=phases[0] if initial_phase is None else initial_phase
    _finite(initial)
    if abs(math.remainder(initial-phases[0],math.tau)) > 1e-9:
        raise ValueError("Initial representative must match the measured phase")
    lifted=[initial]
    for i in range(1,len(times)):
        gap=times[i]-times[i-1]
        if gap <= 0 or omega*gap >= math.pi:
            raise ValueError("Increasing times and less than half-cycle gaps required")
        increment=math.remainder(phases[i]-phases[i-1],math.tau)
        if abs(increment-omega*gap) > 1e-7:
            raise ValueError("Record does not match the specified harmonic clock")
        lifted.append(lifted[-1]+increment)
    return tuple(lifted)


def relational_derivative(state_rate, phase_rate):
    _finite(phase_rate,*state_rate)
    if phase_rate <= 0: raise ValueError("Positive phase speed required for this chart")
    return tuple(x/phase_rate for x in state_rate)


def linear_parameters(g1,g2,kappa,sigma,metric=IdentityMetric()):
    _finite(g1,g2,kappa,sigma)
    a=g1*kappa**2/math.sqrt(metric.beta*metric.gamma)
    b=g2*sigma**2/math.sqrt(metric.beta*metric.delta)
    return a,b,math.hypot(a,b)


def linear_spectrum(a,b):
    _finite(a,b); r=math.hypot(a,b)
    return 1-r,1.,1+r


def linear_modes(a,b):
    """Columns: soft, orthogonal, hard. At r=0 use the coordinate basis."""
    _finite(a,b); r=math.hypot(a,b)
    if r == 0: return ((1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    h=1/math.sqrt(2)
    return ((h,-h*a/r,-h*b/r),(0.,-b/r,a/r),(h,h*a/r,h*b/r))


def linear_solution(a,b,position,velocity,tau):
    if len(position)!=3 or len(velocity)!=3: raise ValueError("Three coordinates required")
    _finite(tau,*position,*velocity)
    out=[0.]*3; rate=[0.]*3
    for lam,e in zip(linear_spectrum(a,b),linear_modes(a,b)):
        q=sum(x*y for x,y in zip(position,e)); v=sum(x*y for x,y in zip(velocity,e))
        if lam > 0:
            f=math.sqrt(lam); c=math.cos(f*tau); s=math.sin(f*tau)
            qt=q*c+v*s/f; vt=-q*f*s+v*c
        elif lam == 0: qt=q+v*tau; vt=v
        else:
            f=math.sqrt(-lam); c=math.cosh(f*tau); s=math.sinh(f*tau)
            qt=q*c+v*s/f; vt=q*f*s+v*c
        for i in range(3): out[i]+=qt*e[i]; rate[i]+=vt*e[i]
    return tuple(out),tuple(rate)


def recurrence_cycles(ratios):
    """Exact common cycles for supplied positive Fraction frequency ratios.

    Floating-point near returns never establish exact commensurability. Include
    only excited modes; integer ratios are required for a one-cycle observable.
    """
    if not ratios: raise ValueError("At least one excited frequency required")
    cycles=1
    for ratio in ratios:
        if not isinstance(ratio,Fraction) or ratio <= 0:
            raise ValueError("Supply exact positive Fraction ratios")
        cycles=math.lcm(cycles,ratio.denominator)
    return cycles


def exact_linear_ratios(r):
    """Rational r certificate relative to a unit-frequency Origin clock."""
    if not isinstance(r,Fraction) or not 0 <= r < 1:
        raise ValueError("Supply exact Fraction r with 0 <= r < 1")
    roots=[]
    for squared in (1-r,Fraction(1),1+r):
        n=math.isqrt(squared.numerator); d=math.isqrt(squared.denominator)
        if n*n!=squared.numerator or d*d!=squared.denominator:
            raise ValueError("Frequency ratio is irrational; no rational certificate")
        roots.append(Fraction(n,d))
    return tuple(roots)


def locking_residual(phase,reference_phase,ratio,offset=0.):
    """Kinematic lift-lock residual; no attracting interaction is implied."""
    _finite(phase,reference_phase,ratio,offset)
    if ratio <= 0: raise ValueError("Positive frequency ratio required")
    return phase-ratio*reference_phase-offset


def periodic_descent_residual(observable,phase):
    """A pointwise diagnostic, not proof of descent on the entire trajectory."""
    _finite(phase)
    return observable(phase+math.tau)-observable(phase)


def rk4_step(rhs,t,state,step):
    _finite(t,step,*state)
    if step <= 0: raise ValueError("Positive integration step required")
    y=tuple(state)
    if not y: raise ValueError("Nonempty state required")
    def evaluate(time,values):
        result=tuple(rhs(time,values))
        if len(result)!=len(y): raise ValueError("RHS dimension mismatch")
        _finite(*result)
        return result
    k1=evaluate(t,y)
    k2=evaluate(t+step/2,tuple(v+step*k/2 for v,k in zip(y,k1)))
    k3=evaluate(t+step/2,tuple(v+step*k/2 for v,k in zip(y,k2)))
    k4=evaluate(t+step,tuple(v+step*k for v,k in zip(y,k3)))
    result=tuple(v+step*(a+2*b+2*c+d)/6 for v,a,b,c,d in zip(y,k1,k2,k3,k4))
    _finite(*result)
    return result


def floquet_monodromy(omega,steps=1024,r0=.5,epsilon=.1):
    """RK4 monodromy for eta''+[1-r0-epsilon*cos(omega*tau)]eta=0."""
    _finite(omega,r0,epsilon)
    if omega <= 0 or not isinstance(steps,int) or steps < 4:
        raise ValueError("Positive drive frequency and integer steps >= 4 required")
    h=math.tau/omega/steps; y=(1.,0.,0.,1.)
    def rhs(t,y):
        stiffness=1-r0-epsilon*math.cos(omega*t)
        return y[2],y[3],-stiffness*y[0],-stiffness*y[1]
    for i in range(steps): y=rk4_step(rhs,i*h,y,h)
    return ((y[0],y[1]),(y[2],y[3]))


def floquet_diagnostics(matrix,period,tolerance=1e-7):
    """Finite precision classification; repeated multipliers are inconclusive."""
    if len(matrix)!=2 or any(len(row)!=2 for row in matrix): raise ValueError("2x2 matrix required")
    a,b=matrix[0]; c,d=matrix[1]; _finite(a,b,c,d,period,tolerance)
    if period <= 0 or tolerance <= 0: raise ValueError("Positive period and tolerance required")
    trace=a+d; determinant=a*d-b*c
    root=cmath.sqrt(trace*trace-4*determinant)
    multipliers=((trace+root)/2,(trace-root)/2)
    radius=max(abs(x) for x in multipliers)
    if abs(determinant-1)>tolerance: classification='unresolved: determinant error'
    elif abs(abs(trace)-2)<=tolerance: classification='boundary: Jordan analysis required'
    elif abs(trace)>2: classification='unstable'
    else: classification='bounded'
    return dict(multipliers=multipliers,spectral_radius=radius,
                growth_rate=math.log(radius)/period,determinant=determinant,
                classification=classification)


def dna_origin_energies(q,velocity,alpha,penalty,mismatch,mismatch_rate=0.):
    """DNA norm penalty includes Origin. Modified energy is conserved only
    for constant mismatch. The metric and penalty coefficient are constant.
    """
    _finite(q,velocity,alpha,penalty,mismatch,mismatch_rate)
    if penalty < 0: raise ValueError("Nonnegative penalty required")
    base=origin_energy(q,velocity,alpha)
    shift=2*penalty*mismatch*mismatch
    return dict(base=base,modified=base+alpha*shift*q*q/2,
                base_rate=-alpha*shift*q*velocity,
                modified_rate=2*alpha*penalty*mismatch*mismatch_rate*q*q,
                acceleration=-(1+shift)*q)


def dna_ratio_mismatch(kappa,sigma):
    _finite(kappa,sigma)
    if sigma == 0: raise ValueError("Helical ratio is undefined at zero torsion")
    return kappa/sigma-math.pi/((1+math.sqrt(5))/2)


def dna_linear_spectrum(a,b,penalty,mismatch):
    _finite(penalty,mismatch)
    if penalty < 0: raise ValueError("Nonnegative penalty required")
    shift=2*penalty*mismatch*mismatch
    return tuple(lam+shift for lam in linear_spectrum(a,b))
