import math
import random
from fractions import Fraction
import pytest
from core.time_dynamics import (
    IdentityMetric, Coupling, IdentityModel, origin_energy, clock_phase,
    phase_history, relational_derivative, linear_parameters, linear_spectrum,
    linear_modes, linear_solution, recurrence_cycles, exact_linear_ratios,
    locking_residual, periodic_descent_residual, rk4_step,
    floquet_monodromy, floquet_diagnostics, dna_origin_energies,
    dna_ratio_mismatch, dna_linear_spectrum,
)


def model(mu=Coupling.linear(.5)):
    return IdentityModel(IdentityMetric(),mu,Coupling.linear(0))


def test_base_energy_derivative_and_driven_work():
    m=IdentityModel(IdentityMetric(2,3,4,5),Coupling.quadratic(.2),Coupling.linear(-.4))
    y=(.3,.7,-.2,.8,-.5,.4,.3,-.1); k,s=1.2,.9; kd,sd=.13,-.2
    rate=m.rhs(y,k,s); h=1e-6
    plus=tuple(a+h*b for a,b in zip(y,rate)); minus=tuple(a-h*b for a,b in zip(y,rate))
    derivative=(m.energy(plus,k+h*kd,s+h*sd)-m.energy(minus,k-h*kd,s-h*sd))/(2*h)
    assert derivative==pytest.approx(m.prescribed_power(y,k,s,kd,sd),abs=2e-9)
    assert y[4]*(rate[4]+y[0])==0


def test_completed_square_and_nonlinear_saddle():
    m=model(Coupling.quadratic(.5))
    assert m.effective_potential(math.sqrt(2),1,0)==pytest.approx(.5)
    center=m.equilibrium(0,1,0); saddle=m.equilibrium(math.sqrt(2),1,0)
    assert saddle[2]==pytest.approx(-1)
    assert min(m.equilibrium_spectrum(0,1,0))>0
    assert min(m.equilibrium_spectrum(math.sqrt(2),1,0))<0
    assert max(abs(v) for v in m.rhs(saddle,1,0))<1e-14
    a=.8; state=(0,a,.4,-.7,0,0,0,0)
    completed=m.effective_potential(a,1,0)+(.4+.5*a*a)**2/2+(-.7)**2/2
    assert m.energy(state,1,0)==pytest.approx(completed)
    with pytest.raises(ValueError): m.equilibrium(.5,1,0)


def test_eigenbasis_and_exact_solution_against_integrator():
    rng=random.Random(42)
    for _ in range(100):
        a,b=rng.uniform(-2,2),rng.uniform(-2,2)
        modes=linear_modes(a,b)
        for lam,e in zip(linear_spectrum(a,b),modes):
            be=(e[0]+a*e[1]+b*e[2],a*e[0]+e[1],b*e[0]+e[2])
            assert be==pytest.approx(tuple(lam*x for x in e),abs=2e-15)
        for i,e in enumerate(modes):
            for j,f in enumerate(modes): assert sum(x*y for x,y in zip(e,f))==pytest.approx(float(i==j),abs=1e-15)
    for r in (0,.5,1,1.2):
        p=(.2,.7,-.3); v=(-.4,.1,.8); y=p+v
        def rhs(t,y): return y[3:]+(-y[0]-r*y[1],-r*y[0]-y[1],-y[2])
        for i in range(200): y=rk4_step(rhs,i*.01,y,.01)
        pp,vv=linear_solution(r,0,p,v,2)
        assert y==pytest.approx(pp+vv,abs=2e-9)


def test_clock_record_multiple_cycles_and_alias_rejection():
    times=[i*.2 for i in range(150)]; offset=6.1
    q=[math.cos(t+offset) for t in times]; v=[-math.sin(t+offset) for t in times]
    lifted=phase_history(times,q,v,initial_phase=offset)
    assert lifted==pytest.approx([t+offset for t in times],abs=2e-13)
    assert origin_energy(q[5],v[5],alpha=3)==pytest.approx(1.5)
    with pytest.raises(ValueError): clock_phase(0,0)
    with pytest.raises(ValueError): phase_history([0,math.tau],[1,1],[0,0])
    with pytest.raises(ValueError): phase_history([0,.2],[1,1],[0,0])
    with pytest.raises(ValueError): relational_derivative([1],0)
    assert relational_derivative((1,2),2)==(.5,1)


def test_recurrence_is_exact_certificate_not_float_fit():
    ratios=exact_linear_ratios(Fraction(24,25))
    assert ratios==(Fraction(1,5),Fraction(1),Fraction(7,5))
    assert recurrence_cycles(ratios)==5
    p,v=linear_solution(.96,0,(.2,.3,.4),(.5,.6,.7),10*math.pi)
    assert p+v==pytest.approx((.2,.3,.4,.5,.6,.7),abs=2e-14)
    with pytest.raises(ValueError): exact_linear_ratios(Fraction(1,2))
    with pytest.raises(ValueError): recurrence_cycles([.2])
    assert periodic_descent_residual(math.cos,.8)==pytest.approx(0,abs=1e-15)
    assert abs(periodic_descent_residual(lambda x: math.cos(x/5),.8))>.1
    assert locking_residual(7,3,2,1)==0


def test_floquet_counterexample_and_convergence():
    matrices=[floquet_monodromy(math.sqrt(2),n) for n in (128,256,512,1024)]
    errors=[max(abs(a-b) for ar,br in zip(matrices[i],matrices[i+1]) for a,b in zip(ar,br)) for i in range(3)]
    assert math.log2(errors[0]/errors[1])==pytest.approx(4,abs=.08)
    matrix=matrices[-1]
    assert matrix[0]==pytest.approx((-1.0123334577078713,-.2248503040205365),abs=2e-14)
    d=floquet_diagnostics(matrix,math.tau/math.sqrt(2))
    assert d['classification']=='unstable'
    assert d['growth_rate']==pytest.approx(.03531402295026381,abs=2e-14)
    assert d['determinant']==pytest.approx(1,abs=3e-14)
    stable=floquet_diagnostics(floquet_monodromy(2.2),math.tau/2.2)
    assert stable['classification']=='bounded'
    assert floquet_diagnostics(((1,1),(0,1)),1)['classification'].startswith('boundary')
    assert floquet_diagnostics(((2,0),(0,2)),1)['classification'].startswith('unresolved')


def test_dna_modified_energy_derivatives_and_static_threshold():
    q,v,alpha,penalty,delta,delta_rate=.7,-.4,2,.8,.3,.2
    e=dna_origin_energies(q,v,alpha,penalty,delta,delta_rate); h=1e-6
    plus=dna_origin_energies(q+h*v,v+h*e['acceleration'],alpha,penalty,delta+h*delta_rate)
    minus=dna_origin_energies(q-h*v,v-h*e['acceleration'],alpha,penalty,delta-h*delta_rate)
    assert (plus['modified']-minus['modified'])/(2*h)==pytest.approx(e['modified_rate'],abs=1e-10)
    assert (plus['base']-minus['base'])/(2*h)==pytest.approx(e['base_rate'],abs=1e-10)
    assert dna_origin_energies(q,v,alpha,penalty,0)['base_rate']==0
    assert min(dna_linear_spectrum(1.2,0,1,.5))>0
    assert dna_ratio_mismatch(21*math.pi,34)==pytest.approx(-.0012155762141235)
    with pytest.raises(ValueError): dna_ratio_mismatch(1,0)


@pytest.mark.parametrize('values',[(0,1,1,1),(-1,1,1,1),(1,float('nan'),1,1)])
def test_invalid_metrics(values):
    with pytest.raises(ValueError): IdentityMetric(*values)


def test_metric_normalization():
    a,b,r=linear_parameters(2,3,2,1,IdentityMetric(1,4,9,16))
    assert (a,b)==pytest.approx((8/6,3/8))
    assert r==pytest.approx(math.hypot(a,b))


def test_integrator_rejects_malformed_and_nonfinite_derivatives():
    with pytest.raises(ValueError): rk4_step(lambda t,y: (1,),0,(1,2),.1)
    with pytest.raises(ValueError): rk4_step(lambda t,y: (float('nan'),),0,(1,),.1)
    with pytest.raises(ValueError): Coupling(1,2,3)
