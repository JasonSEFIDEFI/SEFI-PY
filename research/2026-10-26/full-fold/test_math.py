"""Checks for the new fold expansion and charge-constraint identities."""
import numpy as np
import sympy as sp

def test_fold_first_order_expansion():
    d,A,B,C,D=sp.symbols('d A B C D',nonzero=True)
    # Local a-a0=x. q0'(a0)=0, q0''(a0)=A; q1(a0)=C, q1'(a0)=D.
    x=sp.symbols('x')
    q0=A*x*x/2+B*x**3/6
    q1=C+D*x
    a1=-D/A
    derivative=sp.diff(q0+d*q1,x).subs(x,d*a1)
    assert sp.expand(derivative).coeff(d,1)==0
    assert sp.expand((q0+d*q1).subs(x,d*a1)).coeff(d,1)==C

def test_charge_chain_rule_and_rank_one_coefficient():
    # Differentiate Q=omega*I and the nonlocal charge term independently.
    c=sp.symbols('c')
    w=sp.Function('w')(c);I=sp.Function('I')(c)
    Q=w*I
    amplitude_rhs=2*w*sp.diff(w,c)+2*w*w*sp.diff(I,c)/I
    assert sp.simplify(amplitude_rhs-2*w*sp.diff(Q,c)/I)==0
    eps,q,i,a,b=sp.symbols('eps q i a b')
    charge=q*q/(2*(i+2*eps*a+eps*eps*b))
    exact=sp.diff(charge,eps,2).subs(eps,0)
    assert sp.simplify(exact-(-q*q*b/i**2+4*q*q*a*a/i**3))==0

def test_radial_family_null_direction_scaling():
    # I_scaled=4*pi*J. This catches a spurious 4*pi in the spectral RHS.
    w,J,wc,Jc,d=sp.symbols('w J wc Jc d',nonzero=True)
    Ic=4*sp.pi*Jc;I=4*sp.pi*J;Qc=wc*I+w*Ic
    assert sp.simplify(2*w*Qc/(d*I)-(2*w*wc/d+2*w*w*Jc/(d*J)))==0
