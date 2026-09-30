import unittest
import numpy as np
from full_z_spectral import operators
from branch_audit import physical,discrete_functional
from characteristic_audit import potential_hessian,coupled_phase_speeds,exact_coupled_frequencies


class AuditChecks(unittest.TestCase):
    def test_stress_generator_and_charge_algebra(self):
        rng=np.random.default_rng(18);n=8;h=.4
        fields=[rng.normal(size=(n+1,2*n+1))*.1 for _ in range(3)];fields[0]+=1
        d=physical(fields,n,h,.6,.3,25)
        self.assertLess(abs(d['generator_algebra_residual']),1e-10)
        self.assertLess(abs(d['charge_algebra_residual']),1e-12)

    def test_discrete_action_variation_matches_field_equations(self):
        rng=np.random.default_rng(19);n=8;h=.4;c=.6;nu=.3;lam=25
        fields=[rng.normal(size=(n+1,2*n+1))*.1 for _ in range(3)];fields[0]+=1
        for j,f in enumerate(fields):f[-1,:]=1 if j==0 else 0;f[:,0]=f[:,-1]=1 if j==0 else 0
        u,w,s=[f[:n,1:-1].ravel() for f in fields];lap,dz,bc=operators(n,h)
        db=np.zeros((n,2*n-1));db[:,0]=-1/(2*h);db[:,-1]=1/(2*h)
        a=1-u*u-w*w-4*s*s
        residual=[lap@u+bc-c*(dz@w)+a*u,lap@w+c*(dz@u+db.ravel())+a*w,
            lap@s+(3+nu*nu-4*(u*u+w*w)-lam*s*s)*s]
        r=np.arange(n,dtype=float);r[0]=1/8;weights=2*np.pi*np.repeat(r*h**3,2*n-1)
        directions=[rng.normal(size=(n,2*n-1)) for _ in range(3)]
        expected=-sum(np.dot(weights*f,d.ravel()) for f,d in zip(residual,directions))
        eps=1e-6;plus=[f.copy() for f in fields];minus=[f.copy() for f in fields]
        for fp,fm,d in zip(plus,minus,directions):fp[:n,1:-1]+=eps*d;fm[:n,1:-1]-=eps*d
        got=(discrete_functional(plus,n,h,c,nu,lam)['G_discrete']-discrete_functional(minus,n,h,c,nu,lam)['G_discrete'])/(2*eps)
        self.assertAlmostEqual(got,expected,delta=1e-7)

    def test_cartesian_potential_hessian(self):
        x=np.array([.6,-.2,.12,.04]);lam=25
        def force(v):
            p,s=v[:2],v[2:];return np.r_[(1+p@p+4*(s@s))*p,(lam*(s@s)+4*(p@p)-3)*s]
        eps=1e-6
        got=np.column_stack([(force(x+eps*np.eye(4)[i])-force(x-eps*np.eye(4)[i]))/(2*eps) for i in range(4)])
        np.testing.assert_allclose(got,potential_hessian(x[:2],x[2:]),atol=1e-9)

    def test_ir_speeds_and_uv_front_limit(self):
        speeds,_,_=coupled_phase_speeds(1.,.1)
        self.assertGreater(abs(speeds[1]-speeds[0]),.1)
        k=.001;freq=exact_coupled_frequencies(k,1.,.1)
        np.testing.assert_allclose((freq[:2]/k)**2,speeds,rtol=1e-5)
        k=10000.;np.testing.assert_allclose(exact_coupled_frequencies(k,1.,.1)/k,1.,atol=.001)


if __name__=='__main__':unittest.main()
