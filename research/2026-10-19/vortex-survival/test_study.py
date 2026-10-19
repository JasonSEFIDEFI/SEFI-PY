import unittest
from study import *
from azimuthal import tangent
from balance import energy,charge,charge_rate

class AuditTests(unittest.TestCase):
    def test_core_and_negative_control(self):
        r=np.linspace(0,5,31);z=np.linspace(-3,3,37);rr,zz=np.meshgrid(r,z,indexing='ij')
        found,_=cores((rr-2.31)+1j*(zz-.13),r,z)
        self.assertEqual(len(found),1);self.assertEqual(found[0]['winding'],1)
        found,_=cores(np.ones_like(rr,dtype=complex),r,z);self.assertEqual(found,[])

    def test_core_boundary_ambiguity(self):
        r=np.arange(5.);z=np.arange(5.);a=np.ones((5,5),complex);a[2,2]=0
        _,ambiguous=cores(a,r,z);self.assertEqual(ambiguous,4)

    def test_tangent_matches_nonlinear_directional_derivative(self):
        f,p=physics.load_profile(CONSOLE/'data/charged_ring_refined.npz');p['epsilon']=0
        e=physics.Evolution(f,p);rng=np.random.default_rng(572)
        d=rng.normal(size=e.y.shape)+1j*rng.normal(size=e.y.shape)
        d[:,0,:]=d[:,-1,:]=d[:,:,0]=d[:,:,-1]=0
        eps=1e-6
        finite=(e.rhs(e.y+eps*d)-e.rhs(e.y-eps*d))/(2*eps)
        # Add exact angular Laplacian absent from the axisymmetric RHS.
        for a in range(2):finite[a+2,1:-1,1:-1]-=4*d[a,1:-1,1:-1]/e.r[1:-1,None]**2
        finite[:,0,:]=0
        actual=tangent(e,d,2)
        np.testing.assert_allclose(actual,finite,atol=3e-8,rtol=3e-8)

    def test_energy_and_charge_balance(self):
        f,p=physics.load_profile(CONSOLE/'data/charged_ring_refined.npz');p['epsilon']=.01
        e=physics.Evolution(f,p);rng=np.random.default_rng(1)
        e.y[2:]=.001*(rng.normal(size=e.y[2:].shape)+1j*rng.normal(size=e.y[2:].shape))
        e.y[2:,-1,:]=e.y[2:,:,0]=e.y[2:,:,-1]=0
        rhs=e.rhs(e.y);h=1e-5
        edot=(energy(e,e.y+h*rhs)-energy(e,e.y-h*rhs))/(2*h)
        self.assertLess(abs(edot),1e-7)
        qdot=(charge(e,e.y+h*rhs)-charge(e,e.y-h*rhs))/(2*h)
        self.assertAlmostEqual(qdot,charge_rate(e,e.y),places=7)

if __name__=='__main__':unittest.main(verbosity=2)
