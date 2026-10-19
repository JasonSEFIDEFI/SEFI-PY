import unittest
from cylindrical3d import *
sys.path.insert(0,str(PRIOR))
from azimuthal import tangent

class ExtensionTests(unittest.TestCase):
    def test_axisymmetric_reduction(self):
        c=Cylinder(ntheta=8);rng=np.random.default_rng(3)
        y=c.axis.y.copy();d=.001*(rng.normal(size=y.shape)+1j*rng.normal(size=y.shape))
        d[:,-1,:]=d[:,:,0]=d[:,:,-1]=0;y+=d
        yy=np.repeat(y[:,:,None,:],8,axis=2)
        actual=c.rhs(yy);expected=c.axis.rhs(y)
        np.testing.assert_allclose(actual, np.repeat(expected[:,:,None,:],8,axis=2),atol=1e-12)

    def test_angular_laplacian(self):
        c=Cylinder(ntheta=16);rr,tt,zz=np.meshgrid(c.r,c.theta,c.z,indexing='ij')
        q=(rr**2*np.exp(2j*tt)).astype(complex)
        # r^2 exp(2i theta) is harmonic in the transverse plane.
        np.testing.assert_allclose(c.lap(q)[1:-1,:,1:-1],0,atol=1e-11)

    def test_tangent_reduction(self):
        c=Cylinder(ntheta=16);rr,zz=np.meshgrid(c.r,c.z,indexing='ij')
        g=(rr/c.axis.radius)**2*np.exp(-((rr-c.axis.radius)**2+zz**2)/2)
        g[0,:]=g[-1,:]=g[:,0]=g[:,-1]=0
        d=np.array([g.astype(complex),(.4+.3j)*g,np.zeros_like(g),np.zeros_like(g)])
        angular=np.cos(2*c.theta)[None,None,:,None];full=d[:,:,None,:]*angular
        h=1e-6;actual=(c.rhs(c.y+h*full)-c.rhs(c.y-h*full))/(2*h)
        expected=tangent(c.axis,d,2)[:,:,None,:]*angular
        np.testing.assert_allclose(actual,expected,atol=3e-8,rtol=3e-8)

    def test_energy_charge_differential(self):
        c=Cylinder(ntheta=8);c.perturb(.01)
        c.y[2]=.002j*(c.y[0]-c.base[0]);c.y[3]=.003*(c.y[1]-c.base[1])
        y=c.y.copy();rhs=c.rhs(y);h=1e-5
        c.y=y+h*rhs;ep=c.energy();qp=c.charge()
        c.y=y-h*rhs;em=c.energy();qm=c.charge()
        self.assertLess(abs((ep-em)/(2*h)),1e-7)
        self.assertLess(abs((qp-qm)/(2*h)),1e-7)

if __name__=='__main__':unittest.main(verbosity=2)
