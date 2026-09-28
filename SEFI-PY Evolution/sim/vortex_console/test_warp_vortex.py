import unittest
import numpy as np
import warp_vortex as w

class WarpTests(unittest.TestCase):
    def setUp(self):
        self.axes=(np.linspace(-1,1,7), np.linspace(-1,1,8), np.linspace(-1,1,9))
        self.x,self.y,self.z=np.meshgrid(*self.axes,indexing='ij')
        self.D=np.stack((-self.y,self.x,np.zeros_like(self.x)),axis=-1)
        self.shape=self.x.shape

    def status(self, level=None, weight=None):
        return w.warp_vortex_status(self.D,self.axes,np.ones(self.shape),np.zeros(self.shape),{},
            lambda *a:np.zeros(self.shape+(1,)),lambda *a:np.full(self.shape,2.) if weight is None else weight,
            lambda *a:np.full(self.shape,-1.) if level is None else level)

    def test_rotation_and_integral(self):
        r=self.status()
        np.testing.assert_allclose(r['omega_D'],np.broadcast_to([0,0,2],self.D.shape),atol=1e-12)
        np.testing.assert_allclose(r['hessian'],0,atol=1e-12)
        self.assertAlmostEqual(r['weighted_vortex'],64)
        self.assertAlmostEqual(r['envelope_volume'],8)
        self.assertAlmostEqual(r['rms_envelope'],2)

    def test_hessian_nonuniform(self):
        axes=(np.array([-1,-.6,-.1,.3,1]),)*3
        x,y,z=np.meshgrid(*axes,indexing='ij')
        G,H=w.displacement_derivatives(np.stack((x*x,y*y,z*z),-1),axes)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    np.testing.assert_allclose(H[...,i,j,k],2 if i==j==k else 0,atol=1e-11)
        np.testing.assert_allclose(w.curl_displacement(G),0,atol=1e-12)

    def test_envelope_and_sigma(self):
        r=self.status(level=self.x)
        self.assertTrue(np.all(r['omega_envelope'][self.x>0]==0))
        self.assertAlmostEqual(r['weighted_vortex'],4*r['envelope_volume']*2)
        G,H=w.displacement_derivatives(self.D,self.axes)
        m=w.check_stability_manifold(self.D,G,H,np.ones(self.shape),np.zeros(self.shape),{},
            lambda *a:np.stack((self.x,np.zeros(self.shape)),-1),.01)
        np.testing.assert_array_equal(m,abs(self.x)<=.01)
        r=w.warp_vortex_status(self.D,self.axes,np.ones(self.shape),np.zeros(self.shape),{},
            lambda *a:np.ones(self.shape+(1,)),lambda *a:np.ones(self.shape),lambda *a:-np.ones(self.shape))
        self.assertEqual(r['weighted_vortex'],0)

    def test_empty_and_invalid(self):
        r=self.status(level=np.ones(self.shape))
        self.assertTrue(np.isnan(r['rms_envelope']))
        self.assertEqual(r['weighted_vortex'],0)
        with self.assertRaises(ValueError):self.status(weight=-np.ones(self.shape))
        with self.assertRaises(ValueError):w.displacement_derivatives(self.D,(self.axes[0][::-1],*self.axes[1:]))
        with self.assertRaises(ValueError):self.status(level=np.array(0))

    def test_worldline(self):
        s=np.linspace(0,2*np.pi,1001); zero=np.zeros_like(s); one=np.ones_like(s)
        x=np.stack((np.cos(s),np.sin(s),zero),-1)
        T=np.stack((-np.sin(s),np.cos(s),zero),-1); N=-x; B=np.stack((zero,zero,one),-1)
        G=np.broadcast_to(np.array([[0,-1,0],[1,0,0],[0,0,0.]]),(len(s),3,3))
        r=w.worldline_warp_vortex_profile(s,x,one,zero,zero,T,N,B,one,zero,{'a':2},G)
        np.testing.assert_allclose(r['derivative'],N,atol=1e-12)
        np.testing.assert_allclose(r['derivative_residual'],0,atol=1e-12)
        self.assertLess(np.max(abs(r['tangent_residual'])),2e-5)

if __name__=='__main__':unittest.main()
