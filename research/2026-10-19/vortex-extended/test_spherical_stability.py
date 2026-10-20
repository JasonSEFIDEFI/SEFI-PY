import unittest
from spherical_stability import *

class CurvatureTests(unittest.TestCase):
    def test_real_fixed_charge_hessian(self):
        self.assertLess(check_hessian(),1e-7)

    def test_complex_phase_hessians(self):
        m=Radial(h=.4);f,s=m.initial();_,hp,hs=operators(m,f,s,0)
        rng=np.random.default_rng(18);root=np.sqrt(m.V)
        for which,op in [('phi',hp),('sigma',hs)]:
            q=rng.normal(size=len(f));q/=np.linalg.norm(q);d=q/root
            def gradient(t):
                F=f.astype(complex);S=s.astype(complex)
                if which=='phi':F+=1j*t*d
                else:S+=1j*t*d
                w=m.Q/np.sum(m.V*abs(F)**2)
                gf=m.K(F)/m.V+(1+abs(F)**2+4*abs(S)**2-w*w)*F
                gs=(m.K(S)-m.b)/m.V+(-3+4*abs(F)**2+m.lam*abs(S)**2)*S
                return root*np.imag(gf if which=='phi' else gs)
            finite=(gradient(1e-4)-gradient(-1e-4))/(2e-4)
            np.testing.assert_allclose(finite,op@q,rtol=1e-7,atol=1e-8)

    def test_higher_angular_sector_increment(self):
        m=Radial(h=.4);f,s=m.initial();a=operators(m,f,s,1);b=operators(m,f,s,2)
        rng=np.random.default_rng(5)
        for op1,op2 in zip(a,b):
            q=rng.normal(size=op1.shape[0]);diag=4/m.r**2
            if len(q)==2*len(f):diag=np.r_[diag,diag]
            np.testing.assert_allclose(op2@q-op1@q,diag*q,rtol=1e-10,atol=1e-10)

if __name__=='__main__':unittest.main()
