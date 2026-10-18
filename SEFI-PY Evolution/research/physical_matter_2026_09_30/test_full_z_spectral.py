"""Independent checks of linearization, regularity, and unrestricted parity."""
import unittest
import numpy as np
from full_z_spectral import operators,spectral_operator,load_profile,stationary
from pathlib import Path


class SpectralChecks(unittest.TestCase):
    def test_angular_axis(self):
        n,h=12,.4
        for ell in (0,1,2):
            lap,_,_=operators(n,h,ell)
            i0=0 if ell==0 else 1
            r=np.arange(i0,n)*h
            f=np.repeat((r**ell)[:,None],2*n-1,axis=1)
            residual=(lap@f.ravel()).reshape(f.shape)
            # r^ell is harmonic; finite difference error for ell=2 is exact here.
            self.assertLess(abs(residual[:-1,1:-1]).max(),2e-12)

    def test_bogoliubov_variation(self):
        n,h=8,.5
        rng=np.random.default_rng(7)
        fields=[rng.normal(size=(n+1,2*n+1))*.1 for _ in range(3)]
        fields[0]+=1
        full,b,g,modes,sl=spectral_operator(fields,n,h,.6,.3,25,0)
        p=(fields[0]+1j*fields[1])[sl].ravel(); s=fields[2][sl].ravel().astype(complex)
        xi=rng.normal(size=p.size)+1j*rng.normal(size=p.size)
        eta=rng.normal(size=p.size)+1j*rng.normal(size=p.size)
        lap,dz,_=operators(n,h)
        def force(p,s):
            return np.r_[lap@p+1j*.6*(dz@p)+(1-abs(p)**2-4*abs(s)**2)*p,
                         lap@s+(3+.3**2-4*abs(p)**2-25*abs(s)**2)*s]
        eps=1e-6
        fd=(force(p+eps*xi,s+eps*eta)-force(p-eps*xi,s-eps*eta))/(2*eps)
        y=np.r_[xi,xi.conj(),eta,eta.conj()]; out=b@y
        np.testing.assert_allclose(np.r_[out[:p.size],out[2*p.size:3*p.size]],fd,atol=2e-8,rtol=1e-8)
        np.testing.assert_allclose(out[p.size:2*p.size],out[:p.size].conj(),atol=1e-12)
        np.testing.assert_allclose(out[3*p.size:],out[2*p.size:3*p.size].conj(),atol=1e-12)
        # All Z interior nodes exist for each field, including opposite parity.
        self.assertEqual(full.shape,(8*n*(2*n-1),)*2)
        gv=.6/(2*np.sqrt(2)); gamma=np.sqrt(1+gv*gv)
        expected=np.r_[2*gv*(dz@xi)-2j*np.sqrt(2)*xi,
                       2*gv*(dz@xi.conj())+2j*np.sqrt(2)*xi.conj(),
                       2*gv*(dz@eta)-2j*.3*gamma*eta,
                       2*gv*(dz@eta.conj())+2j*.3*gamma*eta.conj()]
        np.testing.assert_allclose(g@y,expected,atol=1e-12)

    def test_retained_profile_phase_identity(self):
        source=Path(__file__).resolve().parents[2]/'sim/vortex_console/data/charged_ring_refined.npz'
        if not source.exists(): source=Path(__file__).with_name('reference_profiles')/'charged_ring_refined.npz'
        n,h,f,c,nu,lam,winding=load_profile(source)
        f,hist=stationary(f,n,h,c,nu,lam)
        self.assertLess(hist[-1],1e-9)
        _,b,_,modes,_=spectral_operator(f,n,h,c,nu,lam,0)
        y=modes['carrier_phase']
        self.assertLess(abs(b@y).max(),1e-9)


if __name__=='__main__': unittest.main()
