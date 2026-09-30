"""Independent dense oracle and block-inverse checks for the broad search."""
import unittest
import numpy as np
from scipy import sparse as sp
from scipy.linalg import eigvals
from scipy.sparse.linalg import eigs
from full_z_spectral import spectral_operator
from broad_scan import weights, structure_bound, schur_inverse, target_grid


class BroadChecks(unittest.TestCase):
    def test_weighted_structure_and_dense_growth_bound(self):
        n,h = 4,.5
        rng = np.random.default_rng(12)
        fields = [rng.normal(size=(n+1,2*n+1))*.2 for _ in range(3)]
        fields[0] += 1
        for m in (0,1,2,3,4):
            a,b,g,_,_ = spectral_operator(fields,n,h,.6,.3,25,m)
            bound = structure_bound(b,g,weights(n,h,m))
            ev = eigvals(a.toarray())
            growing = ev[abs(ev.real)>1e-7]
            if growing.size:
                self.assertLessEqual(max(abs(growing)**2),bound['lambda_max_B_estimate']+1e-8)
            self.assertLess(bound['B_hermitian_defect'],1e-12)
            self.assertLess(bound['G_skew_defect'],1e-12)

    def test_schur_inverse_against_dense_block_solve(self):
        n,h = 4,.5
        rng = np.random.default_rng(31)
        fields = [rng.normal(size=(n+1,2*n+1))*.1 for _ in range(3)]
        fields[0] += 1
        a,b,g,_,_ = spectral_operator(fields,n,h,.6,.3,25,2)
        for shift in (.17+.31j,.4-.2j):
            rhs = rng.normal(size=a.shape[0])+1j*rng.normal(size=a.shape[0])
            got = schur_inverse(b,g,shift)@rhs
            expected = np.linalg.solve(a.toarray()-shift*np.eye(a.shape[0]),rhs)
            np.testing.assert_allclose(got,expected,rtol=1e-11,atol=1e-11)

    def test_zero_target_can_miss_growth_and_broad_targets_recover(self):
        # Near-zero neutral oscillations mask the pair +/-0.8+/-0.6j.
        b = sp.diags([-.0001,-.0004,1.,1.,-4.,-9.],dtype=complex,format='csc')
        g = sp.diags([0,0,1.2j,-1.2j,0,0],format='csc')
        a = sp.bmat([[None,sp.eye(6)],[b,g]],format='csc')
        near = eigs(a,k=4,sigma=1e-4,OPinv=schur_inverse(b,g,1e-4),return_eigenvectors=False)
        self.assertLess(max(near.real),1e-8)
        found = []
        for shift in target_grid(1.):
            found.extend(eigs(a,k=4,sigma=shift,OPinv=schur_inverse(b,g,shift),return_eigenvectors=False))
        self.assertTrue(any(abs(v-(.8+.6j))<1e-8 for v in found))
        right = eigs(a,k=4,which='LR',return_eigenvectors=False)
        self.assertGreater(max(right.real),.79)


if __name__ == '__main__':
    unittest.main()
