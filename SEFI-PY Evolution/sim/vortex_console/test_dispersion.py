"""Analytic limits, independent eigensystem and live HTTP regression checks."""
import json
import threading
import unittest
import urllib.request
import numpy as np
import dispersion
import server

class DispersionTests(unittest.TestCase):
    def test_static_limit(self):
        k=np.linspace(0,3,41); low,high=dispersion.frequencies(k,0)
        np.testing.assert_allclose(low,k,atol=1e-14)
        np.testing.assert_allclose(high**2,k*k+2,atol=1e-14)

    def test_rotating_limits(self):
        low,high=dispersion.frequencies(np.array([0.,1e-7,1e5]),np.sqrt(2))
        self.assertEqual(low[0],0)
        self.assertAlmostEqual(high[0],np.sqrt(10))
        self.assertAlmostEqual(low[1]/1e-7,1/np.sqrt(5),places=12)
        self.assertAlmostEqual(low[2]/1e5,1,places=4)

    def test_independent_eigenvalues(self):
        for omega in [0,.3,np.sqrt(2)]:
            for k in [.05,.25,1,3]:
                low,high=dispersion.frequencies(k,omega)
                actual=np.sort(np.imag(np.linalg.eigvals(dispersion.linear_matrix(k,omega))))
                np.testing.assert_allclose(actual,[-high,-low,low,high],atol=1e-12)

    def test_invalid_inputs(self):
        for k,omega in [(-1,0),(float('nan'),0),(1,-1),(1,float('inf'))]:
            with self.assertRaises(ValueError):dispersion.frequencies(k,omega)

    def test_http_and_page(self):
        old_port=server.PORT
        http=server.ThreadingHTTPServer(('127.0.0.1',0),server.Handler)
        server.PORT=http.server_port
        thread=threading.Thread(target=http.serve_forever,daemon=True);thread.start()
        try:
            base=f'http://127.0.0.1:{server.PORT}'
            with urllib.request.urlopen(base+'/api/dispersion',timeout=10) as response:
                data=json.load(response)
            self.assertEqual(data['model_version'],dispersion.MODEL_VERSION)
            self.assertEqual(len(data['k']),241)
            self.assertAlmostEqual(data['speed_reduction_percent'],55.27864045)
            with urllib.request.urlopen(base+'/comparison',timeout=10) as response:
                self.assertIn('Background wave comparison',response.read().decode())
            with urllib.request.urlopen(base+'/dispersion.js',timeout=10) as response:
                self.assertIn('/api/dispersion',response.read().decode())
        finally:
            http.shutdown();http.server_close();thread.join();server.PORT=old_port

if __name__=='__main__':unittest.main(verbosity=2)
