import json,threading,unittest,urllib.request
import server

class EvidenceTests(unittest.TestCase):
    def test_read_only_evidence_routes(self):
        old_port=server.PORT;http=server.ThreadingHTTPServer(('127.0.0.1',0),server.Handler)
        server.PORT=http.server_port;t=threading.Thread(target=http.serve_forever,daemon=True);t.start()
        try:
            base=f'http://127.0.0.1:{server.PORT}'
            with urllib.request.urlopen(base+'/api/research-evidence',timeout=10) as response:data=json.load(response)
            self.assertEqual(len(data['axisymmetric']),8)
            self.assertIn('No independent Authorship',data['scope'])
            for route,expected in [('/evidence','Vortex research evidence'),('/evidence.js','/api/research-evidence')]:
                with urllib.request.urlopen(base+route,timeout=10) as response:self.assertIn(expected,response.read().decode())
        finally:http.shutdown();http.server_close();t.join();server.PORT=old_port

if __name__=='__main__':unittest.main(verbosity=2)
