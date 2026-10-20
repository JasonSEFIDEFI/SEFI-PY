import unittest
from run_integrity import eligible

class IntegrityTests(unittest.TestCase):
    def case(self):return dict(reason='completed',T_requested=12,samples=[dict(time=11.9999999999)])
    def test_complete_roundoff(self):self.assertTrue(eligible(self.case()))
    def test_partial(self):
        r=self.case();r['samples'][-1]['time']=8;self.assertFalse(eligible(r))
    def test_stops(self):
        for reason in ['user_interrupted','running','wall_limit','nonfinite','amplitude_guard']:
            r=self.case();r['reason']=reason;self.assertFalse(eligible(r))
    def test_invalid(self):
        r=self.case();r['samples'][-1]['time']=float('nan');self.assertFalse(eligible(r))

if __name__=='__main__':unittest.main(verbosity=2)
