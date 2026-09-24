"""Physical limits, data boundaries, observability, and reproducibility checks."""
import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from workbench import Diagnostics, design, load_observations, response, save_observations, simulate, truth_for


class WorkbenchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Diagnostics()
        cls.model.calibrate(np.random.default_rng(100), runs=8)

    def test_resonant_pi_pulse_and_zero_drive(self):
        settings = np.array([[0, np.pi], [0, 2*np.pi]])
        np.testing.assert_allclose(response(1, 0, settings, 0), [1, 0], atol=1e-14)
        np.testing.assert_allclose(response(0, 0, settings), [.02, .02])

    def test_probabilities_bounded(self):
        p = response(np.linspace(.1, 2, 100), np.linspace(-3, 3, 100), design())
        self.assertTrue(np.all((p >= .02) & (p <= .98)))

    def test_signed_detuning_observable_with_offset_scan(self):
        for a, d in ((.04, .06), (-.03, -.08)):
            data = response(1+a, d, design())[None, :]
            estimate, _, _ = self.model.estimate(data)
            np.testing.assert_allclose(estimate[0, 0], [a, d], atol=1e-8)

    def test_resonant_only_scan_abstains(self):
        model = Diagnostics(settings=np.array([[0, 1.1], [0, 2.2], [0, 3.3]]))
        model.thresholds = np.ones(3)
        model.residual_limit = 10
        np.testing.assert_allclose(response(1, .1, model.settings), response(1, -.1, model.settings))
        self.assertEqual(model.diagnose(model.reference[None, :])["status"], "insufficient evidence")

    def test_tangent_matches_small_physical_perturbation(self):
        shift = np.array([1e-4, -2e-4])
        z = (response(1+shift[0], shift[1], design())-self.model.reference)/self.model.sigma
        np.testing.assert_allclose(z @ self.model.projector.T, shift, atol=1e-7)

    def test_recovery_requires_three_complete_acceptable_scans(self):
        data = np.tile(self.model.reference, (4, 1))
        result = self.model.diagnose(data)
        self.assertEqual(result["recovery"].tolist(), [False, False, True, True])
        faulty = np.tile(response(1.08, .12, design()), (5, 1))
        self.assertFalse(self.model.diagnose(faulty)["recovery"].any())

    def test_mismatch_is_not_confident_diagnosis(self):
        result = self.model.diagnose(np.zeros((4, 12)))
        self.assertTrue(result["ambiguous"].all())
        self.assertFalse(result["recovery"].any())

    def test_replay_roundtrip_has_no_truth_columns(self):
        observed = simulate(truth_for("mixed"), self.model, np.random.default_rng(4))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"scan.csv"
            save_observations(path, observed, self.model)
            np.testing.assert_array_equal(load_observations(path, self.model), observed)
            header = path.read_text().splitlines()[0]
            self.assertNotIn("truth", header)
            self.assertNotIn("scenario", header)
            before = self.model.diagnose(observed)["estimates"]
            after = self.model.diagnose(load_observations(path, self.model))["estimates"]
            np.testing.assert_array_equal(before, after)

    def test_bad_replay_rejected(self):
        observed = simulate(truth_for("healthy"), self.model, np.random.default_rng(4))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"scan.csv"
            save_observations(path, observed, self.model)
            original = path.read_text()
            for broken in (original+original.splitlines()[1]+"\n", "\n".join(original.splitlines()[:-1])+"\n",
                           original.replace("1024", "0")):
                path.write_text(broken)
                with self.assertRaises(ValueError):
                    load_observations(path, self.model)

    def test_invalid_probabilities_and_config_rejected(self):
        for value in (np.nan, np.inf, -1, 2):
            with self.assertRaises(ValueError):
                self.model.estimate(np.full((2, 12), value))
        with self.assertRaises(ValueError):
            Diagnostics(shots=0)
        with self.assertRaises(ValueError):
            Diagnostics(alpha=0)

    def test_reproducibility_and_heldout_seed_independence(self):
        truth = truth_for("mixed")
        a = simulate(truth, self.model, np.random.default_rng(1))
        b = simulate(truth, self.model, np.random.default_rng(1))
        c = simulate(truth, self.model, np.random.default_rng(2))
        np.testing.assert_array_equal(a, b)
        self.assertFalse(np.array_equal(a, c))


if __name__ == "__main__":
    unittest.main(verbosity=2)
