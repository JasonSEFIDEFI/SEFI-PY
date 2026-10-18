"""Paired, seed-cluster percentile bootstrap. Jason D. Dutton."""
import numpy as np


def interval(seed_means, replicates=10000, seed=927123):
    values = np.asarray(seed_means, float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("finite means from at least two independent seeds required")
    rng = np.random.default_rng(seed)
    bootstrap = values[rng.integers(0, len(values), (replicates, len(values)))].mean(axis=1)
    lo, hi = np.quantile(bootstrap, [0.025, 0.975])
    return float(values.mean()), float(lo), float(hi)
