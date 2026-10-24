"""Opt-in SEFI-PY integration with the standalone sefi-quantum-benchmark package.

No dependency is imported until a function is called. Existing dynamics and
repository module paths are not modified.
"""


def seeded_state(state_class, configs, seed=42):
    """Build a Haar state using the actual caller-supplied QuantumState class."""
    import numpy as np
    from sefi_benchmark.states import haar
    from sefi_benchmark.adapters import to_sefi_py
    count = len(configs)
    if count < 2 or count & (count-1):
        raise ValueError("configuration count must be a power of two >= 2")
    return to_sefi_py(haar(count.bit_length()-1, np.random.default_rng(seed)), state_class, configs)


def generate_qec(output, **kwargs):
    from sefi_benchmark.datasets import export_qec
    return export_qec(output, **kwargs)
