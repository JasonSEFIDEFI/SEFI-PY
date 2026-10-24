"""Run from repository root after installing the standalone benchmark package."""
import importlib.util
from pathlib import Path
import sys
import numpy as np

root = Path(__file__).resolve().parents[1]
# Resolve this repository's legacy imports explicitly in this example only.
sys.path.insert(0, str(root/"SEFI-PY"))
# Legacy modules have an order-dependent circular import. Initialize the
# defining module first rather than changing or mocking the existing code.
import core.field_sovereignty
from core.quantum import QuantumState

spec = importlib.util.spec_from_file_location("benchmark_integration", Path(__file__).with_name("quantum_benchmark.py"))
integration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(integration)
state = integration.seeded_state(QuantumState, ["basis_0", "basis_1"], seed=42)
assert np.isclose(np.linalg.norm(state.amplitudes), 1)
from sefi_benchmark.adapters import from_sefi_py
assert np.array_equal(from_sefi_py(state), state.amplitudes)
print("SEFI-PY QuantumState adapter: normalized state and round-trip passed")
