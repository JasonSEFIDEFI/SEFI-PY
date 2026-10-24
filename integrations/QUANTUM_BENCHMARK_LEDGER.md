# Quantum benchmark integration record — 2026-10-23

Base: SEFI-PY `847a071351ed2b1c54e240fd9a520a1cce238178`.
Reference implementation: SEFI_QEC_Stack
`aaad1aa1adfe358f468cfd7c371fb1f29259fcf5`, `benchmark_pipeline/`.

Added an opt-in wrapper, a real QuantumState round-trip example, documentation
and CI pinned to the reference implementation commit. Existing code, dynamics,
licensing and history are preserved. The example passed locally on Python
3.14.7 and NumPy 2.5.4 after installing sefi-quantum-benchmark 0.1.0.

An initial direct `core.quantum` import failed on the existing circular import.
The example now initializes `core.field_sovereignty` first, without changing
or mocking the legacy class. The successful example checks normalization and
exact amplitude round-trip. This is an integration check, not evidence for
SEFI physics or a claim that the complete SEFI-PY test suite passes.

The reference package passed 37 new tests and 25 existing Stack tests; four
10,000-test-trial seeded QEC examples and five noisy-state export kinds ran.
Its public research ledger includes baseline counts, uncertainty, checksums,
failed attempts, correction steps and model limitations. Published seeded
examples are reproducibility demos, not blind evaluations. No geometric QEC
advantage was measured. CI status must be checked independently.
