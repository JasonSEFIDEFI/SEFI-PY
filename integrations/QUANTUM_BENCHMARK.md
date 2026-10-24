# Standard quantum benchmark integration

The standalone MIT-licensed pipeline lives in
[SEFI_QEC_Stack/benchmark_pipeline](https://github.com/JasonSEFIDEFI/SEFI_QEC_Stack/tree/feat/quantum-benchmark-v1/benchmark_pipeline).
It supplies seeded Haar/Ginibre/Bell/GHZ states, valid CPTP channels, small
circuits and ideal syndrome-only repetition-code datasets. These are standard
quantum mechanics references; no new SEFI physics or QEC advantage is claimed.

Install from a sibling checkout, then run from this repository root:

```sh
python -m pip install -e ../SEFI_QEC_Stack/benchmark_pipeline
python integrations/quantum_benchmark_example.py
sefi-benchmark qec runs/example --seed 42 --code bit3 --shots 1000
sefi-benchmark benchmark runs/example
```

The adapter takes the actual QuantumState class explicitly, validates its input
and copies amplitudes. It does not alter existing simulation behavior. The
legacy module has an import-order dependency; the example initializes
`core.field_sovereignty` before importing `core.quantum` and records any failure
honestly in the research ledger. Install the benchmark's NumPy dependency
before running. See the pipeline README for API, noise conventions, data
schema, metrics, limitations and independent validity tests.

QEC decoder inputs contain only code and syndrome. Seeds, injected errors and
test labels belong to the trusted evaluator; never expose its folder to an
untrusted decoder. Public seeded examples are reproducibility demos, not blind
evaluations. Training uses a separate random stream and only training labels.

Optional SEFI-inspired state descriptors are explicitly labeled conventional
purity and spectrum entropy, and never supplied to QEC decoders. Existing
legacy geometric identity bits are not treated as physical stabilizer data.
