# SEFI QEC Stack — Research Software

**Status: computational research prototype.** This project explores software components for quantum-error-correction workflows and proposed geometric diagnostics. It does **not** establish fault-tolerant quantum computation, hardware-level correction, physical error suppression, or real-time decoder performance.

## Scope

The code is organized around stabilizer/parity evaluation, syndrome extraction and transport, decoder experiments, Pauli-frame bookkeeping, and an exploratory geometric residual-analysis module. The geometric terminology is part of a proposed research interpretation; it is not an experimentally validated quantum observable or demonstrated correction mechanism.

The project is not evidence that SEFI, GWFM, or DEFI is a validated physical theory.

## Modules

The documented software organization includes:

- `qec_stack/stabilizers/`: stabilizer and parity utilities
- `qec_stack/syndromes/`: syndrome extraction and transport
- `qec_stack/decoders/`: majority-vote and energy-minimization approaches
- `qec_stack/pauli_frame/`: Pauli-frame routing
- `qec_stack/warp_geometry/`: exploratory residual analysis
- `qec_stack/stack.py`: pipeline orchestration

Module names and interfaces should be checked against the checkout before use. The repository contains historical copies of this README; their presence does not establish independent implementations.

## Using and evaluating the code

Run from the directory containing the `qec_stack` package in an appropriately configured Python environment. Inspect `stack.py` and the available test files for the actual supported API before invoking the orchestrator. Earlier documentation included an unverified `engine.run()` example; it is intentionally not repeated here.

Any claimed improvement must be measured against a matched conventional baseline, with explicit error model, code parameters, decoding metric, sample count, random seed, and uncertainty. Software unit tests establish only the tested software properties.

## Evidence boundaries

| Claim | Current documentation status |
| --- | --- |
| Source code organized into QEC-related modules | Software architecture described here |
| Correctness of a specific decoder or pipeline | Requires reproducible tests for that implementation |
| Lower logical error rates than a baseline | Not established by this README |
| Real-time or low-latency operation | Requires timing benchmarks on specified hardware |
| Physical quantum-error-correction performance | Requires appropriate experimental evidence |
| SEFI-specific physical mechanism | Hypothesis, not demonstrated |

**Research principle:** Keep mathematical definitions, executable behavior, numerical comparisons, and physical claims separate.
