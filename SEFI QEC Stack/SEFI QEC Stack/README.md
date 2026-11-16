# SEFI QEC Stack — Architecture Notes and Research Archive

**Status: historical architecture notes in this directory; the implementation is maintained in a separate repository.**

This folder preserves a QEC software design associated with SEFI/GWFM research. Its earlier README described a complete modular engine and provided a `QECStack().run()` example. The cited `qec_stack/stack.py` path is not present in this directory. The actual Python orchestrator was located at [`JasonSEFIDEFI/SEFI_QEC_Stack/stack.py`](https://github.com/JasonSEFIDEFI/SEFI_QEC_Stack/blob/main/stack.py). The historical `QECStack().run()` example does not match that source. See the [authoritative QEC README](https://github.com/JasonSEFIDEFI/SEFI_QEC_Stack) for the real API and scientific limitations.

## Proposed architecture (not an implementation inventory)

The design calls for stabilizer/parity handling, syndrome transport, decoder comparison, Pauli-frame bookkeeping, and exploratory geometric residual analysis. These are *intended components*. Their existence, correctness, integration, and performance must be established from actual source files and reproducible tests before any operational claim.

The geometric vocabulary is an exploratory interpretation, not a demonstrated physical error-correction mechanism.

## What would be required to claim an executable QEC stack

1. Use the separate [SEFI_QEC_Stack repository](https://github.com/JasonSEFIDEFI/SEFI_QEC_Stack) as the authoritative implementation and verify its dependencies.
2. Document the actual callable interfaces and installation requirements.
3. Add deterministic tests with explicit stabilizer definitions, noise models, expected syndromes, and Pauli-frame updates.
4. Compare decoder outcomes with conventional baselines under matched information and seeds.
5. Report logical-error metrics, uncertainties, code parameters, and any timing measurements on named hardware.

No real-time performance, quantum advantage, fault tolerance, hardware operation, or physical SEFI validation is established by this folder.

## Historical documentation

The original architecture narrative remains accessible in Git history. A separate `sefi_qec_stack.tex` document also describes the earlier design; it should be read as a proposal unless its described source and tests are independently verified.

For the broader research program, see the [repository README](../../README.md). The nested copy of this README provides the same status warning.
