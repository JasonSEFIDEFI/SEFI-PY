# SEFI-PY

## Interactive vortex research console

The [Vortex Research Console](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/main/SEFI-PY%20Evolution/sim/vortex_console) provides a locally run 3D view of computed field profiles, adjustable viewing angles, actual profile recalculation, axisymmetric evolution, saved-frame playback, independent parameter sweeps, and checkpoints. It has labels for general readers and specialists, plus mathematical explanations and explicit limits.

- [Start and run the console](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/main/SEFI-PY%20Evolution/sim/vortex_console/README.md)
- [Reviewer questions and GitHub Copilot guide](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/main/SEFI-PY%20Evolution/sim/vortex_console/REVIEWER_GUIDE.md)
- [Executed tests and numerical evidence](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/main/SEFI-PY%20Evolution/sim/vortex_console/TEST_REPORT.md)
- [How future mathematical findings enter the model](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/main/SEFI-PY%20Evolution/sim/vortex_console/MODEL_HISTORY.md)

Reviewers can use their own GitHub Copilot access with the SEFI-PY repository as context, or read the same questions without Copilot. No personal chat history is shared. This is a classical field research prototype: physical matter and universal spacetime geometry remain open goals. Existing engine and QEC modules are unchanged.

A modular scientific engine implementing the Single Entity Field Interpretation (SEFI), Dynamic Entity Field Integration (DEFI), and the Geometric Waveform Model (GWFM).  
SEFI-PY provides a unified geometric framework for worldlines, warp modes, curvature, torsion, collapse behavior, and quantum error correction.

---

## Core Architecture

### SEFI (Single Entity Field Interpretation)
A geometric field model describing:
- Worldline stability
- Field origin, authorship, and sovereignty layers
- Warp modes (tangent, normal, binormal)
- Collapse geometry and measurement behavior

### DEFI (Dynamic Entity Field Integration)
A dynamic realignment model used for:
- Correction
- Stabilization
- Worldline restoration
- Error integration and geometric consistency

### GWFM (Geometric Waveform Model)
A waveform interpretation built on geometric invariants.

---

## SEFI-QEC: Quantum Error Correction Subsystem

SEFI-QEC integrates quantum error correction directly into the geometric field engine.  
Logical qubits are represented as SEFI worldlines, and physical qubits are geometric samples of that worldline.

### Features
- Real Pauli operations (X, Z, Y)
- Stabilizer parity checks (Z1Z2, Z2Z3)
- Syndrome extraction
- DEFI-based correction
- Warp-mode → error-mode mapping
- Warp-residual geometric alignment
- Pauli-frame consistency checking
- Majority-vote logical consistency
- Stabilizer-energy minimization decoder
- Multi-angle QEC benchmarking

### QEC Checking Suite (5 independent decoders)
1. **Stabilizer Parity Check**  
   Classical repetition-code stabilizers.

2. **Majority Vote Check**  
   Independent logical consistency decoder.

3. **Pauli-Frame Consistency**  
   Frame-based mismatch detection.

4. **Warp-Residual Geometric Check**  
   SEFI-native geometric deviation analysis.

5. **Stabilizer-Energy Minimization**  
   Physics-inspired energy-based decoder.

These decoders run independently and can be cross-validated, providing a multi-angle correction suite.

---

## Benchmarking

SEFI-QEC includes a benchmarking module that:
- Injects random X/Z/Y errors
- Runs all decoders
- Applies DEFI correction
- Measures recovery success rate

The full suite currently passes at **100%**.

---

## Project Status
SEFI-PY remains fully green across all modules:
- SEFI core
- DEFI integration
- GWFM
- Warp simulation
- Quantum warp simulation
- SEFI-QEC subsystem
- Multi-decoder QEC suite
- Benchmarking

All tests pass at 100%.

---

## Running the QEC Demo

```python
from sefi_qec.runner import run_sefi_qec_demo
run_sefi_qec_demo()
```

## Isolated field research checkpoint

The [23 September 2026 research addition](SEFI-PY%20Evolution/research/physical_matter_2026_09_23/README.md) contains coupled vortex/carrier solvers, reference profiles, nonlinear and azimuthal perturbation probes, an analytical rest-state audit, and the research-direction note.

It is opt-in and separate from the existing engine: no engine modules, entry points, configuration or dependency files are changed. The goal of deriving physical matter remains unachieved. Finite-domain classical candidates and limited disturbance tests must not be interpreted as observed particles or universal gravity. The two unsuccessful revised rest-loop searches are retained in the evidence.

See the addition's README for reproduction commands, model restrictions and validation limits.
