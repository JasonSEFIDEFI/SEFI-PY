# SEFI-PY — Python Research Implementation

**Status: exploratory mathematical and computational software.**

This directory contains Python implementations and demonstrations associated with the proposed Single Entity Field Interpretation (SEFI) and related geometric modeling. Its identity-layer terminology describes a **modeling framework**, not an experimentally verified structure of matter or consciousness.

## Conceptual components

The program uses proposed labels including **Field Origin**, **Field Authorship**, **Field Sovereignty**, **Warp Expression**, and **DEFI** to organize model states and transformations. Their names do not establish that the corresponding structures are independent physical observables, universally stable, or derived from established physics.

The software should be assessed module by module: implemented functions, tested behavior, mathematical assumptions, and physical interpretation are separate questions.

## Time dynamics

The [time dynamics module](core/time_dynamics.py) investigates model-specific clock phase, stability, recurrence, and Floquet calculations associated with the Time manuscript. Consult the [mathematical scope and reproduction guide](docs/TIME_DYNAMICS.md) for definitions and limitations.

Suggested local checks, subject to the dependencies and interfaces documented in the checkout:

```bash
python demo_time_dynamics.py
python -m pytest tests/test_time_dynamics.py
```

These are mathematical/software checks, not experimental validation of the manuscript's physical interpretation.

## Software organization

The directory includes a Python core, utilities, tests, demonstrations, and presentation material. For the authoritative current file inventory, inspect the checkout; older directory diagrams and feature lists are not guaranteed to reflect the implementation.

In particular, claims of complete multi-entity dynamics, universally stable identity layers, or a fully validated simulation engine require independent module-level verification and are **not made here**.

## Presentation rendering

Historical presentation-rendering tools may be available in the checkout, including a PowerShell-based renderer. Their outputs are visual communication artifacts and do not constitute scientific validation. Check scripts and dependencies before use.

## Research standards

- State mathematical assumptions before interpreting output.
- Keep exact proofs distinct from numerical experiments and visualizations.
- Record versions, seeds, parameters, tests, and failure cases.
- Do not infer physical behavior from a software demonstration.
- Preserve historical work while clearly identifying superseded claims.

For the broader program and verified mathematical checkpoints, see the [repository overview](../README.md).
