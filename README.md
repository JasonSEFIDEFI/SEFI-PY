# SEFI-PY — SEFI Research Program

**Geometric foundations · Mathematical analysis · Reproducible computational research**

**Researcher:** Jason Duran Dutton  
**Status:** Active mathematical and computational research. Proposed physical interpretations remain unvalidated.

## Central research question

> Can a minimally constrained mathematical structure produce geometry, distinguishable physical entities, and their interactions as emergent properties, rather than requiring those structures to be specified beforehand?

The SEFI Research Program examines this question through mathematical derivation, reproducible computation, and critical comparison with established physics. The purpose is to discover what follows from the assumptions—not to presume a particular physical conclusion.

## Frameworks under investigation

| Framework | Research direction |
| --- | --- |
| **GWFM — Geometric Worldline Foundations Model** | Worldline-based geometric foundations |
| **SEFI — Single Entity Field Interpretation** | Proposed mathematical interpretation of physical entities |
| **DEFI — Dynamic Entity Field Integration** | Coupled dynamics and interacting structures |
| **Warp geometry** | Exploratory geometric transformations and constraints |

These are developing research proposals, **not established physical theories**.

## Research standards

- State assumptions and dependencies explicitly; use exact derivations wherever feasible.
- Distinguish mathematical proofs, model-dependent calculations, numerical experiments, physical measurements, and hypotheses.
- Stop at unresolved mathematical dependencies rather than inventing missing structure.
- Retain unsuccessful experiments, counterexamples, and historical records.
- Compare results with established literature before claiming novelty.
- Do not treat passing software tests as evidence of experimental physical validation.

## Verified mathematical result: robust stability

A completed research checkpoint establishes an **exact stability threshold for one specified uncertain switched linear system**:

```text
x[k+1] = M[k] x[k]

M[k] in {B + E, C + E : |E[i,j]| <= epsilon}

B = [[4/5, -2],      C = [[2/5, 6/5],
     [0,   3/10]]         [0,   7/10]]
```

For arbitrary switching, the family is uniformly exponentially stable **if and only if**

```text
0 <= epsilon < 7/145
```

The proof employs an exact rational common quadratic Lyapunov certificate, verifies the 32 vertices of the critical uncertainty boxes, and identifies a boundary matrix with eigenvalue 1. The threshold is a mathematical result for **this model**, not a validation of SEFI or a general robust-control algorithm. Novelty relative to existing literature has not been established.

- [Exact theorem and certificate](https://github.com/JasonSEFIDEFI/SEFI-PY/blob/research/emergent-consciousness-math/research/emergent_consciousness/math/README.md)
- [Verification scripts](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/research/emergent-consciousness-math/research/emergent_consciousness/math)

**Open checkpoint:** Exact characterization of the full admissible boundary-metric family.

## Repository navigation

| Area | Scope |
| --- | --- |
| [Current mathematical investigations](https://github.com/JasonSEFIDEFI/SEFI-PY/tree/research/emergent-consciousness-math/research/emergent_consciousness/math) | Stability proofs, symbolic verification, and research checkpoints (research branch) |
| [October 18 research](research/2026-10-18/) | Dated dynamics, coupling, and visualization investigations |
| [October 19 research](research/2026-10-19/) | Dated vortex investigations |
| [SEFI-PY Evolution](SEFI-PY%20Evolution/) | Earlier simulations, educational demonstrations, and engineering prototypes |
| [SEFI-PY implementation](SEFI-PY/) | Software and benchmarks |
| [QEC Stack](SEFI%20QEC%20Stack/) | Quantum information and error-correction software investigations |

These folders represent different stages of development. Their inclusion does not mean that all earlier models, claims, or experimental interpretations remain endorsed. Project-level documentation is undergoing review.

## Reproducibility and limitations

Use each project's own instructions, inputs, and verification criteria. For the exact stability checkpoint, activate a Python environment with SymPy and run from the research branch:

```bash
python research/emergent_consciousness/math/verify_exact_threshold_theorem.py
```

A verified property of a selected mathematical model does not establish that the model describes nature. Numerical simulations, software unit tests, and physical experiments have different evidentiary roles.

## Manuscripts and academic work

The wider program includes manuscripts submitted for scholarly evaluation. **Submission, peer review, acceptance, and publication are different stages.** This README makes no claim of journal acceptance or experimental confirmation. The suitability and novelty of the exact stability result remain subjects for literature review.

Jason Duran Dutton is pursuing Information Technology studies at American InterContinental University, with the stated objective of continuing through an MSIT. An intended degree is not represented as awarded. His professional background includes industrial electrical maintenance, controls, troubleshooting, and safety-regulated technical work.

## Related repositories

- [SEFI Manuscripts](https://github.com/JasonSEFIDEFI/SEFI-Manuscripts)
- [SEFI QEC Stack](https://github.com/JasonSEFIDEFI/SEFI_QEC_Stack)
- [Doctoral research portfolio](https://github.com/JasonSEFIDEFI/PhD)
- [Professional portfolio](https://github.com/JasonSEFIDEFI/CV-)

---

**Research principle:** Mathematical results determine which conclusions are justified. Where assumptions, derivations, or validation remain incomplete, the research remains open.
