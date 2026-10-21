# SEFI-PY adaptive recovery benchmark

This self-contained research experiment evaluates conventional adaptive estimators for pilot-assisted bosonic displacement recovery. It does not implement a protecting QEC code, prove a SEFI-specific advantage, or connect to quantum hardware.

## Reproduce

Use Python 3.12+ and install `requirements.txt`. From this folder:

```text
python -m unittest -v test_adaptive
python benchmark.py development
python benchmark.py evaluation
python benchmark.py timing
python validate_saved_results.py
```

Development selects baseline parameters from its own seeds. Evaluation uses the saved selection and separate seeds. Each stage records completion status; failed/interrupted stages are not valid final evidence. Re-running overwrites the corresponding results; preserve the original results directory to compare reproductions. The protocol is an internal freeze, not a registered study.

To rebuild the PDF and vector figures, additionally install `report-requirements.txt` and run `python build_report.py`. Rendering report pages for visual review requires a PDF renderer such as Poppler; it is not needed for the numerical experiment. The scientific tests use the bundled physical reference and require no prior checkout.

Start with `ALGORITHM_AND_PROTOCOL.md`, `protocol.json`, and the generated `RESULTS.md`. `CLAIMS.md` separates supported statements from unsupported ones. `SUBMISSION_ROUTE.md` describes the current release and assessment route.

`physical_reference.py` is copied unchanged from Jason D. Dutton's previous displacement-recovery benchmark for independent overlap validation. Its analytical oscillator convention is [a,a†]=1. Raw truth is saved for evaluation and never passed to the estimators.

The adaptive mixture uses standard empirical Bayes and predictive expert weighting. The SEFI research motivation is not presented as a mathematical derivation of a new estimator. Original research-console equations and earlier benchmark files remain unchanged.

License: this addition follows the existing SEFI-PY repository license. See that repository's root LICENSE. No separate license grant is inferred for unrelated materials.
