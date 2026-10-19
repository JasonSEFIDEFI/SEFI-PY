# Verification - 19 October 2026

- Exact dispersion and numerical eigenvalues agree in 12 benchmark mode cases; RK4 step refinement reduces error.
- `python -m unittest test_dispersion test_console test_rest_console test_warp_vortex -v`: 20 tests passed.
- After UTF-8 integration cleanup, the five new dispersion/API tests passed again.
- `node --check static/dispersion.js`: passed.
- `git diff --check`: passed.
- Browser automation failed to initialize twice. The page was not visually verified; live HTTP page/script/API delivery was tested.
- This is a far-field Fourier-mode comparison, not a nonlinear ring collision or stability simulation.
