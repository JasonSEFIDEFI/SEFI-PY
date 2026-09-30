# N=0 full-Z falsification campaign — 30 September 2026

Starting commit: `bc68d639fe9b175f2e7c4fade05f991c0276157b` on
`research/action-audit-2026-09-30`. All existing tracked files, analytical
checkpoints, retained backgrounds, and v0.1 remain unchanged. This extension
imports the committed operator rather than changing it.

**Checkpoint result:** no accepted instability and no stability claim. The
completed L=16, h=0.5 campaign comprises 55 complex shift targets, 12 returned
pairs per target, and m=0–4 (660 returned pairs, including repeated modes).
The only positive-real values above the 1e-6 screening threshold belong to the
previous longitudinal translation splitting, approximately 0.00270686543,
with tangent overlap 0.98924 after removing carrier phase. Refinement and
independent targeting runs are in progress; this report will be updated.

## Audit and broader strategy

The committed B and G match the retained full-Z quadratic formulation,
including conjugate couplings, gyroscopic terms, time/axial advection, and
regular cylindrical axis spaces. No perturbation Z parity is imposed. N!=0
remains rejected. The original three operator tests pass.

A single shift-invert target returns eigenvalues near that target, not the
eigenvalues with largest real part. With a complex matrix the transformed
values are 1/(sigma-target); largest transformed magnitude therefore selects
nearest original eigenvalues. See the primary [SciPy eigs documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html).
The new adversarial test explicitly constructs an unstable oscillatory pair
missed by a near-zero search and recovered by broader targets and unshifted LR.

`broad_scan.py` uses 11 targets per sector: 0.0001 plus two positive real
offsets and five imaginary offsets. With extent E=max(0.25,1.1 R), targets
are (0.15 or 0.6) E + i(-1,-0.5,0,0.5,1) E. Additional explicit targets extend
the independent factorization search to imaginary offsets +/-1 and +/-2.
The Schur implementation factors B+target G-target^2 I, then reconstructs the
first-order inverse by exact block elimination. The independent `full`
strategy factors the original first-order matrix. Both are ARPACK methods;
they are independent factorizations, not independent eigensolver libraries.
An unshifted largest-real-part option is implemented, with nonconvergence
recorded rather than interpreted as absence of growth.

## Energy structure and search bound

Cylindrical discrete weights are proportional to i h^2 at r=i h; the m=0
axis weight is h^2/8. These satisfy the axis detailed-balance identity
w0*4=w1/2. In this inner product B is Hermitian and G is skew-Hermitian;
the measured defects at L=16, h=0.5 are at most 2.7e-15.

Taking the weighted inner product of sigma^2 y=B y+sigma G y, write
beta=<y,B y>/<y,y> real and <y,G y>/<y,y>=i g with g real. If sigma=a+i b
and a!=0, its imaginary part gives g=2b. Its real part then gives
|sigma|^2=beta<=lambda_max(B). Thus all nonimaginary eigenvalues of the exact
finite discretization lie inside this disk when its upper endpoint is known.

The code estimates lambda_max(B) with an extremal Hermitian solve and records
its residual. This is **not a certified upper bound**: a small residual does
not prove the Ritz value is the largest eigenvalue. It informs target placement
only. The L=16, h=0.5 m=0 estimate is 0.04196872147, giving R≈0.2048626893.
The estimate is negative for m=1 and higher in the runs completed so far;
this is not promoted to a stability conclusion. Neither the finite target
mesh nor ARPACK convergence certifies completeness, multiplicities, or absence
of unreturned right-half-plane eigenvalues.

## Acceptance gates and artifacts

Every returned pair records first-order, quadratic absolute/scaled, and
velocity-consistency residuals; cylindrical boundary norm fraction; peak and
RMS location; neighboring-node variation; and weighted symmetry overlaps.
The carrier phase is removed before measuring longitudinal translation overlap.
Representative target-00 eigenvectors and every nonsymmetry growing candidate
are retained; all eigenvalues and diagnostics are retained at every target.
`saved_indices` maps each saved vector column to its diagnostic record.

No candidate is automatically accepted. Acceptance requires small first-order
and quadratic residuals, independent h refinement, increasing L, a resolved
interior-localized eigenfunction, exclusion of symmetry/finite-box splitting,
and repeated targeting. Growth below the 1e-6 reporting threshold remains an
unresolved possibility; raw eigenvalues are retained. Carrier phase direct-null
checks and translation grid/domain calibration remain mandatory.

The broad-search tests additionally compare the Schur inverse to a dense block
solve and test weighted structure and the growth-disk identity against complete
dense spectra in m=0–4 on a tiny grid. All six tests pass. The tiny-grid oracle
tests numerical machinery; it does not represent the physical profile.

Reproduce from this directory with Python and `requirements-spectral.txt`:

```text
python -m unittest test_full_z_spectral test_broad_scan -v
python broad_scan.py ../../sim/vortex_console/data/charged_ring_refined.npz --output broad_runs/L16_h050 --n 32 --h .5 --m 0 1 2 3 4 --k 12
python broad_scan.py ../../sim/vortex_console/data/charged_ring_refined.npz --output broad_runs/L16_h050_independent --n 32 --h .5 --m 0 2 3 --k 12 --shifts .01 .12+.12j .12-.12j .1+1j .1-1j .2+2j .2-2j --strategies full
python run_broad_suite.py
python compact_broad_vectors.py broad_runs
python summarize_broad.py broad_runs
```

`broad_runs/` contains backgrounds, per-target machine-readable spectra,
selected eigenvectors, and run metadata. `convergence.json` summarizes scope,
residual maxima, candidates, symmetry trends, and artifact hashes. Original
full-Z calibration and v0.1 artifacts are not overwritten.
