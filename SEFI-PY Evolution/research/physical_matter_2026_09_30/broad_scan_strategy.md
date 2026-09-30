# Broad-search method and limits

The committed full-Z B and G agree with the retained quadratic formulation,
including conjugate couplings, gyroscopic frequencies, time/axial advection,
and the regular angular axis spaces. No perturbation Z parity is imposed.
N!=0 remains rejected. The original operator and analytical checkpoints are
imported unchanged; the new code supplies targeting and diagnostics.
The old report references external full-Z run archives, which are not committed
in this checkout. The audit used committed code, tests, calibration JSON, and
fresh reproduced calculations, rather than assuming those absent archives
had been independently inspected. Reproduced neutral-frequency differences
are recorded in convergence.json; carrier phase uses its direct null identity.

## Why the previous search can miss growth

For a complex matrix, shift-invert maps an eigenvalue sigma to
1/(sigma-target). Selecting largest transformed magnitude selects eigenvalues
nearest the target, rather than those with largest original real part. The
primary [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html)
specifies this transformation and the unshifted LR option. The adversarial
test constructs an unstable oscillatory pair masked by neutral near-zero
oscillations, then recovers it with multiple targets and LR.

The default search uses 0.0001 and ten complex targets:
(0.15 or 0.6) E + i(-1,-0.5,0,0.5,1) E, E=max(0.25,1.1 R).
R is an estimated growth-disk radius below. Additional independent targets
reach imaginary offsets +/-1 and +2. A finite k-nearest target mesh is not
complete; closely packed box modes, Arnoldi nonconvergence, and unreturned
eigenvalues remain possible. Returned spectral ranges are not filled regions
of certified exclusion.

The Schur targeting strategy factors Q=B+target G-target^2 I. For
(A-target I)(y,v)=(f,z), v=f+target y and
Q y=z-G f+target f. The full strategy instead factors A-target I directly.
The dense inverse test verifies the block identity. Both strategies use
ARPACK, so agreement checks two factorizations, not two eigensolver libraries.
The implemented unshifted LR attempt is retained even when it does not
converge. Its failure does not exclude growth.

## Weighted energy identity

Radial discrete weights are i h^2 at r=i h, with m=0 axis weight h^2/8.
The axis satisfies w0*4=w1/2. In this inner product B is Hermitian and G
skew-Hermitian; measured transformed defects are stored for every sector.
The conjugate blocks and cylindrical weights are essential to this check.

Take the weighted inner product of sigma^2 y=B y+sigma G y. Define real
beta=<y,B y>/<y,y> and real g by <y,G y>/<y,y>=i g. For sigma=a+i b with
a!=0, the imaginary equation gives g=2b; the real equation then gives
|sigma|^2=beta<=lambda_max(B). Thus the exact finite discretization's
nonimaginary eigenvalues lie inside this disk if its endpoint is known.

The largest Hermitian eigenvalue is estimated with eigsh; its residual is
recorded. A small residual does not certify the largest eigenvalue, so the
estimate is explicitly not a certified upper bound. It guides targeting only.
Negative Ritz estimates in nonzero m sectors are not promoted to stability
claims. Complete dense tiny-grid spectra test the identity independently of
the physical-profile runs. The early campaign used ARPACK's default bound
starting vector; the final implementation fixes its seed. Reproduction is
numerical within tolerances, not promised byte-for-byte across machines.

## Acceptance and saved diagnostics

Every returned pair stores first-order and absolute/scaled quadratic residuals,
velocity consistency, weighted symmetry overlaps, boundary norm fraction,
peak/RMS location, and neighboring-node variation. The boundary quantity is
sqrt(outer-strip weighted energy / total weighted energy), with the outer
10% radial or axial strip. It must not be confused with an energy fraction.
Longitudinal overlap removes the carrier phase from both vectors.

All eigenvalues and diagnostics are retained. Target-00 eigenvectors and every
nonsymmetry positive candidate are retained; saved_indices maps vector columns
to records. The 1e-6 real-part threshold is a reporting threshold, not an
instability acceptance criterion or exclusion of slower growth. Near-zero
scaled quadratic residuals can be sensitive to a vanishing normalization;
absolute residuals and direct carrier-null checks are retained alongside them.

Acceptance requires small first-order and quadratic residuals, independent h
refinement, increasing L, a resolved interior-localized eigenfunction,
exclusion of symmetry/finite-box splitting, and repeated targeting where
practical. No automated acceptance rule is supplied. Exact infinite-domain
translation neutrality remains unresolved; its finite-box splittings must
continue to shrink before any claimed physical growth is accepted.

## Reproduction

Use the pinned original requirements-spectral.txt and Python. Run:

```text
python -m unittest test_full_z_spectral test_broad_scan -v
python broad_scan.py ../../sim/vortex_console/data/charged_ring_refined.npz --output broad_runs/L16_h050 --n 32 --h .5 --m 0 1 2 3 4 --k 12
python run_broad_suite.py
python broad_scan.py ../../sim/vortex_console/data/charged_ring_refined.npz --output broad_runs/L16_h050_independent --n 32 --h .5 --m 0 --k 12 --shifts .01 .12+.12j .12-.12j .1+1j .1-1j .2+2j --strategies full
python broad_scan.py ../../sim/vortex_console/data/charged_ring_refined.npz --output broad_runs/L16_h050_independent_m23 --n 32 --h .5 --m 2 3 --k 12 --shifts .01 .12+.12j .12-.12j --strategies full --maxiter 200
python broad_scan.py ../../sim/vortex_console/data/charged_ring_refined.npz --output broad_runs/L16_h050_LR --n 32 --h .5 --m 0 --k 8 --strategies LR --maxiter 200
python finalize_broad.py broad_runs
python write_broad_report.py broad_runs
```

Output directories are restartable and reject different source/grid/target/k
settings. Existing incomplete records are retained, not silently replaced;
use a new output directory to retry failed calculations. Per-target JSON
records the actual settings; settings missing from early records are explicitly
reconstructed from the campaign execution log. The -2 imaginary-offset target
was interrupted twice without returned eigenpairs and is not in the successful
scope. A resumed L24 calculation retained completed targets. Artifact hashes
and the additive-only baseline preservation check are in verification.json
and convergence.json.
Text hashes normalize CRLF to LF to survive Git line-ending conversion;
binary profile/eigenpair hashes use raw bytes. The scheme is explicit in JSON.
