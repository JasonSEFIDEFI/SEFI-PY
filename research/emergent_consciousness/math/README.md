# Mathematical Research: Discrete-Time Stability and Robustness

**Research program:** SEFI / GWFM / DEFI
**Researcher:** Jason Duran Dutton
**Status:** Exact robust stability theorem established for a specified
two-dimensional switched linear system. Generalization remains open.

## 1. Research scope

This directory contains mathematical investigations of discrete-time
dynamical systems, including convergence, retention, switching,
perturbation robustness, and common quadratic Lyapunov certificates.

The calculations are mathematical model studies. They do not, by
themselves, establish physical properties of SEFI, GWFM, DEFI,
consciousness, or quantum systems.

## 2. Exact robust stability theorem

Consider

    x[k+1] = M[k] x[k]

where M[k] is chosen arbitrarily at every step from

    F(epsilon) = {B + E, C + E : |E[i,j]| <= epsilon}

with

    B = [[4/5, -2],
         [0,   3/10]]

    C = [[2/5, 6/5],
         [0,   7/10]]

and epsilon >= 0.

The family is uniformly exponentially stable under arbitrary
switching if and only if

    0 <= epsilon < 7/145.

The exact critical radius is

    epsilon_critical = 7/145
                     = 0.048275862068965517...

## 3. Exact common quadratic certificate

At the critical radius, define

    P = [[3/50,     -562/3625],
         [-562/3625, 47/50]]

Exact rational verification establishes:

1. P is symmetric positive definite.
2. P - M.T P M is positive semidefinite at all 32 vertices
   of the two entrywise uncertainty boxes.
3. Consequently, every matrix in the critical uncertainty
   family is nonexpansive in the induced P-norm.
4. Both unperturbed matrices B and C strictly contract
   in the same P-norm.

For 0 <= epsilon < epsilon_critical, write

    theta = epsilon / epsilon_critical

and

    A + epsilon D
      = (1-theta) A
        + theta (A + epsilon_critical D).

Convexity of the induced matrix norm, strict contraction of
the nominal matrices, and nonexpansiveness of the critical
uncertainty family imply a common contraction factor below one.

This proves uniform exponential stability under arbitrary switching.

## 4. Exact spectral obstruction

At epsilon = 7/145, the admissible matrix

    M_star = B + (7/145) [[ 1, -1],
                          [-1,  1]]

has eigenvalues

    1 and 57/290.

Its eigenvalue 1 prevents uniform asymptotic stability.

Because the uncertainty sets are nested as epsilon increases,
the same obstruction applies for every epsilon >= 7/145.

Together, the contraction certificate and spectral obstruction
establish the exact threshold.

## 5. Verification

Principal final verification scripts:

    python retry_boundary_metric_exact.py
    python verify_exact_threshold_theorem.py
    python verify_active_vertex_spectral_bound.py

Additional scripts document earlier numerical searches,
diagnostic investigations, and intermediate mathematical checks.

The exact theorem depends on rational symbolic verification,
not on floating-point solver accuracy.

Suggested Python dependencies include:

    sympy
    numpy
    scipy
    cvxpy

Not every script requires every dependency.

## 6. Research limitations

The exact threshold theorem applies to the specified matrices,
uncertainty structure, and arbitrary-switching dynamics.

It does not establish a general exact-threshold algorithm for
arbitrary uncertain switched systems.

The novelty of the result relative to published robust-control
literature has not yet been established.

No physical interpretation is inferred from this theorem.

## 7. Checkpoint 041: complete trace-one boundary-metric interval

For the boundary family

    P(r) = [[1-r, (608*r-594)/145],
            [(608*r-594)/145, r]],

the exact set satisfying P(r) positive definite and all 32 boundary
vertex inequalities P(r)-M.T*P(r)*M positive semidefinite is the
**closed interval** [r_minus, r_plus], where

    r_minus = (1585543261 - 145*sqrt(163411650609))/1654143502
    r_plus  = (4871334677 + 1595*sqrt(12061387921))/5143279214

Approximately: [0.923092923378325734, 0.981184241759181940].

The exact positive-definite range is (484/533, 729/733).
Sufficiency follows by checking both diagonal entries and the
determinant of each of 32 vertex margins: exact real-root isolation
and algebraic sign comparisons establish nonnegativity over the
closed interval. The matrix-norm convexity argument extends the
vertex inequalities to each full uncertainty box.

Necessity within the trace-one family follows from two B-box
vertex determinants, proportional respectively to the negatives of

    827071751*r**2 - 1585543261*r + 758855474
    2571639607*r**2 - 4871334677*r + 2303901468.

Their sign constraints, together with P(r) positive definite,
restrict r to precisely the stated interval. At either endpoint
a margin becomes singular; P(r) remains positive definite.

Reproduce with:

    python verify_checkpoint_041_exact_interval.py

The verifier was independently executed using SymPy and passed
96 exact scalar checks. Its decimal output is for presentation only;
the proof uses exact rational/algebraic arithmetic.

**Scope:** This characterizes the specified trace-normalized
one-parameter boundary family, not arbitrary Lyapunov matrices.
It does not establish a physical SEFI, GWFM, DEFI, or consciousness
interpretation. Literature novelty remains unassessed.

## 7a. Next research stop-gate

Before generalizing, determine whether the unit-eigenvector
constraint exhausts all common quadratic boundary metrics up to
positive scaling, and distinguish this from general (nonquadratic)
stability certificates. Do not infer physical conclusions.

## 8. Publication status

This is a reproducible mathematical research result.

Journal-level novelty, generality, and suitability remain
subject to literature review and further investigation.
