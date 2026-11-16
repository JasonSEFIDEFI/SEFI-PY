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

## 7. Next research checkpoint

Checkpoint 041:

Characterize symbolically the complete set of positive-definite
boundary metrics satisfying all 32 semidefinite inequalities.

The boundary metric family obtained from the unit-eigenvector
constraint and trace normalization is

    P(r) = [[1-r, (608*r-594)/145],
            [(608*r-594)/145, r]].

Determine its admissible parameter set using exact algebraic
inequalities rather than a numerical grid.

Do not advance to broader claims until the required
mathematical dependencies are established.

## 8. Publication status

This is a reproducible mathematical research result.

Journal-level novelty, generality, and suitability remain
subject to literature review and further investigation.
