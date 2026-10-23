# Research notebook: critical binding in the repository's two-scalar model

This is mathematical exploration, not a manuscript or a claim of new physical matter. The model comes from SEFI-PY/research/2026-10-19/vortex-extended/vacuum_spherical.py. Its relationship to nuclear interactions or the SEFI foundational action has not been derived. All quantities below are dimensionless; lambda > 0.

## Exact calculation

For real amplitudes f,s and fixed global charge Q, let I=integral f^2 and

E = integral [ (|grad f|^2+|grad s|^2)/2 + U ] + Q^2/(2I),

U = f^4/4 + f^2/2 + 2 f^2 s^2 + lambda (s^2-3/lambda)^2/4.

The vacuum is (0,sqrt(3/lambda)); the free charged-field mass obeys m^2=1+12/lambda. Direct completion of the square gives

U - m^2 f^2/2 =
lambda/4 [s^2+(4 f^2-3)/lambda]^2
+ (lambda-16)f^4/(4 lambda).

SymPy verifies this identity exactly in critical.py.
For lambda >= 16 this proves E >= m |Q| by positivity of gradients and the arithmetic-geometric mean inequality. Thus no state with energy below separated free quanta exists in this parameter range. This statement does not exclude above-threshold stationary states.

Minimizing 2U/f^2 over s^2>=0 gives, with r=f^2:

- r<=3/4: m^2 + (lambda-16)r/(2 lambda).
- r>=3/4: 1+r/2+9/(2 lambda r).

For lambda<16 the minimum is omega_min^2=1+3/sqrt(lambda), at r=3/sqrt(lambda),s=0. Large-radius constant-core trial states have E/Q approaching omega_min<m: volume energy grows as radius cubed and interface energy as radius squared. This is an energetic trial-state argument, not a compactness/existence proof or proof of dynamical stability.

## Formal near-critical reduction

At lambda=16, U-m_c^2 f^2/2=(f^2/2+2s^2-3/8)^2, m_c=sqrt(7/4).
The zero valley can be parameterized by f=sqrt(3)/2 sin(theta), s=sqrt(3)/4 cos(theta).

Set delta=16-lambda>0, omega^2=m^2-delta alpha, rho=sqrt(delta) R.
The induced gradient metric and leading effective potential are

G(theta)=3(1+3cos^2(theta))/16,
V_alpha(theta)=3 alpha sin^2(theta)/8 - 9 sin^4(theta)/1024.

Both expressions are checked symbolically, including the lambda-dependent vacuum mass subtraction. The reduced action is

S_red = 4 pi integral rho^2 [G theta'^2/2 + V_alpha] d rho,

and its radial equation is G(theta''+2theta'/rho)+G' theta'^2/2=V_alpha'.
Boundary conditions: theta'(0)=0, theta(infinity)=0.
The localized frequency window is 0<alpha<3/128.

The formal scaling predictions are

Q ~ delta^(-3/2) q(alpha),
q=m_c 4 pi integral rho^2 (3/4)sin^2(theta) d rho,

(E-mQ) ~ delta^(-1/2) [S_red-alpha q/(2m_c)].

These are asymptotic predictions conditional on convergence to the valley and the localized solution branch. They are not yet rigorous singular-limit theorems. In particular the scaling of a minimum charge requires branch/fold control, beyond one fixed-center profile.

## Numerical protocol and limitations

critical.py fixes the central angle and solves for alpha with SciPy solve_bvp, including its singular-origin term. It records all solves and rejects nonpositive alpha, nodes/overshoot, failed residual convergence, or central angles beyond the chosen valley segment. The initial central-angle 0.3 run converged numerically but has negative alpha and is excluded.

validate.py refines domain and tolerance and solves the original two-field equations at several delta values, independently of the reduced field equation. Initial direct guesses for the refinement run failed (overflow/singular Newton Jacobian); validation-failed-initialization.log preserves those failures. The corrected run uses continuation from accepted small-center solutions.

Residual tolerances are numerical diagnostics, not statistical confidence intervals. Stability needs the fixed-charge second variation and time evolution; the sign of E-mQ tests only energetic dispersion into free quanta.

## Prior art and novelty status

Two-field charge-stabilized solitons are established physics:
- Heeck and Sokhashvili, Revisiting the Friedberg-Lee-Sirlin soliton model (2023), https://link.springer.com/article/10.1140/epjc/s10052-023-11710-9
- Hamada et al., Q-Balls in the presence of attractive force (2024), https://arxiv.org/abs/2407.11115

The repository includes an additional charged-field quartic term compared with the simplest FLS model. No novelty is claimed for charge stabilization, the frequency criterion, thin/thick-wall branches, or minimum charges. The exact parameterized critical reduction and its quantitative corrections are candidates for further comparison with prior work. A preliminary search does not establish novelty.

## Next discriminating work

1. Verify full-model profile/charge convergence to the reduced curve.
2. Locate the fold and energetic crossing with domain and tolerance checks.
3. Compute fixed-charge Hessian eigenvalues, including nonspherical perturbations.
4. Derive the first correction in delta and test its predicted coefficient on held-out delta values.
5. Check prior literature specifically for the quartic-extended model and critical flat-valley limit.
6. Only if warranted, derive a connection to measurable matter. No mapping to elements or planets follows from this calculation.
