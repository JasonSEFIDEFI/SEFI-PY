# Stationary obstruction and a formal critical scaling problem

Working derivations, 25 October 2026. Exploratory / Not Promoted. These arguments require independent technical review and a precise literature comparison before submission.

## 1. Stated model and frequency window

With canonical real-component kinetic normalization, rho=|Phi|^2, eta=|Sigma|^2 and lambda>0:

    V = rho^2/4 + rho/2 + 2 rho eta + lambda(eta-3/lambda)^2/4.
    m^2 = 1 + 12/lambda.
    omega_min^2 = 1 + 3/sqrt(lambda),   lambda<16.
    delta = m^2 - omega_min^2 = 3(4-sqrt(lambda))/lambda.

For the ansatz Phi=exp(i omega t) f(x), Sigma=s(x), the stationary reduced functional is

    S_omega = integral [(|grad f|^2+|grad s|^2)/2 + W_omega(f,s)] d^3x,
    W_omega = V - omega^2 f^2/2.

The boundary is f->0, s->sqrt(3/lambda), with fixed asymptotic phase. Coordinates and background spacetime are assumed.

## 2. Nonexistence within this stationary class for lambda>=16

The exact decomposition is

    W_omega
      = (lambda/4)(s^2-3/lambda+4f^2/lambda)^2
        + (1/4-4/lambda) f^4
        + (m^2-omega^2) f^2/2.

If lambda>=16 and omega^2<=m^2, every term is nonnegative.

Assume f,s are regular stationary critical points with finite gradient energy K and integrable W_omega, and that dilation/integration by parts produces no surface contribution. Under f_a(x)=f(x/a), s_a(x)=s(x/a),

    S_omega(a) = a K + a^3 U_omega,
    0 = dS_omega/da at a=1 = K + 3 U_omega.

K>=0 and U_omega>=0 imply K=0; the prescribed boundary then forces the trivial vacuum. Thus there is no nontrivial localized finite-energy stationary solution in this class under those assumptions. Spherical symmetry is not needed for this dilation argument.

This is a standard Derrick/Pohozaev application to the displayed potential. It does not exclude arbitrary time-dependent states, charged media, different boundary conditions, additional time-dependent Sigma charge or a different action. No conclusion is asserted here about embedded states outside the stated frequency range.

## 3. Deriving a candidate limit as lambda approaches 16 from below

Set

    omega^2 = m^2 - delta beta,  0<beta<1,
    z = sqrt(delta) r.

The squared term penalizes deviations from the curve

    s^2 = (3-4f^2)/lambda.

Assuming a smooth limiting localized branch remains close to that curve with 0<=f^2<=3/4, restrict W_omega to it. Exactly on the curve:

    W_omega/delta = beta rho/2 - (4+sqrt(lambda)) rho^2/12.

Consequently its limiting potential is beta rho/2 - 2 rho^2/3. Parameterize the limiting curve by

    f = (sqrt(3)/2) sin(theta),
    s = (sqrt(3)/4) cos(theta),
    0<=theta<=pi/2.

The induced kinetic coefficient is

    a(theta) = (df/dtheta)^2+(ds/dtheta)^2
             = (3/16)(1+3 cos^2(theta)).

The candidate reduced radial functional is

    A_beta = 4 pi integral z^2 [
                   a(theta) theta_z^2/2
                   + (3/8) sin^2(theta)(beta-sin^2(theta)) ] dz.

Its Euler equation is

    a(theta)(theta_zz+2 theta_z/z)
       + a'(theta) theta_z^2/2
       - (3/4) sin(theta)cos(theta)(beta-2 sin^2(theta)) = 0.

Boundary conditions are theta_z(0)=0 and theta(infinity)=0. The code uses x=sqrt(beta) z with a finite tail regulator X.

This is a formal constrained-curve reduction. Uniform control of the massive transverse field, possible layers, endpoint behavior, existence, and a quantified convergence theorem are not supplied. Numerical agreement cannot replace those missing estimates.

The potential identities and induced metric are checked symbolically in summarize_frequency.py.

## 4. Charge and binding predictions

Since Q=omega integral f^2 d^3r and m approaches m_c=sqrt(7)/2:

    Q delta^(3/2) -> Qbar(beta)
       = m_c 4 pi integral z^2 (3/4) sin^2(theta) dz.

Let K_z denote the reduced gradient contribution and R4=4 pi integral z^2 rho^2 dz. Expanding omega-m = -delta beta/(2m)+O(delta^2) gives

    (E/Q-m)/delta -> [K_z-(2/3)R4]/Qbar.

The beta rho/2 contribution in A_beta cancels the leading (omega-m)Q term. Therefore the limiting binding crossing is found by K_z=(2/3)R4, separately from minimizing Qbar(beta).

Computed reduced predictions at X=40,tol=1e-9:

    beta_turn ~= 0.12437038,     Qbar_turn ~= 38.53903.
    beta_binding ~= 0.26513103,  Qbar_binding ~= 50.59922.

The full-model values approach these numbers over the tested coupling range; coupling_table.csv reports discrepancies. These coefficients are numerical predictions of the formal reduced problem, not exact constants or experimentally measured quantities.

At those predicted beta values, enlarging the reduced domain from X=40 to 48 and tightening tolerance from 1e-9 to 1e-10 changes the scaled charges by less than 1.3e-10. This is consistency of numerical quadrature/collocation at selected parameters; it is not a certified bound on either the minimizer or the singular-limit approximation.

Because delta ~ 3(16-lambda)/128, the same formal calculation predicts growth proportional to (16-lambda)^(-3/2). This explains why approaching the no-solution separator can require increasingly large charge and spatial extent. It is not a proof that every branch obeys this scaling.

## 5. Equal-charge energy character

The two Q=1000 profiles at lambda=12 are checked with the original finite-volume fixed-Q Hessian, including its positive rank-one term. Newton polishing preserves access to the saddle; a minimization-only algorithm could discard it.

The high-frequency profile has a resolved negative radial amplitude curvature. The low-frequency comparison has positive computed radial amplitude curvatures. The slope of Q(omega) and the Hessian comparison are consistent in this sample, but a general spectral-stability theorem requires its own hypotheses and analysis.

Energetic saddle, local energetic minimum, global binding and nonlinear physical-time stability are distinct statements.

## 6. Literature context

Multi-field charge thresholds and thick-wall issues are established subjects: [Lennon (2021)](https://arxiv.org/abs/2112.14263). Coupled-field branch and stability analyses also appear in [Hamada et al. (2024)](https://arxiv.org/abs/2407.11115). Those references motivate a detailed comparison; they neither validate this computation nor establish its novelty.

A potential contribution is the specific separator, formal limiting variational problem and checked branch classification for this potential family. Whether this combination is original remains unresolved.

