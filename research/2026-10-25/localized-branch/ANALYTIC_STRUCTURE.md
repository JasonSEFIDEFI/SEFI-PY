# Analytic structure and conditional phase-sector result

Working derivations, 25 October 2026. These are model-specific statements; they do not establish novelty or the full stability of the computed branch.

## Exact potential separator

For rho,eta>=0 and lambda>0, with m_Phi^2=1+12/lambda:

    V - (m_Phi^2/2) rho
      = (lambda/4)(eta - 3/lambda + 4rho/lambda)^2
        + (1/4 - 4/lambda) rho^2.

The identity is checked exactly by SymPy in summarize.py. For lambda>=16 the right side is nonnegative. Combining it with kinetic energy and |Im(Phi*Phi_t)|<=|Phi||Phi_t| yields:

    E >= m_Phi |Q_Phi|.

This excludes strictly below-free-threshold Phi-charge binding over the stated vacuum, including Q_Sigma=0. It does not exclude all stationary, metastable or medium-supported structures. This separator was already present in the programme's public background analysis; it is not a new discovery here.

For lambda<16, minimization over rho>0,eta>=0 of 2V/rho gives the subthreshold value:

    omega_min^2 = 1 + 3/sqrt(lambda),
    rho = 3/sqrt(lambda), eta=0.

The vacuum mass threshold is 1+12/lambda. A subthreshold well permits but does not prove a localized solution. Gradients, boundary conditions and the charge constraint must still be solved.

## A conditional analytic check of the phase blocks

Let a smooth stationary solution on R^3 have strictly positive real amplitudes f and s, with adequate decay/regularity. Its equations imply:

    L_Phi f = 0,
    L_Sigma s = 0,

where

    L_Phi = -Delta + 1+f^2+4s^2-omega^2,
    L_Sigma = -Delta -3+4f^2+lambda s^2.

For a positive function u obeying (-Delta+W)u=0 and a smooth compactly supported perturbation z, integration by parts gives the standard ground-state representation:

    integral [|grad z|^2 + W z^2] d^3x
      = integral u^2 |grad(z/u)|^2 d^3x >= 0.

To see this, put z=u v, expand |grad(uv)|^2, and integrate the cross term together with v^2|grad u|^2 into -u Delta u v^2. The stationary equation cancels W u^2 v^2. Compact support removes the boundary term.

Apply it to u=f and u=s separately. Under these hypotheses, both imaginary-field/phase energy quadratic forms are nonnegative. For spherical harmonics, the centrifugal contribution l(l+1)/r^2 is nonnegative as well. Extending the statement beyond compactly supported perturbations requires the appropriate quadratic-form domain and control of boundary terms; it is not silently assumed.

This is a standard identity applied conditionally to this model, not a new stability theorem. Floating-point profiles do not prove strict continuum positivity or existence. The Sigma phase is gapless in the vacuum, so nonnegativity does not provide a positive spectral gap. Finite-box Sigma phase eigenvalues must not be called localized internal modes without domain/localization checks.

The coupled real-amplitude block, fixed-charge constraint, physical-time dynamics and finite-amplitude perturbations require separate analysis. In particular, phase-block nonnegativity does not establish global energetic optimality or nonlinear stability.

## Charge first law

On a differentiable stationary fixed-charge branch, differentiating E_Q[f(Q),s(Q)] gives:

    dE/dQ = partial E_Q/partial Q = Q/I = omega.

The profile-derivative terms vanish by stationarity and the boundary conditions. Thus a centered finite-difference derivative of independently computed branch energies supplies a consistency check. It is not an independent proof of stability.

## Fragmentation comparison

For a parent candidate with positive charge Q and two infinitely separated candidate daughters at Q1 and Q2=Q-Q1, define:

    Delta E = E_branch(Q1)+E_branch(Q2)-E_branch(Q).

Positive Delta E disfavors only that compared channel. Our computed branch energies are candidate energies, not certified global minima at each charge. Other branches, additional charges, nonspherical configurations, finite-separation effects, radiation and more complicated partitions are not exhausted.

Consequently all 26 positive sampled margins are encouraging branch evidence, not a theorem excluding fission.

