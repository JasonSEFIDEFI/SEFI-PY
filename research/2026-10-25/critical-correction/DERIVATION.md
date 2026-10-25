# A formal first correction with conditional charge formulas

Working derivation, 25 October 2026. Exploratory / Not Promoted.

## 1. Definitions and chart

Use the same canonical model and real stationary profiles as the prior record. Let

    m_lambda^2=1+12/lambda,
    delta=12/lambda-3/sqrt(lambda),
    omega^2=m_lambda^2-delta beta,
    z=sqrt(delta) r,    m_c=sqrt(7)/2.

Write rho=f^2 and

    f=(sqrt(3)/2) sin(theta),
    s0_lambda=sqrt(3/lambda) cos(theta),
    s=s0_lambda+delta chi.

This is a fixed-f transverse-amplitude chart, not an orthogonal normal-coordinate claim. It is regular on the selected profiles where cos(theta) stays positive. The full Sigma phase is not being eliminated; this is a stationary real-amplitude calculation.

The scaled stationary action is A_delta=delta^(1/2) S_omega. With
v=1/sqrt(lambda)=(1+sqrt(1+16 delta/3))/8,

    v=1/4+delta/3+O(delta^2),
    m_lambda^2=7/4+2 delta+O(delta^2),
    omega=m_c+delta(2-beta)/(2m_c)+O(delta^2).

The potential decomposition is

    W_omega=(lambda/4)(s^2-s0_lambda^2)^2
             +(1/4-4/lambda) f^4 + delta beta f^2/2.

## 2. Leading problem

The leading angular metric and potential are

    a(theta)=(3/16)(1+3 cos^2(theta)),
    U_beta=beta rho/2-(2/3)rho^2.

Thus

    A0[theta]=integral [a(theta)|grad_z theta|^2/2+U_beta] d^3z.

The leading solution theta0(beta) satisfies the radial Euler equation, regularity at zero and decay at infinity. In this record it is computed by collocation with a scaled outer regulator, not proved to exist by interval methods.

## 3. First-order local functional

Before eliminating chi, direct Taylor expansion gives

    A_delta = A0
       + delta integral [
            (1/4)sin^2(theta)|grad theta|^2
            +(4/9)rho^2
            +grad s00 . grad chi
            +3 cos^2(theta) chi^2 ] d^3z
       + O(delta^2),

where s00=(sqrt(3)/4)cos(theta). The first two terms account for the changing vacuum-curve metric and restricted potential. Omitting them would give an incorrect charge correction.

After integration by parts, the chi terms are

    -chi Delta s00+3 cos^2(theta)chi^2.

Their stationary value occurs at

    chi0=Delta s00/[6 cos^2(theta)].

Completing the square yields the formal first correction

    A1[theta]=integral [
        (1/4)sin^2(theta)|grad theta|^2
        +(4/9)rho^2
        -(Delta s00)^2/[12 cos^2(theta)] ] d^3z.

The Laplacian on a leading solution can be evaluated without differentiating noisy data:

    Delta theta = U_beta,theta/a - a_theta |grad theta|^2/(2a),
    Delta s00=-(sqrt(3)/4)[sin(theta) Delta theta+cos(theta)|grad theta|^2].

The symbolic code checks the pre-elimination coefficient and completion-of-square identity exactly. The numerical code also inserts s0_lambda+delta chi0 into the original finite-delta functional.

The negative higher-gradient term is used perturbatively on theta0. We do not solve its fourth-order Euler equation as a fundamental theory, infer an unbounded physical energy of the full two-field model from it, or claim a new stable ultraviolet completion.

## 4. Envelope relation and charge correction

At fixed delta, stationarity gives the exact on-shell relation

    Q delta^(3/2)=2 omega partial_beta A_delta,on-shell.

Assume a smooth branch with a sufficiently differentiable uniform expansion on a beta neighborhood:

    A_delta,on-shell(beta)=A0(beta)+delta A1(beta)+O(delta^2).

Here A1(beta) means the functional above evaluated on theta0(beta). The unknown first profile correction does not enter this first on-shell action coefficient because the first variation of A0 vanishes.

Set I0=2 A0'(beta) and Q0=m_c I0. Then

    Q delta^(3/2)=Q0(beta)+delta Q1(beta)+O(delta^2),
    Q1=2m_c A1'(beta)+(2-beta)I0/(2m_c).

We evaluate A1' in two ways: a five-point centered action derivative at four step sizes, and a linear sensitivity problem.

For the sensitivity, use x=sqrt(beta) z and define

    b(theta)=a_theta/(2a),
    g(theta,beta)=[(3/4)sin(theta)cos(theta)/a]
                  [1-2 sin^2(theta)/beta].

The leading equation is theta_xx+2 theta_x/x+b theta_x^2-g=0. For u=partial_beta theta at fixed x,

    u_xx+2u_x/x+2b theta_x u_x
      +(b_theta theta_x^2-g_theta)u=g_beta,

with u_x(0)=0, u(infinity)=0. Differentiating the x-coordinate quadrature includes the beta^(-3/2) measure and every explicit beta term. The independent sensitivity result agrees with the frozen turn coefficient within about 6.2e-11.

## 5. Conditional formulas at a minimum and at binding

These formulas require a smooth local branch and nondegenerate roots; they do not establish those hypotheses for every continuum solution.

Suppose beta_t is a nondegenerate minimum of Q0(beta). The shift in its location contributes only at second order to its charge value. Therefore

    C0_t=Q0(beta_t),
    C1_t=Q1(beta_t).

For binding, the exact condition E/Q=m_lambda is equivalent to

    A_delta,on-shell=beta [Q delta^(3/2)]/(m_lambda+omega).

At its leading root beta_b, A0=beta_b Q0/(2m_c). Suppose this root is simple. Expanding the binding condition and including its shifted beta gives

    C0_b=Q0(beta_b),
    C1_b=(2m_c/beta_b)A1(beta_b)
           +Q0(beta_b)(4-beta_b)/(4m_c^2).

The Q1 terms cancel between the direct charge correction and the root displacement. Using Q1 at fixed beta_b alone would therefore give the wrong binding charge correction.

Computed coefficients:

    C0_t=38.5390285307,  C1_t=66.2516913306,
    C0_b=50.5992194999,  C1_b=76.6692699530.

The minimum is at beta_t approximately 0.12437038 and the binding root at beta_b approximately 0.26513103. The chart's minimum cos(theta) is about 0.449 at the turn and 0.150 at binding, so this calculation does not approach the chart's singular endpoint.

## 6. What the checks establish

[AMPLITUDE_BOUND.md](AMPLITUDE_BOUND.md) gives an explicit conditional estimate for the exact transverse equation at prescribed f. The finite-lambda approximation has an exact O(delta^2) residual formula; a positive lower bound on the mean reaction derivative yields a maximum-principle error bound. This controls one amplitude step under its hypotheses, not existence or uniform control of the coupled angular branch.

As a secondary calculation after the charge tests, expanding the extrema/root locations gives beta_t,1=-Q1'(beta_t)/Q0''(beta_t) and beta_b,1=(C1_b-Q1(beta_b))/Q0'(beta_b). Reduced-action differentiation and sensitivity give approximately 0.01183051 and 0.02868270. The sampled reduced curvature Q0'' at the turn is about 1751.74 and Q0' at binding about 159.50. These support, but do not rigorously certify, the numerical nondegeneracy assumptions. The position checks were not frozen before the full-model charge tests and are labeled secondary rather than prospective validation.

The held-out full-model discrepancies after correction are consistent with O(delta^2) over the two declared couplings. The correction improves charge predictions by factors of roughly 1,200 to 2,400. Larger-domain and tighter-solver checks are much smaller than those prediction errors.

These checks support the algebra and the selected asymptotic regime. They do not prove a uniform singular-perturbation theorem, a continuum existence theorem, uniqueness, global minimization, or nonlinear dynamical stability. A rigorous result would need function spaces, invertibility/coercivity estimates, decay control and a uniform remainder estimate.

