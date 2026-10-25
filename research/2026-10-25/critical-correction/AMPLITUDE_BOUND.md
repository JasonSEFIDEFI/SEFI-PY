# Conditional error bound for transverse-amplitude elimination

Working mathematical note, 25 October 2026. Exploratory / Not Promoted.

This is an a priori estimate for a prescribed charged-field profile. It is not a proof of the full coupled critical-limit expansion.

## Setting

Fix lambda>0 and delta>0. Let f be a prescribed sufficiently smooth real function and define

    s0^2=(3-4f^2)/lambda,
    c^2=1-4f^2/3.

Assume s0>0 and c>=c_*>0 throughout the domain. The exact static amplitude equation in z coordinates is

    -delta Delta s+g(s)=0,
    g(s)=lambda(s^2-s0^2)s.

Define

    chi_lambda=Delta s0/(2lambda s0^2),
    s_app=s0+delta chi_lambda.

Here chi_lambda uses the finite-lambda s0. Its leading delta->0 value is the chi used in the correction functional.

## Exact residual

Direct substitution gives

    R=-delta Delta s_app+g(s_app)
      =delta^2[-Delta chi_lambda+3lambda s0 chi_lambda^2]
        +delta^3 lambda chi_lambda^3.

No asymptotic truncation is used in this identity. check_amplitude_bound.py checks its algebra exactly.

## Conditional maximum-principle estimate

Suppose an exact real solution s and the approximation both obey

    s>=3s0/4,    s_app>=3s0/4.

Every intermediate value q between them then satisfies

    g'(q)=lambda(3q^2-s0^2)
          >=(11/16)lambda s0^2
          =(33/16)c^2 >= mu,
    mu=(33/16)c_*^2>0.

For e=s-s_app, the mean-value formula gives

    -delta Delta e+a(z)e=-R,
    a(z)=integral_0^1 g'(s_app+t e)dt >= mu.

On a bounded smooth domain with prescribed boundary error, the classical maximum principle yields

    ||e||_infinity
       <=max( ||e||_boundary,infinity, ||R||_infinity/mu ).

Indeed, a positive interior maximum exceeding the right side would have -Delta e>=0 and a e>||R||, contradicting the equation. Apply the same argument to -e for a negative minimum. On an unbounded domain, the same conclusion requires the stated decay and attainment/limiting conditions; it is not asserted for arbitrary bounded solutions.

With zero boundary error and sufficient integrability, multiplying the error equation by e and integrating gives

    delta ||grad e||_2^2+integral a e^2
       =-integral R e
       <=||R||_2 ||e||_2,

and hence ||e||_2<=||R||_2/mu.

## What this controls

If lambda, s0, chi_lambda and Delta chi_lambda have uniform bounds, c_* stays positive, the lower-bound assumptions hold, and the boundary error is O(delta^2), the displayed estimate gives O(delta^2) amplitude error for the prescribed f.

This supplies an explicit conditional control of one elimination step. Establishing the existence of the coupled branch, its uniform regularity, the positive-chart conditions and the inverse estimates for its angular equation remain separate tasks. No interval-certified verification of these hypotheses on the computed full branch is claimed here.

The estimate concerns the real stationary amplitude. It does not remove the full Sigma phase, establish a dynamical mass gap, or prove nonlinear stability.

