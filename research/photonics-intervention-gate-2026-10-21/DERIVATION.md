# Correct only when correction improves expected retention

Date: October 21, 2026. Proposed engineering extension to Jason Duran Dutton's SEFI photonics program; not a new physical law or a uniquely SEFI consequence. Novelty has not been established.

## Why this rule

The photonics manuscript describes ideal restoration to a stability surface, but does not specify the disturbance introduced by the recovery actuator. A noisy correction can damage a state that was already close to its target. We can replace unconditional restoration by a decision based on a physically defined retention score. This is a mathematical engineering extension, not a derivation of the manuscript's curvature/torsion syndrome.

## Measured input and posterior

Use the previous pilot-assisted bosonic displacement model. beta is the unknown displacement and y is the mean of sixteen auxiliary heterodyne probe measurements after subtraction of their known input amplitude. For a calibrated Gaussian prior beta~N(0,S) and y=beta+epsilon with R=I/32:

    m = S(S+R)^(-1)y,
    C = (S^(-1)+R^(-1))^(-1).

Thus beta conditional on y is N(m,C). Neither beta nor the unknown payload is observed by the controller. Exact knowledge of S,R and actuator noise is a feasibility assumption shared by every policy.

Hold: residual distribution N(m,C). Correct: apply D(-m), but the actuator adds independent zero-mean displacement noise with covariance Q; residual distribution N(0,C+Q). Ideal bypass is assumed. Insertion loss, measurement backaction on probes/payload, drift and correlated actuator disturbances need separate treatment.

## Exact optical-encoding score

For the unprotected encoding span{|0>,|1>} and residual real vector z=(Re beta,Im beta), unconditional entanglement fidelity is

    f(z) = exp(-z^T z) (1-z^T z/2)^2.

This is a quantum overlap score, counts leakage as failure, and is not a protecting-code logical error rate. Define G(d,V)=E[f(z)] for z~N(d,V). For A(t)=I+2tV,

    g(t) = det(A(t))^(-1/2) exp[-t d^T A(t)^(-1) d],
    G(d,V) = g(1) + g'(1) + g''(1)/4.

The identity follows by differentiating the Gaussian Laplace transform under the integral: g'=-E[|z|² exp(-t|z|²)] and g''=E[|z|⁴ exp(-t|z|²)]. The implementation computes the derivatives via an eigendecomposition and checks them against independent Gauss-Hermite quadrature.

## Decision rule

    Apply correction if and only if G(0,C+Q) > G(m,C).

Otherwise bypass correction. Conditional on the calibrated model this maximizes expected fidelity over the two specified actions. Consequently its expected score cannot be below either fixed action under that model: it averages their pointwise maximum. This statement is conditional, not a universal superiority theorem. Finite-sample measured means may fluctuate; covariance mismatch can invalidate the benefit.

This is standard Bayesian decision theory applied to this encoding and actuator model. A conventional implementation with identical likelihood, score and actions is mathematically equivalent. There is no legitimate claim of superiority over that equivalent conventional rule.

An equal-information conventional baseline minimizes posterior residual mean-square error: it corrects iff ||m||² > tr(Q). Other baselines always correct or always bypass. All four actuator settings in protocol.json must be reported, including zero actuator disturbance. No favorable setting is selected as a confirmatory endpoint.

## Laboratory bridge

A group with auxiliary-probe sensing and displacement feed-forward can calibrate correction disturbance and compare alternate policies under controlled common-mode noise. The immediate hardware question is whether actuator-added noise is large enough for bypassing some corrections to help, after including bypass-path loss and controller overhead. If not, this intervention rule has little practical value on that apparatus.

A classical coherent-beam test can assess sensing, feed-forward and action fraction; it cannot establish retention of unknown quantum states. A quantum test needs specified source states, calibration of the auxiliary-probe/common-noise assumption, detector/actuator imperfections and unconditional tomography or another valid fidelity estimator. Lost trials must not be silently discarded.

The experiment would test cost-aware photonic control, not confirm the entire SEFI framework. A genuinely SEFI-specific experiment still requires operational identity variables and a distinct prediction not reducible to the established optical/control model.
