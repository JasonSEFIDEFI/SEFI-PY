# Adaptive pilot-assisted displacement recovery

Author project: Jason Duran Dutton. Protocol frozen October 20, 2026 before development/evaluation output inspection. This is an internal protocol freeze, not an independent preregistration.

## Physical and information setting

The payload is an unknown qubit in span{|0>,|1>}, not a protecting QEC code. Sixteen auxiliary coherent probe modes share the payload displacement beta. Their heterodyne sample mean gives y=beta+epsilon, with R=I/32. Every estimator receives identical noisy calibration records and current/past observations. Truth and regime metadata are used exclusively in the evaluator. Corrections are the ideal oscillator displacement D(-estimate). No actuator noise, loss channel, detector hardware, or deadline is simulated.

## Conventional baseline estimators

Fixed: calibration sample mean and sample covariance minus R, projected onto positive covariance matrices. ML: estimate=y. EW: causal exponentially weighted observation mean/covariance, with covariance deconvolution by subtracting R. Rolling: the analogous estimate from a trailing window. All adaptive moments are updated after estimating the current displacement. Covariance eigenvalues have floor 1e-8.

For prior mean m and covariance S, estimate=m+S(S+R)^(-1)(y-m). This is conventional Gaussian empirical Bayes. With adaptive estimated moments, it is an approximate model-based estimator, not an exact posterior for every scenario, and is not proven optimal for fidelity.

## Candidate

Three EW experts have fixed rates 0.01, 0.05, 0.15. Before updating their moments, each produces an estimate and Gaussian predictive likelihood p_j(y)=N(y;m_j,S_j+R). Posterior expert weights are proportional to previous weights times p_j(y). The correction estimate is their weighted mean. Weights are then smoothed with fixed share 0.005 toward uniform. All experts subsequently update using y.

This is a conventional predictive mixture. The SEFI project motivates studying identity retention and reference geometry, but no unique SEFI deduction, new physical law, or algorithmic novelty is asserted. Introducing additional identity variables would require separate operational definitions and ablation experiments.

## Selection and evaluation

The EW rate and rolling-window baseline are selected only on four development seeds over five equally weighted scenarios. The stronger selected family becomes the primary comparator. The fixed mixture candidate is not tuned. Evaluation uses twenty separate seeds and 2400 observations per scenario per seed. Known scenario families are shared across development/evaluation; generalization to unknown noise families is not tested.

Primary comparison: paired difference in unconditional entanglement fidelity in joint_drift, averaged within each seed, then a 10,000-resample paired seed bootstrap. CI is 95%; remaining comparisons are descriptive and not multiple-testing-adjusted. A positive lower bound supports only this bounded comparison. It does not establish superiority to all conventional methods.

Fidelity for residual r=beta-estimate is exp(-|r|^2)(1-|r|^2/2)^2. This equals |Tr(V^dagger D(r)V)|^2/4; leakage counts as failure. Secondary metrics include residual mean-square error and residual-radius failure fraction. These are simulation scores, not observed hardware logical error rates.

Timing includes estimator initialization and every prediction/update on a complete stream, excludes simulation and disk I/O, and uses seven alternating-order repetitions. Methods share inputs and machine; they do not receive equal computation by construction. Thus fidelity superiority alone cannot establish superiority at a matched compute budget. No real-time or faster-at-matched-quality claim follows without additional measurements.

## Interruptions and reproducibility

Each stage records started/completed/interrupted/failed status, source and protocol hashes, completed cases, and environment. Incomplete stages cannot support a final claim. Full truth, observations, calibrations, and estimates are saved for all evaluation cases. Original console equations and prior benchmarks are unchanged.
