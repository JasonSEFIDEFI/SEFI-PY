# Completed conventional controls — 23 October 2026

The extended retrospective benchmark completed in 766.72 seconds on the recorded environment. It evaluated six model families over 17 partitions, producing 306 interval cells and 107,820 prediction records (interval methods repeat the same point predictions). Eight analysis/GP tests passed. Protocol, explicit split IDs, selections, predictions and environment are in `results/publication/`.

Mean held-out binding-energy RMSE in MeV (equal weight per partition within each row):

| Partition | SEMF | Gaussian residual smoother | GP RBF | GP Matern32 | GP RBF parity |
|---|---:|---:|---:|---:|---:|
| Random, 5 splits | 3.3395 | 1.1834 | 0.4282 | 0.4582 | 0.5212 |
| Whole element, 5 folds | 3.3268 | 1.3073 | 0.5488 | 0.5448 | 0.6909 |
| Proton extrapolation | 12.7649 | 12.7116 | 12.7351 | 12.7411 | 12.7551 |
| Neutron frontier | 5.0808 | 3.4806 | 1.0236 | 1.2103 | 1.0851 |
| Historical 2016-to-2020, 5 splits | 5.2079 | 2.8737 | 0.9631 | 0.9849 | 0.8764 |

Conventional Gaussian processes substantially outperform the pilot Gaussian residual smoother on several partitions. This does not establish a unique SEFI advantage. Heavy-proton extrapolation remains poor for every model. Hyperparameters were selected on tuning labels, not test labels; the protocol was fixed before this run but is not an external preregistration. Data and historical split limitations from the pilot remain applicable.

The older manuscript and plots in this archive predate these stronger controls and are preserved for provenance. They have not been revised into a submission-ready paper. Manuscript preparation was paused at the author's request in favor of mathematics. The current critical-binding derivation and numerical checks are in the sibling `critical-binding/` directory.

Uncertainty intervals and their observed coverage are recorded in the raw metrics. These average RMSEs are descriptive summaries, not confidence intervals. Full correlated model comparisons and external mass-model comparisons remain future work. Timing fields for the known-shell SEMF variant contain a zero placeholder; do not interpret that as a zero computational cost. GP covariance intervals omit uncertainty in the fitted SEMF mean.
