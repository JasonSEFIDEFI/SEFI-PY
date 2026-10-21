# Results: adaptive displacement recovery

In a frozen synthetic pilot-assisted displacement-recovery benchmark, a conventional adaptive mixture improved mean entanglement fidelity over a development-selected 128-sample rolling estimator by 0.210 percentage points in the prespecified drifting-noise scenario (paired seed-bootstrap 95% interval 0.177 to 0.244 points; 20 independent streams). It did not outperform the tested EW estimator in that scenario and was slower in the measured Python implementation. This does not establish a SEFI-specific advantage.

## Fidelity

| Condition | Fixed | EW | Rolling | Mixture |
|---|---:|---:|---:|---:|
| Stationary | 0.964754 | 0.960610 | 0.962630 | 0.963237 |
| Mean steps | 0.826256 | 0.958464 | 0.957767 | 0.960683 |
| Rotation | 0.958703 | 0.960915 | 0.962702 | 0.962949 |
| Joint drift | 0.798441 | 0.949550 | 0.946370 | 0.948473 |
| Bursts | 0.940898 | 0.936684 | 0.938391 | 0.940183 |

## Timing

| Method | Median microseconds per observation |
|---|---:|
| Fixed prior | 5.62 |
| EW (rate 0.03) | 15.37 |
| Rolling (128) | 75.57 |
| Adaptive mixture | 108.08 |

All 100 cases completed; 240,000 observations, 20 independent seeds per condition. Seven scientific checks passed and 600 saved seed/method scores independently reproduced. Primary positive result does not establish superiority over the EW baseline or a SEFI-specific mechanism. Full details: Adaptive_Recovery_Report.pdf.
