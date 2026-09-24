"""Synthetic microwave calibration diagnostics. Requires Python >=3.10 and NumPy.

No hardware control, proprietary parameters, hidden-state access, or QEC claims.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path

import numpy as np

METHODS = ("Independent fit", "Conventional EWMA", "SEFI-inspired tangent EWMA")


def response(amplitude, detuning, settings, readout_error=0.02):
    """Driven two-level transition; time in 1/Omega0, detuning in Omega0.

    H/hbar = [(detuning-offset)*sigma_z + amplitude*sigma_x]*Omega0/2.
    Symmetric readout error is an explicit synthetic assumption.
    """
    a = np.asarray(amplitude)[..., None]
    d = np.asarray(detuning)[..., None] - settings[:, 0]
    w = np.sqrt(a * a + d * d)
    p = np.divide(a*a, w*w, out=np.zeros_like(w), where=w != 0)
    p = p * np.sin(w * settings[:, 1] / 2)**2
    return readout_error + (1 - 2*readout_error) * p


def design():
    return np.array([(offset, time) for offset in (-0.7, 0.0, 0.7)
                     for time in (1.1, 2.2, 3.3, 4.4)])


def smooth(values, alpha):
    state = np.zeros(values.shape[1])
    result = []
    for value in values:
        state = (1-alpha)*state + alpha*value
        result.append(state.copy())
    return np.asarray(result)


class Diagnostics:
    """Receives measurements only. Simulator truth is never passed here."""
    def __init__(self, shots=1024, alpha=0.3, settings=None):
        if not isinstance(shots, int) or shots < 1:
            raise ValueError("shots must be a positive integer")
        if not 0 < alpha <= 1:
            raise ValueError("alpha must be in (0,1]")
        self.shots, self.alpha = shots, alpha
        self.settings = design() if settings is None else np.asarray(settings, float)
        if (self.settings.ndim != 2 or self.settings.shape[1] != 2
                or not len(self.settings) or not np.isfinite(self.settings).all()
                or np.any(self.settings[:, 1] <= 0)):
            raise ValueError("settings require finite offset and positive pulse duration")
        self.reference = response(1.0, 0.0, self.settings)
        self.sigma = np.sqrt(self.reference*(1-self.reference)/shots)
        eps = 1e-5
        self.jacobian = np.stack([
            (response(1+eps, 0, self.settings)-response(1-eps, 0, self.settings))/(2*eps),
            (response(1, eps, self.settings)-response(1, -eps, self.settings))/(2*eps),
        ], axis=1)
        self.whitened = self.jacobian / self.sigma[:, None]
        self.identifiable = (np.linalg.matrix_rank(self.whitened) == 2
                             and np.linalg.cond(self.whitened) < 100)
        self.projector = np.linalg.pinv(self.whitened)
        self.covariance = np.linalg.pinv(self.whitened.T @ self.whitened)
        a, d = np.meshgrid(np.linspace(-.20, .20, 81), np.linspace(-.30, .30, 121))
        self.grid = np.column_stack((a.ravel(), d.ravel()))
        self.predictions = response(1+self.grid[:, 0], self.grid[:, 1], self.settings)
        self.weighted_grid = (self.predictions-self.reference)/self.sigma
        self.thresholds = None
        self.residual_limit = None

    def estimate(self, measured):
        measured = np.asarray(measured, float)
        if (measured.ndim != 2 or measured.shape[1] != len(self.settings)
                or not np.isfinite(measured).all()
                or np.any((measured < 0) | (measured > 1))):
            raise ValueError("measurements must be finite probabilities with one column per setting")
        z = (measured-self.reference)/self.sigma
        tangent = z @ self.projector.T
        # Distance expansion avoids a time x grid x setting allocation.
        distance = (np.sum(z*z, axis=1)[:, None]
                    + np.sum(self.weighted_grid**2, axis=1)[None, :]
                    - 2*z @ self.weighted_grid.T)
        fitted = self.grid[np.argmin(distance, axis=1)]
        estimates = np.stack((fitted, smooth(fitted, self.alpha), smooth(tangent, self.alpha)))
        scales = np.sqrt(np.maximum(np.diag(self.covariance), 1e-30))
        scores = np.sqrt(np.sum((estimates/scales)**2, axis=2))
        residual = np.sum((z-tangent @ self.whitened.T)**2, axis=1)
        return estimates, scores, residual

    def calibrate(self, rng, runs=40, steps=45):
        scores, residuals = [], []
        for _ in range(runs):
            observations = rng.binomial(self.shots, np.tile(self.reference, (steps, 1)))/self.shots
            _, score, residual = self.estimate(observations)
            scores.append(score)
            residuals.extend(residual.tolist())
        # Method-specific thresholds, identical healthy calibration data and target.
        self.thresholds = np.quantile(np.concatenate(scores, axis=1), .99, axis=1)
        self.residual_limit = float(np.quantile(residuals, .99))

    def diagnose(self, measured):
        if self.thresholds is None:
            raise ValueError("calibrate thresholds before diagnosis")
        if not self.identifiable:
            return dict(status="insufficient evidence", reason="The scan cannot locally distinguish amplitude and detuning.",
                        next_test="Add measurements at both positive and negative frequency offsets.")
        estimates, scores, residual = self.estimate(measured)
        alerts = scores > self.thresholds[:, None]
        boundaries = (np.abs(estimates[0, :, 0]) >= .199) | (np.abs(estimates[0, :, 1]) >= .299)
        ambiguous = (residual > self.residual_limit) | boundaries
        # Intervals are approximate local uncertainty, not calibrated coverage guarantees.
        local_sigma = np.sqrt(np.diag(self.covariance))
        recovery = []
        consecutive = 0
        for k in range(len(measured)):
            within = bool(np.all(np.abs(estimates[0, k])+2*local_sigma <= [.03, .04]))
            consecutive = consecutive+1 if within and not ambiguous[k] else 0
            recovery.append(consecutive >= 3)
        return dict(status="evaluated", estimates=estimates, scores=scores, alerts=alerts,
                    ambiguous=ambiguous, residual=residual, recovery=np.array(recovery),
                    local_sigma=local_sigma)


def truth_for(scenario, steps=45):
    if steps != 45:
        raise ValueError("this demonstration uses 45 fixed scan times")
    state = np.zeros((steps, 2))
    ramp = np.linspace(.1, 1, 20)
    if scenario in ("amplitude", "mixed"):
        state[10:30, 0] = .08*ramp
    if scenario in ("detuning", "mixed"):
        state[10:30, 1] = .12*ramp
    if scenario not in ("healthy", "amplitude", "detuning", "mixed"):
        raise ValueError("unknown scenario")
    return state


def simulate(state, model, rng):
    p = response(1+state[:, 0], state[:, 1], model.settings)
    return rng.binomial(model.shots, p)/model.shots


def benchmark(model, rng, runs):
    rows = []
    for scenario in ("healthy", "amplitude", "detuning", "mixed"):
        accum = [dict(squared=[], flags=[], delays=[], misses=0, recovery=[]) for _ in METHODS]
        truth = truth_for(scenario)
        for _ in range(runs):
            result = model.diagnose(simulate(truth, model, rng))
            for m in range(len(METHODS)):
                a = accum[m]
                a["squared"].extend(np.sum((result["estimates"][m]-truth)**2, axis=1).tolist())
                if scenario == "healthy":
                    a["flags"].extend(result["alerts"][m].tolist())
                else:
                    idx = np.flatnonzero(result["alerts"][m, 10:30])
                    if len(idx):
                        a["delays"].append(int(idx[0]))
                    else:
                        a["misses"] += 1
                    idx = np.flatnonzero(~result["alerts"][m, 30:])
                    if len(idx):
                        a["recovery"].append(int(idx[0]))
        for m, a in enumerate(accum):
            rows.append(dict(scenario=scenario, method=METHODS[m], runs=runs,
                vector_rmse=float(np.sqrt(np.mean(a["squared"]))),
                healthy_scan_false_alarm_rate=float(np.mean(a["flags"])) if a["flags"] else None,
                median_detection_delay_scans=float(np.median(a["delays"])) if a["delays"] else None,
                missed_fault_runs=a["misses"] if scenario != "healthy" else None,
                median_first_clear_delay_scans=float(np.median(a["recovery"])) if a["recovery"] else None))
    return rows


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_observations(path, observations, model):
    rows = [dict(scan=k, setting=j, offset_omega0=float(s[0]), duration_omega0=float(s[1]),
                 shots=model.shots, excited=int(round(p*model.shots)))
            for k, scan in enumerate(observations) for j, (s, p) in enumerate(zip(model.settings, scan))]
    write_csv(path, rows)


def load_observations(path, model):
    with Path(path).open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("empty measurement CSV")
    scans = {}
    for row in rows:
        k, j, n, count = (int(row[x]) for x in ("scan", "setting", "shots", "excited"))
        if k < 0 or not 0 <= j < len(model.settings) or n != model.shots or not 0 <= count <= n:
            raise ValueError("invalid scan, setting, shot count, or excited count")
        actual = [float(row["offset_omega0"]), float(row["duration_omega0"])]
        if not np.allclose(actual, model.settings[j], rtol=0, atol=1e-10):
            raise ValueError("measurement settings differ from calibrated design")
        if j in scans.setdefault(k, {}):
            raise ValueError("duplicate setting in scan")
        scans[k][j] = count/n
    if sorted(scans) != list(range(len(scans))) or any(len(s) != len(model.settings) for s in scans.values()):
        raise ValueError("scans must be contiguous, complete, and start at zero")
    return np.array([[scans[k][j] for j in range(len(model.settings))] for k in range(len(scans))])


def service_report(result, name):
    if result["status"] != "evaluated":
        return f"# {name}\n\nInsufficient evidence. {result['reason']}\n\nNext test: {result['next_test']}\n"
    lines = [f"# {name}", "", "SYNTHETIC DEMONSTRATION — illustrative parameters; not eleQtron hardware data.", "",
             "No automatic correction is applied. Recommendations concern measurements, not component replacement.", "",
             "| Scan | Amplitude error estimate | Detuning / Omega0 | Evidence | Next action |", "|---|---:|---:|---|---|"]
    for k in range(len(result["recovery"])):
        a, d = result["estimates"][0, k]
        if result["ambiguous"][k]:
            status, action = "Insufficient model agreement", "Repeat scan; check readout/reference; broaden physical model before diagnosis"
        elif result["alerts"][1, k]:
            status = "Persistent drift indication"
            action = "Cross-check delivered microwave amplitude and resonance with independent measurements"
        elif result["recovery"][k]:
            status, action = "Within demo acceptance bands for 3 scans", "Record evidence; validate hardware-specific acceptance criteria"
        else:
            status, action = "Monitor / recovery not yet verified", "Acquire another complete scan"
        lines.append(f"| {k} | {a:+.3f} | {d:+.3f} | {status} | {action} |")
    lines += ["", "Acceptance bands: amplitude +/-3%, detuning +/-0.04 Omega0, including approximate 2-sigma local uncertainty; three consecutive scans; residual gate must pass. These are demonstration choices, not device specifications.",
              "", "Recovery assessment uses independent fits to avoid hiding residual faults through smoothing. A clear alarm alone is not proof of recovery. The report cannot localize detuning to a magnetic source or a microwave source."]
    return "\n".join(lines)+"\n"


def chart_svg(truth, result, component):
    colors = ["#9daebc", "#4ea8de", "#eea64c", "#6fe0ba"]
    series = [truth[:, component]] + [x[:, component] for x in result["estimates"]]
    lo, hi = -.06, .16
    paths = []
    for values, color in zip(series, colors):
        points = " ".join(f"{40+i*680/44:.1f},{210-(v-lo)*180/(hi-lo):.1f}" for i, v in enumerate(values))
        paths.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>')
    return '<svg viewBox="0 0 760 250" role="img" aria-label="Synthetic truth and estimated drift over 45 scans">' + \
        '<path d="M40 20V210H720" fill="none" stroke="#8090a0"/>' + \
        '<text x="40" y="235" fill="#bbb">Scan 0</text><text x="620" y="235" fill="#bbb">Scan 44</text>' + \
        f'<text x="2" y="35" fill="#bbb">{hi}</text><text x="2" y="210" fill="#bbb">{lo}</text>' + ''.join(paths) + '</svg>'


def write_html(path, rows, truth, result, config):
    cells = []
    for row in rows:
        values = [row["scenario"], row["method"], f"{row['vector_rmse']:.4f}",
                  "—" if row["healthy_scan_false_alarm_rate"] is None else f"{100*row['healthy_scan_false_alarm_rate']:.2f}%",
                  row["median_detection_delay_scans"], row["missed_fault_runs"]]
        cells.append('<tr>'+''.join('<td>'+html.escape(str(v) if v is not None else '—')+'</td>' for v in values)+'</tr>')
    path.write_text('''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>SEFI calibration workbench</title><style>body{background:#101923;color:#e6edf3;font:17px/1.6 system-ui;margin:40px auto;padding:0 24px;max-width:1150px}h1{font-size:38px;line-height:1.2}h2{margin-top:36px}.tag{color:#6fe0ba;letter-spacing:2px}.box{padding:20px;background:#1b2937;border-radius:12px;margin:20px 0}svg{width:100%;max-width:850px}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:10px;text-align:left;border-bottom:1px solid #344454}.scroll{overflow-x:auto}a{color:#81cfff}small{color:#b9c6d3}</style>
<div class="tag">FIELD ENGINEERING / RESEARCH PROTOTYPE</div><h1>Detect drift. Test the cause.<br>Verify recovery.</h1>
<div class="box"><strong>Synthetic data only.</strong> This demonstrates a diagnostic workflow, not eleQtron performance, a hardware digital twin, or validated quantum error correction.</div>
<p>A controlled microwave-transition model produces repeated measurements. Three methods see those measurements; hidden fault labels are reserved for evaluation.</p>
<h2>Mixed-fault demonstration</h2><p>Fault ramps during scans 10–29; synthetic restoration starts at scan 30. Gray: hidden truth (evaluation only). Blue: independent fit. Orange: conventional EWMA. Green: SEFI-inspired tangent EWMA.</p>
<h3>Fractional amplitude error</h3>''' + chart_svg(truth, result, 0) + '<h3>Detuning / nominal Rabi angular frequency</h3>' + chart_svg(truth, result, 1) + '''
<h2>Held-out benchmark</h2><p>All methods receive the same scans. Each threshold is set from separate healthy calibration runs at a nominal 1% per-scan false-alarm target. Detection delays count scans after fault onset; medians exclude missed runs.</p><div class="scroll"><table><tr><th>Scenario</th><th>Method</th><th>Vector RMSE</th><th>Healthy false alarms</th><th>Median delay</th><th>Missed runs</th></tr>''' + ''.join(cells) + '''</table></div>
<p>RMSE combines two dimensionless coordinates and includes onset and recovery. Temporal averaging introduces lag. Observed differences do not establish novelty or a SEFI advantage.</p>
<h2>When evidence is insufficient</h2><div class="box">A resonant-only scan has no first-order sensitivity to the sign of detuning at the reference. The prototype abstains and recommends scans at positive and negative frequency offsets.</div>
<h2>Service workflow</h2><p>Capture settings and counts → compare with reference → flag drift or ambiguity → cross-check independent measurements → verify three consecutive acceptable scans. An amplitude or detuning estimate does not identify a failed physical component.</p>
<p><a href="service_report.md">Mixed-fault service report</a> · <a href="ambiguous_report.md">Ambiguous-case report</a> · <a href="metrics.json">Machine-readable results</a></p>
<h2>Scope and assumptions</h2><p>Two-level coherent response, fixed 2% symmetric readout error, independent binomial shots, drift constant within each scan, trusted initial reference and configured settings. No phase estimation, leakage, heating, colored noise, real instrument integration, or storage-transfer dynamics. Real use needs reference calibration and model validation.</p><small>''' + html.escape(json.dumps(config, sort_keys=True)) + '</small></html>', encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("demo_output"))
    parser.add_argument("--runs", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20260925)
    parser.add_argument("--shots", type=int, default=1024)
    parser.add_argument("--input", type=Path, help="Replay CSV with the exact demo settings; uses synthetic reference calibration")
    args = parser.parse_args()
    if args.runs < 2:
        parser.error("runs must be at least 2")
    args.output.mkdir(parents=True, exist_ok=True)
    model = Diagnostics(shots=args.shots)
    calibration_rng, evaluation_rng, demo_rng = [np.random.default_rng(s) for s in np.random.SeedSequence(args.seed).spawn(3)]
    model.calibrate(calibration_rng)
    config = dict(seed=args.seed, evaluation_runs_per_scenario=args.runs, calibration_runs=40,
                  scans_per_run=45, shots_per_setting=args.shots, settings_per_scan=12, alpha=model.alpha,
                  data="synthetic", numpy_version=np.__version__, readout_error=.02,
                  thresholds=model.thresholds.tolist(), residual_limit=model.residual_limit)
    if args.input:
        result = model.diagnose(load_observations(args.input, model))
        report = service_report(result, "CSV replay diagnostic report").replace(
            "SYNTHETIC DEMONSTRATION — illustrative parameters; not eleQtron hardware data.",
            "CSV REPLAY — provenance supplied by caller; reference and thresholds are synthetic. Not hardware-validated.")
        (args.output/"replay_report.md").write_text(report, encoding="utf-8")
        print(args.output/"replay_report.md")
        return
    rows = benchmark(model, evaluation_rng, args.runs)
    truth = truth_for("mixed")
    observed = simulate(truth, model, demo_rng)
    result = model.diagnose(observed)
    save_observations(args.output/"measurements.csv", observed, model)
    write_csv(args.output/"evaluation_truth.csv", [dict(scan=k, amplitude_error=a, detuning_omega0=d) for k, (a, d) in enumerate(truth)])
    write_csv(args.output/"metrics.csv", rows)
    (args.output/"metrics.json").write_text(json.dumps(dict(config=config, results=rows), indent=2)+"\n", encoding="utf-8")
    (args.output/"service_report.md").write_text(service_report(result, "Mixed-fault service report"), encoding="utf-8")
    sparse = Diagnostics(settings=np.array([(0, t) for t in (1.1, 2.2, 3.3, 4.4)]))
    sparse.thresholds, sparse.residual_limit = model.thresholds, model.residual_limit
    ambiguous = sparse.diagnose(np.tile(sparse.reference, (3, 1)))
    (args.output/"ambiguous_report.md").write_text(service_report(ambiguous, "Ambiguous measurement design"), encoding="utf-8")
    write_html(args.output/"index.html", rows, truth, result, config)
    print(json.dumps(dict(output=str(args.output), results=rows), indent=2))


if __name__ == "__main__":
    main()
