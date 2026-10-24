"""Plot saved results only; no simulated or fitted observations are added."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
rows = json.loads((root/'results/profiles.json').read_text())
branches = json.loads((root/'results/branch_summary.json').read_text())
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
for h, marker in ((.4, 'o'), (.2, 'x')):
    selected = sorted((r for r in rows if r['h'] == h), key=lambda r: r['Q'])
    axes[0].plot([r['Q']/1000 for r in selected], [r['energy_per_charge'] for r in selected], marker=marker, label=f'h={h}')
    branch = next(b for b in branches if b['h'] == h)
    axes[1].plot([r['parent_Q']/1000 for r in branch['selected_fission']],
                 [r['daughter_energy_minus_parent'] for r in branch['selected_fission']], marker=marker, label=f'h={h}')
axes[0].axhline(np.sqrt(2), ls='--', color='gray', label='free Phi threshold')
axes[0].axhline(np.sqrt(1+3/np.sqrt(12)), ls=':', color='black', label='potential lower bound')
axes[0].set(xlabel='Conserved charge / 1000', ylabel='Energy / charge', title='Specified scalar model, lambda=12')
axes[1].axhline(0, color='gray', ls='--')
axes[1].set(xlabel='Parent charge / 1000', ylabel='Daughter energy minus parent energy', title='Four selected same-branch splits')
for ax in axes:
    ax.grid(alpha=.2); ax.legend(fontsize=8)
fig.savefig(root/'branch-results.png', dpi=180)
