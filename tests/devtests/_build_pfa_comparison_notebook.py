"""One-shot builder for pfa_comparison.ipynb. Run once, then keep for regeneration."""
import nbformat as nbf
from pathlib import Path

nb = nbf.v4.new_notebook()

cells = []

cells.append(nbf.v4.new_markdown_cell("""\
# PFA comparison

Compare seven pair-formation algorithm variants by runtime and network structure.
See `docs/superpowers/specs/2026-05-15-pfa-comparison-design.md` for the spec.

Run `python devtest_pfa_comparison.py` first to produce `pfa_comparison_results.obj`."""))

cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt
import seaborn as sns

results = sc.load('pfa_comparison_results.obj')
print('Keys:', list(results.keys()))
print('Generic variants:', sorted({k[0] for k in results.generic.keys()}))"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 1: wall-clock vs n_agents"))

cells.append(nbf.v4.new_code_cell("""\
rows = []
for (name, n, rep), v in results.generic.items():
    rows.append({'method': name, 'n_agents': n, 'rep': rep, 'wall_time': v.wall_time})
df = pd.DataFrame(rows)
agg = df.groupby(['method', 'n_agents'])['wall_time'].agg(['mean', 'std']).reset_index()

fig, ax = plt.subplots(figsize=(8, 5))
for method in sorted(agg['method'].unique()):
    sub = agg[agg['method'] == method]
    ax.errorbar(sub['n_agents'], sub['mean'], yerr=sub['std'].fillna(0), marker='o', label=method)
ax.set_xlabel('n_agents')
ax.set_ylabel('wall time (s)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
ax.set_title('Wall time vs n_agents')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Figure 2: concurrency rate by age, sex, rel_type

For each (method, sex, age_bin, rel_type) cell, fraction of partnered agents with **≥2
partners** in the last year. Mean is unhelpful here because almost everyone with a
partner has exactly one — the interesting signal is the tail."""))

cells.append(nbf.v4.new_code_cell("""\
rows = []
n_max = max({k[1] for k in results.generic.keys()})
for (name, n, rep), v in results.generic.items():
    if n != n_max:
        continue
    for r in v.partners_last_year:
        rows.append({**r, 'method': name})
df = pd.DataFrame(rows)
df['age_bin'] = pd.cut(df['age'], bins=[0, 20, 25, 30, 40, 50, 100],
                      labels=['<20', '20-25', '25-30', '30-40', '40-50', '50+'])
df['concurrent'] = (df['n_partners'] >= 2).astype(int)

g = sns.catplot(data=df, x='age_bin', y='concurrent', hue='rel_type',
                col='method', row='sex', kind='bar', errorbar='se',
                height=2.5, aspect=1.3)
g.set_axis_labels('Age', 'P(≥2 partners last year)')
g.set(ylim=(0, None))
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Figure 3: male vs female age heatmaps per rel_type

Hexbin density of (male age, female age) for pairs formed in the last year. Red dashed
line = age equality. Drops 'onetime' because at the tested scales it has too few records
to plot meaningfully. Aggregates across reps to get smoother density."""))

cells.append(nbf.v4.new_code_cell("""\
all_methods = sorted({k[0] for k in results.generic.keys()})
n_max = max({k[1] for k in results.generic.keys()})

# Aggregate (age_p1, age_p2) per (method, rel_type) at n_max.
records_by = {(m, rt): [] for m in all_methods for rt in ('stable', 'casual', 'onetime')}
for (m, n, rep), v in results.generic.items():
    if n != n_max:
        continue
    for r in v.pair_age_heatmap:
        key = (m, r['rel_type'])
        if key in records_by:
            records_by[key].append((r['age_p1'], r['age_p2']))

# Drop methods with no records at n_max (e.g. LSA, which is capped at n<n_max).
methods = [m for m in all_methods if any(records_by[(m, rt)] for rt in ('stable', 'casual', 'onetime'))]
rel_types = [rt for rt in ('stable', 'casual', 'onetime')
             if sum(len(records_by[(m, rt)]) for m in methods) >= 50]

# Layout: rel_type rows × method columns -- landscape.
fig, axes = plt.subplots(len(rel_types), len(methods),
                         figsize=(2.5*len(methods), 2.8*len(rel_types)),
                         sharex=True, sharey=True, squeeze=False)
for i, rt in enumerate(rel_types):
    for j, method in enumerate(methods):
        ax = axes[i, j]
        ages = records_by[(method, rt)]
        if not ages:
            ax.text(0.5, 0.5, '(no records)', ha='center', va='center', transform=ax.transAxes)
            continue
        a1 = np.array([a[0] for a in ages])
        a2 = np.array([a[1] for a in ages])
        ax.hexbin(a1, a2, gridsize=25, cmap='viridis', mincnt=1)
        ax.plot([15, 80], [15, 80], 'r--', alpha=0.5, lw=1)
        ax.set_xlim(15, 80)
        ax.set_ylim(15, 80)
        if i == 0:
            ax.set_title(method, fontsize=9)
        if j == 0:
            ax.set_ylabel(f'{rt}\\nfemale age')
        if i == len(rel_types)-1:
            ax.set_xlabel('male age')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 4: lifetime partner distribution"))

cells.append(nbf.v4.new_code_cell("""\
fig, ax = plt.subplots(figsize=(8, 5))
n_max = max({k[1] for k in results.generic.keys()})
for method in methods:
    vals = []
    for (m, n, rep), v in results.generic.items():
        if m == method and n == n_max:
            vals.extend(v.lifetime_partners[v.lifetime_partners > 0])
    if vals:
        ax.hist(vals, bins=np.arange(0, 30), histtype='step', density=True, label=method, lw=2)
ax.set_xlabel('lifetime partners')
ax.set_ylabel('density')
ax.legend()
ax.set_title(f'Lifetime partner distribution (n={n_max})')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Takeaway

Compare:

- **Speed:** which variants stay tractable as n grows? Where does LSA become unusable?
- **Age mixing:** which variants put more mass off the diagonal in fig 3?
- **Concurrency:** which produce wider lifetime-partner tails?

The production default is SortBisect. Departures from it in figs 2-4 indicate the structural
cost (or benefit) of a different algorithmic choice."""))

nb.cells = cells

out = Path(__file__).parent / 'pfa_comparison.ipynb'
with out.open('w') as f:
    nbf.write(nb, f)
print(f'Wrote {out}')
