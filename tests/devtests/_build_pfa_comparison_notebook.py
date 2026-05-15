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

cells.append(nbf.v4.new_markdown_cell("## Figure 2: partners in last year by age, sex, rel_type"))

cells.append(nbf.v4.new_code_cell("""\
rows = []
n_max = max({k[1] for k in results.generic.keys()})
for (name, n, rep), v in results.generic.items():
    if n != n_max:
        continue
    for r in v.partners_last_year:
        rows.append({**r, 'method': name})
df = pd.DataFrame(rows)
df['age_bin'] = pd.cut(df['age'], bins=[0, 20, 25, 30, 40, 50, 100], labels=['<20', '20-25', '25-30', '30-40', '40-50', '50+'])

g = sns.catplot(data=df, x='age_bin', y='n_partners', hue='rel_type',
                col='method', row='sex', kind='bar', errorbar='se',
                height=2.5, aspect=1.3)
g.set_axis_labels('Age', 'Mean partners (last year)')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 3: male vs female age heatmaps per rel_type"))

cells.append(nbf.v4.new_code_cell("""\
methods = sorted({k[0] for k in results.generic.keys()})
rel_types = ['stable', 'casual', 'onetime']
n_max = max({k[1] for k in results.generic.keys()})

fig, axes = plt.subplots(len(methods), len(rel_types),
                         figsize=(3*len(rel_types), 2.5*len(methods)),
                         sharex=True, sharey=True)
for i, method in enumerate(methods):
    records = []
    for (m, n, rep), v in results.generic.items():
        if m == method and n == n_max:
            records.extend(v.pair_age_heatmap)
    for j, rt in enumerate(rel_types):
        ax = axes[i, j]
        ages = [(r['age_p1'], r['age_p2']) for r in records if r['rel_type'] == rt]
        if not ages:
            ax.set_visible(False)
            continue
        a1 = np.array([a[0] for a in ages])
        a2 = np.array([a[1] for a in ages])
        ax.hexbin(a1, a2, gridsize=20, cmap='viridis', mincnt=1)
        ax.plot([15, 80], [15, 80], 'r--', alpha=0.4, lw=0.8)
        if i == 0:
            ax.set_title(rt)
        if j == 0:
            ax.set_ylabel(f'{method}\\nfemale age')
        if i == len(methods)-1:
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
