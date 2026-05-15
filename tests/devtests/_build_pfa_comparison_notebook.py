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

# Production first, then everything else. Unknown methods (if any) go to the end.
METHOD_ORDER = ['SortBisect', 'SortPair', 'LSA', 'BandMatch',
                'DesiredAgeBucket', 'GreedyOldEnough', 'KDTreeNN']

def order_methods(methods):
    known = [m for m in METHOD_ORDER if m in methods]
    extra = [m for m in methods if m not in METHOD_ORDER]
    return known + extra

print('Keys:', list(results.keys()))
print('Generic variants:', order_methods(sorted({k[0] for k in results.generic.keys()})))"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 1: wall-clock vs n_agents"))

cells.append(nbf.v4.new_code_cell("""\
rows = []
for (name, n, rep), v in results.generic.items():
    rows.append({'method': name, 'n_agents': n, 'rep': rep, 'wall_time': v.wall_time})
df = pd.DataFrame(rows)
agg = df.groupby(['method', 'n_agents'])['wall_time'].agg(['mean', 'std']).reset_index()

fig, ax = plt.subplots(figsize=(8, 5))
for method in order_methods(agg['method'].unique()):
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
                col='method', col_order=order_methods(df['method'].unique()),
                row='sex', kind='bar', errorbar='se',
                height=2.5, aspect=1.3)
g.set_axis_labels('Age', 'P(≥2 partners last year)')
g.set(ylim=(0, None))
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Figure 3a: incidence — pair age at formation (MF rel_types only)

Hexbin density of (male age, female age) for **every pair formed over the full
sim**. Each cell shows pair formation density for one method × rel_type. SW
relationships are excluded — they have different age dynamics and are shown
separately in fig 3c. Red dashed line = age equality; expected pattern is a
ridge **above** the diagonal (males ~7-8 yrs older than females per
`age_diff_pars`)."""))

cells.append(nbf.v4.new_code_cell("""\
MF_RT = ('stable', 'casual', 'onetime')

def _collect(records_field, rt_filter, source='generic'):
    src = getattr(results, source)
    all_methods_sorted = sorted({k[0] for k in src.keys()})
    if source == 'generic':
        n_max = max({k[1] for k in src.keys()})
        keep = lambda k: k[1] == n_max
    else:
        n_max = None
        keep = lambda k: True
    by = {(m, rt): [] for m in all_methods_sorted for rt in rt_filter}
    for k, v in src.items():
        if not keep(k):
            continue
        m = k[0]
        for r in v[records_field]:
            key = (m, r['rel_type'])
            if key in by:
                by[key].append((r['age_male'], r['age_female']))
    methods = order_methods([m for m in all_methods_sorted
                             if any(by[(m, rt)] for rt in rt_filter)])
    rel_types = [rt for rt in rt_filter
                 if sum(len(by[(m, rt)]) for m in methods) >= 50]
    return by, methods, rel_types, n_max

def _heatmap_grid(by, methods, rel_types, title, add_aggregate=True):
    # Optionally append an 'all MF' row that pools across rel_types per method.
    rows = list(rel_types) + (['all'] if add_aggregate and len(rel_types) > 1 else [])
    fig, axes = plt.subplots(len(rows), len(methods),
                             figsize=(2.5*len(methods), 2.8*len(rows)),
                             sharex=True, sharey=True, squeeze=False)
    for i, rt in enumerate(rows):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            if rt == 'all':
                ages = [pair for r in rel_types for pair in by[(method, r)]]
            else:
                ages = by[(method, rt)]
            if not ages:
                ax.text(0.5, 0.5, '(no records)', ha='center', va='center', transform=ax.transAxes)
                continue
            a1 = np.array([a[0] for a in ages])  # male
            a2 = np.array([a[1] for a in ages])  # female
            ax.hexbin(a1, a2, gridsize=25, cmap='viridis', mincnt=1)
            ax.plot([15, 80], [15, 80], 'r--', alpha=0.5, lw=1)
            mean_diff = float(np.mean(a1 - a2))
            ax.text(0.03, 0.97, f'Δ={mean_diff:+.1f}', transform=ax.transAxes,
                    ha='left', va='top', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7, lw=0))
            ax.set_xlim(15, 80); ax.set_ylim(15, 80)
            if i == 0: ax.set_title(method, fontsize=9)
            if j == 0: ax.set_ylabel(f'{rt}\\nfemale age')
            if i == len(rows)-1: ax.set_xlabel('male age')
    fig.suptitle(title, y=1.02)
    plt.tight_layout()
    plt.show()

by, methods, rel_types, n_max = _collect('pair_formation_ages', MF_RT)
print(f'n_agents={n_max}, MF rel_types kept:', rel_types)
print('Records per (method, rel_type):',
      {k: len(v) for k, v in by.items() if v})
_heatmap_grid(by, methods, rel_types, f'Incidence: pair formation ages, MF only, n={n_max}')"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Figure 3b: prevalence — currently active pairs at sim end (MF only)

Same hexbin grid but using a *snapshot* of pairs active at the final sim
timestep, with the agents' **current** ages. Long-lived pairs are over-
represented here relative to the incidence view, so age skew may attenuate."""))

cells.append(nbf.v4.new_code_cell("""\
by, methods, rel_types, n_max = _collect('pair_prevalence', MF_RT)
print('Records per (method, rel_type):',
      {k: len(v) for k, v in by.items() if v})
_heatmap_grid(by, methods, rel_types, f'Prevalence: active pairs at sim end, MF only, n={n_max}')"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Figure 3c: sex-work network age mixing (shown only if SW records exist)

The generic block uses bare `MFNetwork_*` variants which have no SW layer, so
this figure is empty for the generic results. For Zimbabwe / Eswatini blocks
the SW network is active and you'd see different (client × FSW) age dynamics
here."""))

cells.append(nbf.v4.new_code_cell("""\
by_inc, methods_inc, rt_inc, _ = _collect('pair_formation_ages', ('sw',))
by_prev, methods_prev, rt_prev, _ = _collect('pair_prevalence', ('sw',))
if rt_inc:
    _heatmap_grid(by_inc, methods_inc, rt_inc, 'SW incidence (formation)')
if rt_prev:
    _heatmap_grid(by_prev, methods_prev, rt_prev, 'SW prevalence (at sim end)')
if not rt_inc and not rt_prev:
    print('No SW edges in this run (generic block uses bare MFNetwork variants).')"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Figure 4: Zimbabwe age mixing (calibrated context)

The generic block runs with bare `MFNetwork_*` variants. Zimbabwe uses the
full `StructuredSexual` (MF + SW) and a calibrated HIV module. The same
incidence vs prevalence pattern should appear but is more pronounced because
the calibrated sim runs longer and has more pairs to draw on."""))

cells.append(nbf.v4.new_code_cell("""\
if not hasattr(results, 'zimbabwe') or not results.zimbabwe:
    print('No Zimbabwe results in this file.')
else:
    by, methods, rts, _ = _collect('pair_formation_ages', MF_RT, source='zimbabwe')
    _heatmap_grid(by, methods, rts, 'Zimbabwe MF incidence (pair formation)')
    by, methods, rts, _ = _collect('pair_prevalence', MF_RT, source='zimbabwe')
    _heatmap_grid(by, methods, rts, 'Zimbabwe MF prevalence (active at sim end)')
    by, methods, rts, _ = _collect('pair_formation_ages', ('sw',), source='zimbabwe')
    if rts:
        _heatmap_grid(by, methods, rts, 'Zimbabwe SW incidence (pair formation)',
                      add_aggregate=False)"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 4: lifetime partner distribution"))

cells.append(nbf.v4.new_code_cell("""\
fig, ax = plt.subplots(figsize=(8, 5))
n_max = max({k[1] for k in results.generic.keys()})
all_methods = order_methods(sorted({k[0] for k in results.generic.keys()
                                    if k[1] == n_max}))
for method in all_methods:
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
