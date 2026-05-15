# PFA comparison: design spec

**Date:** 2026-05-15
**Branch:** `feat/pfa-comparison` off `rc1.5.5`
**Owner:** Robyn

## Goal

Compare seven pair-formation algorithms (PFAs) that can override
`MFNetwork.match_pairs` in `stisim.networks`. Measure runtime scaling
and network-structure diagnostics on (a) a generic STIsim sim with
`StructuredSexual`, and (b) two realistic, calibrated localizations
(`hivsim.demo('zimbabwe')` and `hivsim_eswatini`).

This is dev work; outputs are devtest scripts, a Jupyter notebook, and
saved result objects — not formal docs.

## Methods

All seven are subclasses of `MFNetwork` that override `match_pairs` only.
Female desired ages are sampled the same way in every variant (existing
`get_age_risk_pars` + `age_diffs` logic); methods differ only in how
they pair `f_looking` to `m_eligible`.

| # | Class | Algorithm |
|---|-------|-----------|
| 1 | `MFNetwork_LSA` | Full distance matrix + `scipy.optimize.linear_sum_assignment`. O(n³). |
| 2 | `MFNetwork_SortPair` | `np.argsort` both groups, zip, truncate to `min(len)`. No bisect trim. |
| 3 | `MFNetwork_DesiredAgeBucket` | Round desired age to integer year. For each year-bucket, match women to available men; if M < W in a bucket, sample men **with replacement**. Post-filter (in `match_pairs` before returning): (i) drop matches that put the male over his `concurrency` cap, processing women in random order so the cap is enforced fairly; (ii) drop matches whose risk-group combination would not form a relationship under `p_matched_stable` / `p_mismatched_casual` (i.e. apply the relationship-acceptance Bernoulli draw inside `match_pairs` rather than in `add_pairs`). |
| 4 | `MFNetwork_SortBisect` | **Current production** (`MFNetwork.match_pairs`): argsort + bisect-trim the head and tail to match age supports, then subsample the larger group. |
| 5 | `MFNetwork_GreedyOldEnough` | Sort women by desired age ascending. Walk through: for each woman, pop the youngest available man with `age >= desired_age`. No replacement. Women who can't match are skipped. |
| 6 | `MFNetwork_KDTreeNN` | Build `scipy.spatial.KDTree` on male ages. For each woman, query k=1 nearest. Resolve collisions greedily by age-error ascending; losers are dropped. |
| 7 | `MFNetwork_BandMatch` | Bucket both groups into 5-yr age bands. Within each band, shuffle and zip; truncate to band-min. Fast, coarse. |

### Common variant API

Each variant is a thin subclass:

```python
class MFNetwork_LSA(sti.MFNetwork):
    def match_pairs(self):
        # ... method-specific implementation ...
        return p1, p2  # ss.uids for males and females respectively
```

A `StructuredSexual_LSA` etc. is **not** created — variants are used
via direct instantiation of the MF subclass, since `add_pairs` in the
base class calls `self.match_pairs` polymorphically.

For sex-work behavior, `SWNetwork` is unchanged (orthogonal to MF
matching). Comparison can be MF-only or MF+SW; default is MF-only to
isolate the algorithmic change.

## File layout

```
stisim/
  stisim/
    networks/
      pfa_variants.py        # NEW. Seven MFNetwork subclasses.
  tests/
    devtests/
      devtest_pfa_comparison.py    # NEW. Runs all 7 at n=1k/5k/10k.
      pfa_comparison.ipynb         # NEW. Plots + writeup.
hivsim_eswatini/
  run_pfa_comparison.py            # NEW. Calibrated comparison.
```

Note: `stisim/networks.py` is currently a single file, not a package.
`pfa_variants.py` will be created as a sibling module
(`stisim/stisim/pfa_variants.py`) and imported in `stisim/__init__.py`.
The "networks/" subdirectory above is a logical grouping, not a
required refactor.

## Devtest script (`devtest_pfa_comparison.py`)

Pseudocode:

```python
N_AGENTS = [1_000, 5_000, 10_000]
N_REPS = 3
SIM_YEARS = 20
PFA_CLASSES = [LSA, SortPair, DesiredAgeBucket, SortBisect,
               GreedyOldEnough, KDTreeNN, BandMatch]

results = sc.objdict()
for n in N_AGENTS:
    for pfa_cls in PFA_CLASSES:
        for rep in range(N_REPS):
            net = pfa_cls()
            sim = ss.Sim(
                n_agents=n,
                start='1990-01-01', stop=f'{1990+SIM_YEARS}-01-01',
                networks=net,
                diseases='sis',   # cheap disease for warm sim
                analyzers=[
                    sti.partner_age_diff(year=...),
                    PartnersLastYearAnalyzer(),   # NEW
                    PairAgeHeatmapAnalyzer(),     # NEW
                ],
                rand_seed=rep,
            )
            with sc.timer() as t:
                sim.run()
            results[(pfa_cls.__name__, n, rep)] = dict(
                wall_time=t.elapsed,
                analyzers=sim.analyzers,
                lifetime_partners=net.lifetime_partners.values,
            )

sc.save('pfa_comparison_results.obj', results)
```

### New analyzers (added in same file, not promoted)

1. **`PartnersLastYearAnalyzer`** — at sim end, for each agent, counts
   unique partners in last 365 days by `rel_type` (m / c / o). Stores
   `(age, sex, rel_type, n_partners)` tuples for plotting as
   age-stratified bar charts.
2. **`PairAgeHeatmapAnalyzer`** — at sim end, for each rel_type, store
   `(age_p1, age_p2)` pairs from active edges over the last year.

`partner_age_diff` is reused from sti analyzers.
`net.lifetime_partners` is already a state on `MFNetwork`.

## Jupyter notebook (`pfa_comparison.ipynb`)

Single notebook, cells in order:

1. Intro (1 markdown cell — 4 sentences)
2. Method descriptions (1 markdown cell each — pseudocode + 2-3 lines)
3. Load `pfa_comparison_results.obj`
4. **Figure 1**: wall-clock vs `n_agents` (one line per method)
5. **Figure 2**: partners-in-last-year by age × sex × rel-type
   (7 methods × 3 rel-types small-multiple grid)
6. **Figure 3**: male-vs-female age heatmaps per rel-type, per method
7. **Figure 4**: lifetime-partner distribution (KDE or histogram, 7 lines)
8. Takeaway (1 markdown cell)

## Calibrated comparisons

### `hivsim.demo('zimbabwe')`

Add a small block at the end of `devtest_pfa_comparison.py`:

```python
for pfa_cls in PFA_CLASSES:
    sim = hivsim.demo('zimbabwe', run=False, n_agents=10_000)
    # swap the MF network for the variant
    sim.pars.networks = [pfa_cls(...), ...]
    sim.run()
    # save HIV prevalence + same network diagnostics
```

### `hivsim_eswatini/run_pfa_comparison.py`

Imports the PFA variant classes from stisim. Builds the standard
Eswatini sim, swaps the MF network, runs each variant, saves
network diagnostics + key HIV outcomes (incidence, prevalence, ART
coverage by age/sex). Plots overlay on existing calibration figures
so we can see if PFA choice meaningfully shifts calibrated
quantities.

## Out of scope

- Promoting `pfa_variants.py` to a top-level export
- Tuning method 3's bucket width or post-filter beyond defaults
- Documenting findings as `.qmd` slides (this is dev work)
- Modifying `SWNetwork.match_pairs`
- Touching CHANGELOG (per [feedback_changelog_timing](../../../../../.claude/projects/-Users-robynstuart-gf-stisim/memory/feedback_changelog_timing.md))

## Open risks

- Method 1 (LSA) at n=10k builds a 100M-cell distance matrix and runs
  O(n³). May exceed 1 min per timestep and many GB of RAM. Plan: cap
  LSA at n=5k unless explicitly enabled; document the cap.
- Method 3's post-filter on relationship-type compatibility duplicates
  logic that currently sits in `add_pairs`. Either factor that logic
  into a helper or accept the duplication in the variant.
- Diagnostics depend on `lifetime_partners`, `partner_age_diff`, and
  rel-type tagging being consistent across variants. Variants that
  drop matches (3, 5, 6) will produce fewer total pairs — diagnostics
  must normalize by total pairs formed, not by total agents, when
  comparing structural shape.

## Approval gate

After Robyn reviews this spec, proceed to writing-plans.
