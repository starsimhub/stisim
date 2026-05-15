# PFA comparison: run notes

Branch `feat/pfa-comparison` (off `rc1.5.6`). Re-run 2026-05-15 with sex-aware incidence + prevalence analyzers.

## TL;DR — algorithmic bug in production PFA

**The production matcher (`SortBisect`) and the unrestricted matcher (`SortPair`)
do not honour `age_diff_pars`.** They produce essentially zero average M-F age
gap in MF relationships (~0-3 years instead of the intended +6-8 years). The
strict variants (`DesiredAgeBucket`, `KDTreeNN`, `BandMatch`) match the target;
`GreedyOldEnough` overshoots (~9-10 years). This is most clearly visible in
the calibrated Zimbabwe and Eswatini runs.

### Mean (M-F) age difference at pair formation

**Eswatini (calibrated, n≈10k, 2 reps, 46 yrs, target ~7-8 yrs)**

| Variant            | stable | casual | onetime |
|--------------------|--------|--------|---------|
| GreedyOldEnough    |  7.85  |  9.70  |   9.32  |
| DesiredAgeBucket   |  6.19  |  6.67  |   6.65  |
| KDTreeNN           |  6.13  |  6.54  |   6.14  |
| BandMatch          |  5.73  |  6.19  |   5.83  |
| SortBisect (prod)  |  1.01  |  0.70  |   0.17  |
| SortPair           | -0.31  | -0.13  |  -0.68  |

Zimbabwe shows the same ordering and similar magnitudes.

**Generic block (no demographics/HIV, n=10k, 3 reps, 20 yrs)** — strict
variants converge on ~+7 yrs; SortBisect lands at +2-3 yrs; SortPair lands
near 0 or negative.

### Prevalence (active pairs at sim end) is even more skewed for the strict variants

Long-lived pairs age together so the original gap persists. In Eswatini at
sim end:

| Variant            | stable | casual | onetime |
|--------------------|--------|--------|---------|
| GreedyOldEnough    | 10.06  | 12.37  |  10.40  |
| DesiredAgeBucket   |  8.58  |  8.90  |  10.25  |
| KDTreeNN           |  8.50  |  9.93  |   8.89  |
| BandMatch          |  7.78  |  8.92  |   8.94  |
| SortBisect (prod)  |  1.28  |  0.23  |   0.09  |
| SortPair           |  0.28  | -1.08  |  -1.42  |

### Why `SortBisect` fails on age skew

Both sexes are sorted by their own age (men by `age`, women by `desired_age =
age + age_gap`). The bisect trim only removes head/tail where one group has
no eligible counterpart; it does not enforce the actual age difference. When
the male age distribution and the female desired-age distribution have
similar shapes (which they do, because both come from the same base
population), sort-then-zip pairs by quantile and the average M-F gap
collapses to ~0. The bisect step does not fix this.

`DesiredAgeBucket`, `KDTreeNN`, and `BandMatch` all preserve the gap because
they match on the *actual* desired age (or a discretized proxy) rather than
on sorted index.

## HIV prevalence under each PFA (calibrated runs)

These numbers are from the previous run with the old (sex-agnostic) heatmap
analyzer, but the timings and HIV outcomes are unaffected — only the age
diagnostics needed fixing.

**Zimbabwe (n=10k, 3 reps)**

| Variant            | Final HIV prev. |
|--------------------|-----------------|
| SortPair           | 4.31%           |
| SortBisect (prod)  | 4.17%           |
| GreedyOldEnough    | 4.07%           |
| BandMatch          | 3.80%           |
| DesiredAgeBucket   | 3.68%           |
| KDTreeNN           | 3.64%           |

**Eswatini (n=10k, 2 reps)**

| Variant            | Final HIV prev. |
|--------------------|-----------------|
| SortPair           | 4.18%           |
| SortBisect (prod)  | 3.30%           |
| GreedyOldEnough    | 3.61%           |
| BandMatch          | 3.41%           |
| KDTreeNN           | 2.86%           |
| DesiredAgeBucket   | 2.63%           |

PFA choice moves modelled HIV by ~18% (Zim) to ~59% (Eswatini). The strict
variants produce lower HIV, consistent with the steeper age gap concentrating
incidence among older men → fewer at-risk transmissions.

## Wall time

All variants are within ~10% of each other at the scales tested (n ≤ 10k,
20 yrs generic / 25 yrs Zimbabwe / 46 yrs Eswatini). LSA capped at n ≤ 5k;
otherwise PFA cost is not the deciding factor.

## What changed in the diagnostics

The original `PairAgeHeatmapAnalyzer` was buggy: it read pair keys from
`relationship_durs` (which uses `(min(uid), max(uid))`, **not** sex-ordered)
and looked up ages by raw UID. So plotted `age_p1` and `age_p2` were sex-
agnostic and the resulting heatmaps appeared symmetric regardless of the
real age skew. Fixed by:

- `PairFormationAgesAnalyzer` — hooks `step()`, reads each network's
  `edges.age_p1` (male age at formation) and `edges.age_p2` (female age at
  formation) directly. Records incidence over the full sim.
- `PairPrevalenceAnalyzer` — at `finalize_results`, snapshots active edges
  and records current ages of partnered agents (via `ppl.age[p1]` /
  `ppl.age[p2]`).
- Notebook filters `sw` out of the MF age-mixing panels and renders SW
  pairs separately when present (Zimbabwe / Eswatini).

## Lifetime partner distribution by sex

Means are equal across sexes by construction (each MF edge increments
`lifetime_partners` for one M and one F, populations are balanced). The
differentiation across algorithms shows in the **tails**:

**Generic n=10k, 20 yrs (max lifetime partners by method × sex):**

| Variant            | M mean | M max | F mean | F max |
|--------------------|--------|-------|--------|-------|
| SortBisect (prod)  | 1.11   |  11   | 1.13   |  13   |
| SortPair           | 1.43   |  21   | 1.45   |  15   |
| BandMatch          | 1.10   |  19   | 1.12   |  13   |
| DesiredAgeBucket   | 0.86   |  13   | 0.88   |  10   |
| GreedyOldEnough    | 1.07   |  18   | 1.08   |  12   |
| KDTreeNN           | 1.01   |  16   | 1.03   |  14   |

Observations:

- **DesiredAgeBucket has the tightest tail despite sampling men with
  replacement.** The post-filter on male `concurrency` is binding — duplicate
  draws of the same man get dropped before the match is recorded.
- **SortPair has the heaviest tail on both sexes** (no trim, no replacement,
  but the same low-quantile man gets re-paired across timesteps).
- **All algorithms except DesiredAgeBucket show heavier male tails than
  female.** This is the cross-timestep effect: certain men in the dense
  middle of the age distribution get matched repeatedly while a long tail of
  women with desired ages outside the easy-match band stay unmatched.
- **SortBisect (production) has the shortest male tail** of the algorithms
  that don't post-filter. The bisect trim at each step throws out the
  extremes of the male age distribution, suppressing the "popular middle-aged
  man" effect.

Zimbabwe (with SW) shows max counts in the 100-150 range for both sexes — sex
work dominates the upper tail there.

## Files

- Generic + Zimbabwe results: [pfa_comparison_results.obj](pfa_comparison_results.obj)
- Eswatini results: `/Users/robynstuart/gf/hivsim_eswatini/results/pfa_comparison.obj`
- Notebook: [pfa_comparison.ipynb](pfa_comparison.ipynb) — 8 figures: timing, concurrency, generic-MF incidence + prevalence, Zimbabwe-MF incidence + prevalence, Zimbabwe-SW incidence, lifetime partners by sex (log scale).

## Next steps

- **Decide whether the production PFA should change.** If the model is
  intended to honour `age_diff_pars`, switching to `KDTreeNN` or
  `DesiredAgeBucket` is the natural move. A switch will require recalibrating
  Zimbabwe and Eswatini — HIV prevalence shifts by 15-30%.
- **Sanity-check `age_diff_pars`.** The strict variants produce ~6-8 yr gaps
  even though the user expects 7-8 yrs; that's within tolerance. But it would
  be worth pulling the actual distribution and making sure the loc/scale is
  what we think it is.
- **Onetime pairs sometimes show a small negative age gap** (women slightly
  older) for SortBisect/SortPair. Worth investigating whether that's a
  separate issue with how onetime relationships are classified vs how matches
  are made.
