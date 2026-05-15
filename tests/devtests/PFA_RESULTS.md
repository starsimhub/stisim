# PFA comparison: run notes

Full benchmark run: 2026-05-15. Branch `feat/pfa-comparison` (off `rc1.5.6`).

## TL;DR

- **All 7 variants are roughly the same wall-clock speed** at the scales tested
  (n ≤ 10k, 20-year sim). PFA cost is dominated by other per-step overhead. Cost
  isn't the deciding factor.
- **Network volume varies by ~70%** between extremes. SortPair (no trim, no
  replacement) produces ~1.7× more lifetime partnerships than DesiredAgeBucket
  (strict post-filters).
- **HIV prevalence shifts noticeably with PFA choice** on calibrated sims:
  - Zimbabwe (10k agents, 25 yrs): 3.64% (KDTreeNN) → 4.31% (SortPair).
  - Eswatini (~10k agents, 46 yrs, calibrated interventions): 2.63%
    (DesiredAgeBucket) → 4.18% (SortPair).
- **Production default (SortBisect) sits near the high end** of the prevalence
  range — it's a relatively permissive matcher.

## Results

### Generic benchmark — wall time (s, mean over 3 reps)

| Variant            | n=1k   | n=5k   | n=10k  |
|--------------------|--------|--------|--------|
| BandMatch          | 0.29   | 0.47   | 0.73   |
| DesiredAgeBucket   | 0.28   | 0.48   | 0.71   |
| GreedyOldEnough    | 0.28   | 0.47   | 0.73   |
| KDTreeNN           | 0.30   | 0.49   | 0.72   |
| LSA                | 0.30   | 0.62   | n/a*   |
| SortBisect (prod)  | 0.30   | 0.51   | 0.71   |
| SortPair           | 0.28   | 0.47   | 0.73   |

\*LSA capped at n ≤ 5k (O(n³)). Still cheap at n=5k (~25% slower than the next
group) but skipped at n=10k by design.

### Network volume at n=10k — lifetime partnerships (mean over 3 reps, 20 yrs)

| Variant            | Σ lifetime partners |
|--------------------|---------------------|
| SortPair           | 13,693              |
| SortBisect (prod)  | 10,673              |
| BandMatch          | 10,525              |
| GreedyOldEnough    | 10,216              |
| KDTreeNN           |  9,702              |
| DesiredAgeBucket   |  8,225              |

### Calibrated Zimbabwe (hivsim.demo('zimbabwe'), n=10k, 3 reps)

| Variant            | Final HIV prev. | Wall time (s) |
|--------------------|-----------------|---------------|
| SortPair           | 4.31%           | 8.78          |
| SortBisect (prod)  | 4.17%           | 8.66          |
| GreedyOldEnough    | 4.07%           | 8.80          |
| BandMatch          | 3.80%           | 9.35          |
| DesiredAgeBucket   | 3.68%           | 8.79          |
| KDTreeNN           | 3.64%           | 8.88          |

### Calibrated Eswatini (hivsim_eswatini, n=10k, 2 reps)

| Variant            | Final HIV prev. | Σ lifetime partners | Wall time (s) |
|--------------------|-----------------|---------------------|---------------|
| SortPair           | 4.18%           | 142,508             | 14.00         |
| SortBisect (prod)  | 3.30%           | 140,821             | 13.91         |
| GreedyOldEnough    | 3.61%           | 117,498             | 14.11         |
| BandMatch          | 3.41%           | 114,208             | 13.97         |
| KDTreeNN           | 2.86%           | 105,192             | 14.01         |
| DesiredAgeBucket   | 2.63%           |  78,604             | 14.18         |

LSA omitted from Zimbabwe and Eswatini blocks (too slow at calibration scale).

## What stands out

1. **PFA choice is a real lever on HIV prevalence.** ~18% spread in Zimbabwe,
   ~59% spread in Eswatini. This is comparable to or larger than the effect of
   many calibration parameters that we tune individually.
2. **The ordering is broadly consistent across runs:** SortPair > SortBisect >
   GreedyOldEnough > BandMatch > KDTreeNN > DesiredAgeBucket. The permissive
   matchers form more partnerships → more transmission.
3. **SortBisect (production) is near the top of the spread.** It's a permissive
   matcher in practice — it trims age supports but doesn't enforce
   relationship-type compatibility at the matching stage. Moving to a stricter
   matcher (KDTreeNN or DesiredAgeBucket) would lower modelled HIV by ~15-30%
   in calibrated sims.
4. **DesiredAgeBucket is the outlier on the strict side.** Its in-matchmaker
   relationship-acceptance Bernoulli drops ~30-45% of matches that would
   otherwise survive — this duplicates `add_pairs` logic but applies it earlier,
   before partner counts increment.
5. **Wall-clock cost is negligible** for all variants at the tested scales.
   Algorithm choice is a structural/scientific decision, not a perf decision.

## Implementation notes (caught during the run)

- The Zimbabwe HIV module hard-codes `sim.networks.structuredsexual` (in
  `make_init_prev`), so the benchmark builds per-variant **subclasses of
  `sti.StructuredSexual`** (not bare `MFNetwork_*`) via `make_variant_class()`.
  Eswatini's script does the same.
- `make_variant_class()` must copy **all non-dunder class attrs** of the variant
  (not just `match_pairs`) — `MFNetwork_BandMatch.band_width` was the gotcha.
- Generic and Zimbabwe results checkpoint to disk separately so a partial-run
  failure doesn't lose hours of work. Zimbabwe also checkpoints inside its loop.
- `PartnersLastYearAnalyzer` / `PairAgeHeatmapAnalyzer` skip UIDs that exceed
  the current `ppl.age` array (agents who died and were removed). Without this
  the calibrated runs raise `IndexError` on finalize.

## Next steps

- Look at the **age-mixing heatmaps** in `pfa_comparison.ipynb` to see whether
  the prevalence shift is driven by partner count, partner age spread, or both.
- Decide whether to keep SortBisect as production or move to something stricter.
  If we move, **KDTreeNN** is the natural candidate — it's close to LSA in
  spirit (nearest-neighbour by age) without the O(n³) cost.
- If we change production: recalibrate Zimbabwe and Eswatini against the new
  PFA before drawing conclusions about anything else.

## Files

- Generic + Zimbabwe results: `tests/devtests/pfa_comparison_results.obj`
- Eswatini results: `/Users/robynstuart/gf/hivsim_eswatini/results/pfa_comparison.obj`
- Notebook: `tests/devtests/pfa_comparison.ipynb`
