---
name: result-extraction
description: Use when writing analysis, plotting, or post-processing code that reads sim results. Directs to ss.Result / ss.Results built-in methods (annualize, resample, to_df) instead of hand-rolled groupby aggregation, which silently mishandles the flow-vs-stock distinction.
---

# STIsim result extraction

## When to use

- Any code that reads `sim.results`, a `MultiSim`, or a parquet / csv exported from a sim.
- Aggregating from timestep to annual, monthly, or other cadences.
- Building analyzer outputs, plotting series, or writing summary tables.
- Trigger phrases: "aggregate the results", "annual rates", "roll up to year", "compute per-year", "group by year", "the plot needs annual data", "convert to yearly".

## When NOT to use

- The user is inspecting a single-timepoint value (e.g. `sim.results.hiv.prevalence[-1]`) — no aggregation to do.
- The user is computing a genuinely custom quantity (e.g. incidence with a bespoke denominator) that no `ss.Result` method covers. Use step 4 below rather than hand-rolling from scratch.

## Framing

Every `ss.Result` carries a `summarize_by` field, or falls back to a name-based heuristic (`new_* → sum`, `n_* → mean`, `cum_* → last`, see [`starsim/results.py:149`](../../../../../../starsim/starsim/results.py)). The built-in API — `annualize()`, `resample()`, `to_df(resample=…)` — reads that field and applies the correct aggregation automatically. Hand-rolled `df.groupby('year').mean()` or `.sum()` bypasses this, so the analyst has to re-decide flow-vs-stock at every call site and gets it wrong silently when they guess.

The wrong answer is usually a factor of `n_timesteps_per_year` off (flows aggregated as means, or stocks aggregated as sums). At weekly or monthly timesteps that is a 12–52× error, but it still looks plausible for calibration ranges — the errors don't scream, and they compound with any downstream per-capita or per-year calculation that uses the wrong quantity as its numerator.

Common failure modes this skill prevents:
- `df.groupby('year')['new_infections'].mean()` — silently divides annual incidence by the number of timesteps per year.
- `df.groupby('year')['n_alive'].sum()` — silently multiplies population by the number of timesteps per year.
- Reading a saved parquet, losing the `ss.Result` wrapper, then guessing the aggregation method from column names alone.

## Instructions

1. **Identify what the input is.** If it is `sim.results`, an `ss.Result`, or an `ss.Results` container, use the API directly — do not convert to a raw DataFrame first. If it is a saved parquet / csv that has lost the Result wrapper, note the round-trip loss and prefer reloading the sim or MultiSim if it is still available.

2. **Choose the method that matches the goal:**
   - `result.annualize()` — annual rollup, fast numpy-based, respects `summarize_by` / `summary_method()`. Use this for the common case.
   - `sim.results.annualize()` (container-level) — annualize every result in the container at once with the correct method for each. Use when producing a whole-sim summary table.
   - `result.resample(new_unit='month', …)` — non-annual cadences, pandas-backed, more flexible.
   - `result.to_df(resample='year')` — annual DataFrame ready for plotting or joining, with low/high bounds preserved.

3. **If the result name doesn't match the heuristic, pass `summarize_by=` explicitly.** The heuristic only recognises `new_*`, `n_*` / `_n_*`, and `cum_*`. A result named `treatment_success_rate` or `prev_15_49` will default to `mean`, which may or may not be what you want. Either set `summarize_by='sum'|'mean'|'last'` on the Result at creation time (preferred — the decision travels with the data), or pass it into `annualize()` / `resample()` at the call site.

4. **Never hand-roll aggregation over sim results.** If the code you are about to write contains `df.groupby(…)[<result_col>].mean()` or `.sum()`, stop and use one of the above instead. The only legitimate exception is a genuinely custom denominator that no method covers — and even then, aggregate the numerator and denominator separately using the correct per-Result method, then divide. Do not compute per-capita rates by summing a stock.

5. **When saving results for downstream analysis or dashboarding**, save the annualized (or otherwise-resampled) table rather than the raw timestep table. Downstream code cannot recover the flow-vs-stock decision from column names alone once the Result wrapper is lost. If raw timestep data must be preserved, save both, and label the resampled file so downstream code reaches for it by default.

## Checks before completion

1. No `df.groupby('year').mean()` or `.sum()` over `new_*`, `n_*`, or `cum_*` columns.
2. Any explicit aggregation of a non-standard-named result passes `summarize_by=` rather than relying on the fallback default.
3. If the code exports for downstream use, the exported table is already annualized (or otherwise resampled) — downstream doesn't have to re-derive it.

## How this skill is evaluated

Analysis and plotting code produced under this skill should be readable without the analyst having to re-derive the flow-vs-stock decision at each call site. On review, `git grep "groupby.*year"` in analysis code should return zero hits over sim-result columns.
