# PFA comparison: run notes

Status: implementation complete on `feat/pfa-comparison`; full benchmark not yet run.

## What's here

- Seven `MFNetwork` variants in [stisim/pfa_variants.py](../../stisim/pfa_variants.py):
  `LSA`, `SortPair`, `DesiredAgeBucket`, `SortBisect` (production), `GreedyOldEnough`, `KDTreeNN`, `BandMatch`.
- Smoke tests in [tests/test_pfa_variants.py](../test_pfa_variants.py) — 9 passing.
- Diagnostic analyzers in [tests/devtests/pfa_diagnostics.py](pfa_diagnostics.py).
- Generic + Zimbabwe benchmark in [devtest_pfa_comparison.py](devtest_pfa_comparison.py)
  (`--quick`, `--skip-zimbabwe` flags).
- Notebook builder + notebook in [_build_pfa_comparison_notebook.py](_build_pfa_comparison_notebook.py)
  and [pfa_comparison.ipynb](pfa_comparison.ipynb).
- Calibrated Eswatini comparison in
  `/Users/robynstuart/gf/hivsim_eswatini/run_pfa_comparison.py` (branch `update-coverage`).

## --quick mode observations (n=1k, 1 rep, 5 years)

All seven variants run without error. Wall times are dominated by per-step overhead at
n=1k so the algorithmic differences aren't yet visible — the full n=10k × 20yr run is
needed to see real separation.

Edge counts at sim end (`lifetime_partners.sum`):

| Variant            | total partnerships |
|--------------------|--------------------|
| SortPair           | 974                |
| LSA                | 984                |
| SortBisect (prod)  | 750                |
| GreedyOldEnough    | 720                |
| BandMatch          | 708                |
| KDTreeNN           | 644                |
| DesiredAgeBucket   | 533                |

DesiredAgeBucket is markedly lower because its post-filter on relationship-acceptance
Bernoulli drops a substantial fraction of matches inside `match_pairs`. SortPair and LSA
are highest because they don't trim either side of the age support.

Eswatini smoke run (SortBisect, 1 rep, full timeline): 13 s; 12,562 partner records
and 146,638 lifetime partnerships across all agents — diagnostics fire correctly under
demographics.

## To run the full benchmark

```bash
cd /Users/robynstuart/gf/stisim/tests/devtests
python devtest_pfa_comparison.py             # n=1k/5k/10k × 3 reps + Zimbabwe
# (skip Zimbabwe if you want network-only timings)
python devtest_pfa_comparison.py --skip-zimbabwe
```

Then:

```bash
jupyter nbconvert --to notebook --execute pfa_comparison.ipynb --output pfa_comparison_executed.ipynb
```

For Eswatini (n=10k calibrated):

```bash
cd /Users/robynstuart/gf/hivsim_eswatini
python run_pfa_comparison.py    # ~13s × 6 variants × 2 reps ≈ 3 min
```

## Known gotchas

- `diseases='sis'` doesn't work with `MFNetwork` because `net_beta` reads
  `disease.pars.eff_condom`, undefined on `ss.SIS`. The generic benchmark uses
  `diseases=[]`.
- `PartnersLastYearAnalyzer` and `PairAgeHeatmapAnalyzer` skip UIDs that exceed the
  current `ppl.age` array (agents who died and were removed in a demographics-enabled
  sim). Without this filter, Eswatini's birth/death dynamics raise IndexError on
  finalize.
- LSA is capped at n ≤ 5,000 in the generic block; omitted entirely for Zimbabwe and
  Eswatini (O(n³) cost).
- Eswatini's `run_pfa_comparison.py` swaps in a per-variant subclass of
  `sti.StructuredSexual` rather than a bare `MFNetwork_*` — preserves SW state, risk
  groups, condom data, and HIV result paths.

## Open questions for the full run

- Where does each variant cross from "fast enough" to "unusable"?
- How much do figures 2-3 (partners-last-year, age heatmaps) shift between SortBisect
  and the more permissive variants (SortPair, LSA)?
- Does Eswatini's calibrated HIV prevalence move under different PFAs, or do post-
  pair-formation dynamics (acts, condoms, intervention coverage) wash out the
  algorithmic differences?
