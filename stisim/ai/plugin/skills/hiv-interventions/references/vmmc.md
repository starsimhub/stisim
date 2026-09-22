# VMMC — coverage as a stock target

STIsim's VMMC intervention (`sti.VMMC`, source at
`stisim/interventions/hiv_interventions.py:579`) treats male
circumcision as a **stock target**: the input series is prevalence of
circumcised men, and the intervention adjusts each timestep to hit
that prevalence. It is *not* a flow-of-procedures input.

This shape matters for how you feed data.

## Instantiation

```python
vmmc = sti.VMMC(coverage=vmmc_data)
```

Reference examples:
- Eswatini: `vmmc = sti.VMMC(coverage=vmmc_data)` where `vmmc_data`
  is a stratified Year × Gender × AgeBin table of `p_vmmc`.
- Zimbabwe: `vmmc = sti.VMMC(coverage=n_vmmc)` — aggregate counts by
  year.

## Supported data shapes

Same as ART: the coverage parser (`parse_coverage`) accepts:

- **A scalar** — constant prevalence.
- **A DataFrame with `p_vmmc`** — proportion of men circumcised, by
  the present strata.
- **A DataFrame with `n_vmmc`** — absolute count of circumcised men.

Stratification is optional. Age bins (`AgeBin`) are the most
consequential stratification because circumcision prevalence rises
sharply with age in scaled-up programmes, and coverage is often
tracked in 5-year bins.

## Baseline (traditional) vs programme scale-up

Both reference analyses **combine traditional baseline and programme
VMMC in a single prevalence series** — the input CSV runs across the
full period and the intervention treats whatever prevalence the data
shows as the joint result of both. Eswatini's data starts at 1990
with low prevalence (~0.25–2%), rising to survey-measured levels by
2007+ (4–19%, age-dependent). Zimbabwe's data starts from the
programme launch in 2008.

STIsim exposes a `traditional_prob` parameter (default
`ss.bernoulli(p=0)`, i.e. off) for representing traditional
circumcision as a separate always-on process at a given age
(`traditional_age`, default 15). Neither reference analysis uses it.

Use `traditional_prob` when the analysis needs to distinguish
traditional and programmatic circumcision explicitly — e.g. an
analysis of VMMC programme impact that must not conflate the two.
Otherwise, feeding a combined prevalence series is simpler.

## Age targeting

Expressed via the coverage DataFrame's `AgeBin` column. Eswatini
provides bins [10,15), [15,20), ..., [60,65); the per-stratum
correction hits each bin against its target. Zimbabwe's unaggregated
data applies a single national prevalence to all ages.

For programmes that explicitly targeted a specific age range (VMMC
scale-up commonly targets 10–29 or 15–34), fill zeros in the bins
outside the target range if the observed data supports it. Do not
extrapolate from the target range to older bands unless the data
does.

## Protective effect

**The protective multiplier on male susceptibility lives on the HIV
disease module, not on the VMMC intervention.** The parameter is
`eff_circ` and is applied in `HIV.circumcise()`. VMMC hands off to
this parameter for the actual effect on transmission.

This split matters for two reasons:

1. **Sensitivity analyses on protective effect** need to tune
   `sti.HIV`'s `eff_circ`, not the VMMC intervention.
2. **The intervention itself has no efficacy parameter** — passing
   `coverage=data` to `VMMC` just moves agents into the circumcised
   state; the risk reduction they experience is determined by the
   HIV module's `eff_circ`.

Neither reference analysis overrides `eff_circ`; both use stisim's
default. If the analysis question involves the protective effect
directly, tune it explicitly and record the choice.

## Recap for the interview

- Is VMMC in scope (country + period)?
- What data is available (procedures counted, coverage measured, or
  both)?
- Age-disaggregated?
- Traditional circumcision — combine into a single series (default)
  or model separately (`traditional_prob`)?
- Is the analysis sensitive to the protective effect? If so, plan a
  sensitivity on `eff_circ`.

For non-VMMC settings (or the pre-programme period), leave the
intervention out entirely rather than passing zero coverage — the
model treats no intervention differently from zero coverage in some
edge cases.
