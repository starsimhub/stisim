# ART — the two modes

STIsim's ART intervention (`sti.ART`, source at `stisim/interventions/hiv_interventions.py`) supports two distinct modes of use. Choosing between them is a design decision that changes what the model can tell you.

## Mode 1: match observed coverage (the usual choice)

```python
art = sti.ART(coverage=art_data)   # coverage is a DataFrame or scalar
```

When a `coverage` argument is supplied, `sti.ART` treats it as a constraint. On every timestep it calls its internal `art_coverage_correction()`, which adds or removes agents from ART to exactly match the input each timestep (per stratum if the data is stratified). Both reference analyses (`hivsim_eswatini`, `hivsim_zim`) use this mode.

This is the right choice when observed ART coverage is a data input the model should reproduce rather than something the model should predict — which is the standard framing.

### Supported data shapes

The `coverage` argument accepts:

- **A scalar** — a constant coverage over all time, sexes, and ages. Rarely useful for real analyses.
- **A DataFrame with an `n_art` column** — absolute number on ART, by the strata present in the file (Year, Gender, AgeBin as available).
- **A DataFrame with a `p_art` column** — proportion of PLHIV on ART, by the same strata.
- **A mix over time** — a file whose format changes between periods (e.g. counts pre-2000, proportions after) is supported when the loader is set up for it.

Stratification is optional. Files can carry Year alone (aggregate national), Year+Gender, Year+Gender+AgeBin, etc. STIsim's coverage correction runs per stratum present in the data.

Reference examples:
- Eswatini: `art = sti.ART(coverage=art_data, vls_coverage=art_vls_coverage)` where `art_data` is a stratified Year × Gender × AgeBin table.
- Zimbabwe: `art = sti.ART(coverage=p_art)` where `p_art` is aggregate annual proportions.

## Mode 2: emergent from process (`coverage=None`)

```python
art = sti.ART()   # no coverage argument, or coverage=None explicitly
```

When `coverage` is not supplied, ART is unconstrained. Whether an agent starts (and continues) ART emerges from:

- Whether they've been diagnosed (see `references/testing.md`).
- The `art_initiation` parameter (default `ss.bernoulli(p=0.9)` — probability a newly diagnosed agent starts).
- Any discontinuation dynamics represented in the current STIsim ART implementation.

In this mode, the user reproduces observed ART coverage by tuning the underlying processes — testing rates, initiation probability, any discontinuation parameters — rather than by setting a coverage target. This matters when the analysis is *about* the drivers of ART coverage (what would happen to coverage if testing scaled up differently, for example) rather than assuming the coverage trajectory as an input.

## Diagnosed-pool constraint (bites in both modes)

The ART intervention can only place *diagnosed* people on ART. In Mode 1, this means the coverage correction cannot fill above the diagnosed pool. If the simulated on-ART count runs below the input target even though the correction is active, the constraint is usually upstream: not enough diagnosed people.

**Diagnostic sequence** when this happens:

1. Check the testing configuration. Is there a modality serving the population that's short of ART? Is its input probability sensible?
2. If the testing input looks right but too few people are being diagnosed, check `rel_test` on the testing modality (default 1.0; see `references/testing.md`) — it multiplies the effective testing probability.
3. Only after (1) and (2) look sound is it worth suspecting the ART parameters.

Do not tune ART parameters to compensate for a diagnosis deficit. That papers over the actual mismatch.

## Other ART parameters

- **`vls_coverage`** — viral load suppression coverage among those on ART. Same shape as `coverage` (DataFrame or scalar). Default is `None`, treated as 100% VLS. Eswatini fits this from SHIMS surveys (2016, 2021) and back-fills to 1985; Zimbabwe uses the default.
- **`art_initiation`** — probability a newly diagnosed agent starts ART, default `ss.bernoulli(p=0.9)`. Neither reference project overrides this.
- **`pmtct_efficacy`** — efficacy of maternal ART in preventing vertical transmission, default 0.96. Neither reference project overrides.

For the full parameter list see `stisim/interventions/hiv_interventions.py:179`.

## Which mode to recommend

Default: Mode 1 (match coverage), unless the analysis is specifically about the drivers of coverage. Say this to the user rather than picking silently. If Mode 2 is chosen, note explicitly in `analysis-brief.md` that ART coverage is emergent and record which process-level knobs the user intends to tune.
