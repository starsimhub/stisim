# HIV testing — modalities, construction, tuning

STIsim's testing surface is a few concrete intervention classes
(`sti.HIVTest`, `sti.InfantHIVTest`, `ss.ANCTest`) composed to
represent whatever set of testing modalities an analysis needs.
There is no single "population testing rate" — analyses typically
combine several modalities each targeting a specific population or
condition.

Source: `stisim/interventions/hiv_interventions.py:40` (`HIVTest`),
`hiv_interventions.py:156` (`InfantHIVTest`),
`stisim/interventions/base_interventions.py` (`ANCTest`, multi-disease).

## The standard modality set

Both reference analyses (`hivsim_eswatini`, `hivsim_zim`) use variants
of the following. Not all of them are mandatory; pick what the
analysis question requires.

### General-population testing

One `sti.HIVTest` instance, with `eligibility` set to *undiagnosed
and not on ART, excluding key populations*. Coverage is a per-period
probability, often ramping from near-zero in the early epidemic to
current national HTS coverage.

Reference example (Eswatini and Zimbabwe both use):

```python
other_testing = sti.HIVTest(
    eligibility=lambda sim: ~sim.networks.structuredsexual.fsw
                            & ~sim.diseases.hiv.diagnosed
                            & ~sim.diseases.hiv.on_art,
    test_prob_data=general_pop_ramp,
    years=years,
)
```

### Key-population testing

A separate `HIVTest` instance per key population, each with a
different probability ramp reflecting that population's programmatic
coverage. Both reference analyses model FSW-specific testing at
higher rates than the general population:

```python
fsw_testing = sti.HIVTest(
    eligibility=lambda sim: sim.networks.structuredsexual.fsw
                            & ~sim.diseases.hiv.diagnosed
                            & ~sim.diseases.hiv.on_art,
    test_prob_data=fsw_ramp,   # higher plateau than general
    years=years,
)
```

Other key populations (MSM, PWID, adolescent girls and young women)
follow the same pattern if the model represents them and the analysis
needs them.

### Symptomatic / low-CD4 testing

A `HIVTest` instance whose eligibility function is CD4-gated. Captures
the empirical pattern that people with advanced HIV disease are more
likely to present for testing — either through symptom-driven care
seeking or through provider-initiated testing at clinical
encounters.

Both reference analyses use a CD4<200 threshold:

```python
low_cd4_testing = sti.HIVTest(
    eligibility=lambda sim: (sim.diseases.hiv.cd4 < 200)
                            & ~sim.diseases.hiv.diagnosed,
    test_prob_data=low_cd4_ramp,   # highest plateau of the three
    years=years,
)
```

### ANC testing

A `HIVTest` instance whose eligibility is *pregnant, first-trimester,
undiagnosed*. Eswatini implements ANC testing this way (a `HIVTest`
with pregnancy-scoped eligibility and `dt_scale=False` so the given
probability is per-timestep rather than annualized).

```python
anc_testing = sti.HIVTest(
    eligibility=lambda sim: sim.demographics.pregnancy.tri1_uids[
        ~sim.diseases.hiv.diagnosed[sim.demographics.pregnancy.tri1_uids]
    ],
    test_prob_data=anc_ramp,
    dt_scale=False,   # per-timestep probability, not annualized
    years=years,
)
```

Alternative: `ss.ANCTest` is a multi-disease intervention (HIV +
syphilis + other STIs) that runs at ANC visits. Neither reference
analysis uses it, because they only need HIV testing at ANC — the
`HIVTest` pattern above is simpler and single-purpose. Use `ANCTest`
if the analysis needs coordinated multi-disease testing at ANC.

Zimbabwe (`hivsim_zim`) does **not** currently model ANC testing. If
your analysis includes vertical transmission or maternal / infant
health outcomes, ANC testing is usually essential.

### Infant testing

`sti.InfantHIVTest` — follow-up testing of HIV-exposed infants at
delivery. Scheduled by pregnancy events, not by a rate. Neither
reference analysis currently uses it. Include when the analysis needs
pediatric diagnosis dynamics.

## The `test_prob_data` shape

Usually an array or Series aligned with a `years` axis. Both
reference analyses construct it as concatenated linspace ramps:

- Very low or zero in the pre-testing era.
- Linear ramp from an early scale-up year (~1990 in the reference
  analyses) to a mid-2000s or 2010s target.
- Continued (usually gentler) ramp or plateau to the present.
- Held constant into the projection period.

Reference example:

```python
scaleup_years = np.arange(1990, 2021)
years = np.arange(1990, 2041)
n_years = len(scaleup_years)

fsw_prob = np.concatenate([
    np.linspace(0.00, 0.75, n_years),           # 1990–2020 ramp
    np.linspace(0.75, 0.85, len(years) - n_years),  # 2020–2040 gentle rise
])
```

Linear interpolation is a convention, not a requirement. Alternatives
worth offering the user when their data has structure:

- **Step changes** tied to program introductions or guideline
  changes (e.g. treat-all, universal ART) — replace the linspace with
  piecewise constants at guideline-change years.
- **Custom shape** from a fitted curve through anchor points (spline,
  logistic).

When data are sparse, use the linear-ramp default but label it as an
assumption in `analysis-brief.md` — do not present it as observed
history.

## `dt_scale` — per-timestep vs annualized

`HIVTest.test_prob_data` is interpreted as an annual probability by
default and converted to per-timestep via `ss.probperyear`. Set
`dt_scale=False` to pass the per-timestep probability directly. The
Eswatini ANC testing sets `dt_scale=False` because the input is
already a per-visit probability.

## The `rel_test` tuning knob

`rel_test` is a scalar multiplier on the effective testing probability,
inherited from the `STITest` base class
(`stisim/interventions/base_interventions.py`, default `1.0`).
Multiplies through to every draw.

Neither reference analysis overrides `rel_test`. It is the intended
knob for:

- Reconciling a mismatch between simulated diagnoses and observed
  diagnosis counts when the testing input series is trusted but the
  effective yield doesn't match.
- Sensitivity analyses on testing yield (what if the same input
  produced 20% more / fewer diagnoses in practice?).

Do not tune `rel_test` silently. If used, note the value and the
justification in `analysis-brief.md`.

## Diagnosis-to-ART timing

`dur_dx2tx` on `HIVTest` controls the delay from diagnosis to ART
initiation (default `ss.constant(0)`, meaning ART starts on the next
timestep). Neither reference analysis overrides it. Consider setting
a non-zero delay if the analysis is sensitive to same-day ART policy
or to programmatic delays.

The full mechanism: on diagnosis, `HIVTest` sets `hiv.diagnosed=True`,
`hiv.ti_diagnosed`, samples `dur_dx2tx`, and schedules
`hiv.ti_art = ti + delay`. The ART intervention picks up scheduled
agents on the following timestep.

## Recap for the interview

For each modality the user opts into:
- Confirm the data shape (annual probability, per-visit probability,
  count-based).
- Confirm the population targeting (via eligibility function).
- Confirm the temporal shape (ramp, step, custom) — labeling as
  data-driven or assumption.
- Record the choice in `analysis-brief.md`.
