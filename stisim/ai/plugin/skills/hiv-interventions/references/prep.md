# PrEP — populations, parameters, and a divergence to notice

STIsim's PrEP intervention (`sti.Prep`, source at `stisim/interventions/hiv_interventions.py:739`) parameterizes PrEP across four distinct axes: **efficacy**, **adherence**, **course duration**, and **program coverage**. Do not collapse them into a single number — they combine multiplicatively, and the analysis often has separate observations for each.

## When PrEP is (and is not) in scope

PrEP rolled out in most high-burden settings from the mid-2010s onward — later in many. If the analysis covers only the pre-PrEP period, PrEP is out of scope: leave the intervention out rather than passing zero coverage. If the analysis covers forward projections beyond ~2015, PrEP is usually in scope even if historical coverage was low.

Ask the user directly whether PrEP needs to be represented before launching into parameters.

## Instantiation

```python
prep = sti.Prep(
    prep_eff=0.8,
    prep_adh=1.0,
    prep_dur=ss.months(3),
    coverage=None,          # DataFrame or scalar; default: FSW ramp 2004–2025
    eligibility=None,       # callable(sim) -> uids; default: FSW
)
```

Reference examples:
- **`hivsim_zim`**: `sti.Prep()` with no arguments. This activates the intervention with all defaults, including the implicit `coverage` ramp starting in 2004 and reaching 80% of FSW by 2025 — well before evidence-based PrEP scale-up in most settings.
- **`hivsim_eswatini`**: `sti.Prep()` was explicitly disabled in the production calibration for exactly this reason, retained only for decision-analysis scenarios.

**This divergence is a real caution.** Do not accept `sti.Prep()` with no arguments without confirming that its implicit defaults match the intended history for the analysis. See "Cautions" below.

## The four axes

### `prep_eff` — biological efficacy

Baseline efficacy against acquisition per exposure event, in [0, 1]. Default `0.8`. This is the parameter tied to trial-based estimates (e.g. iPrEx, PROUD) and represents the ceiling under perfect adherence.

### `prep_adh` — adherence multiplier

Multiplier on `prep_eff` representing adherence-adjusted protection, in [0, 1]. Default `1.0` (fully adherent). Real programmatic adherence in many settings is meaningfully below 1.0; the user should supply an empirically justified value if the analysis is sensitive to adherence.

Effective per-exposure protection is `prep_eff * prep_adh` — that's why collapsing them is misleading: you lose the ability to reason about which of the two is driving the result.

### `prep_dur` — course duration before renewal

Default `ss.months(3)`. Governs when an agent's PrEP course expires and needs renewal. Together with adherence, this parameterizes persistence — an agent whose course expires without renewal drops off PrEP.

### `coverage` — program coverage / enrollment target

Same shape as ART and VMMC coverage: DataFrame or scalar. Interpreted as prevalence of PrEP use in the eligible population. The default (when `coverage=None`) is a built-in ramp from 2004 to 2025 reaching 80% of FSW — **not** country-specific and not a historical fact for most settings.

For a country analysis, supply either:

- **A coverage series from national program data** (PEPFAR MER, national HIV program reports, published PrEP scale-up studies).
- **A coverage series constructed from initiation counts** if the data is program-side rather than survey-side — divide initiators by the eligible population size to get a coverage-equivalent.

If no data are available and PrEP must be represented for the analysis to be complete, supply an explicit assumption (starting year, ramp shape, plateau) and record it in `analysis-brief.md`. Do not rely on the implicit default.

### `eligibility` — target population

Callable `(sim) -> uids` that returns the eligible agent IDs each timestep. Default: `sim.networks.structuredsexual.fsw`, i.e. FSW only. To model other populations, provide a custom eligibility function.

Populations worth considering (per national guidelines and programmatic history):
- Female sex workers
- Adolescent girls and young women (AGYW)
- Serodiscordant couples (HIV-negative partner)
- Men who have sex with men (MSM)
- Other high-risk general-population subgroups

Model each population as a separate `sti.Prep` instance with its own `eligibility`, `coverage`, and (if warranted) different `prep_adh` and `prep_dur`. Collapsing multiple populations into a single intervention with a union eligibility function loses the per-population coverage detail.

## Cautions

### The `sti.Prep()` no-args default is not universal

The implicit `coverage` ramp starts in 2004. That predates demonstrable PrEP evidence in most settings and predates program scale-up in every one. `hivsim_eswatini`'s production calibration disabled `sti.Prep()` for exactly this reason.

If the analysis covers 2004 onwards and calls `sti.Prep()` with no arguments, PrEP-attributed protection appears in the results throughout an era in which almost no PrEP was used. This inflates historical PrEP impact.

**Fix**: either pass an explicit `coverage` argument (a DataFrame whose early years are zero, ramping up from the country's actual rollout year), or leave the intervention out entirely for the pre-PrEP period.

### Data source specificity

PrEP is program-specific in a way ART and testing often aren't. There is rarely a single "national PrEP coverage" number; coverage is typically reported per program, per target population, per funding stream (PEPFAR-supported vs government, etc.). When compiling coverage data, the user should be explicit about which programs and populations the number covers, and reflect that in the `eligibility` function.

## Recap for the interview

- In scope for this analysis's time window?
- If yes: which populations?
- For each population: coverage data, efficacy, adherence, persistence — separate values or a single collapsed number?
- If the user proposes `sti.Prep()` with no arguments, push back: what's the implicit default doing to the results, and is that intended?
