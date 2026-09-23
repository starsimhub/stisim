# Authoring checklist

Run this checklist against a Sim drafted by `stisim:model-writer` before
handing back to the user. Each item is a specific question with a
concrete pass/fail signal. Skip nothing silently; if an item does not
apply, note that it was considered.

## Composition

- [ ] **Module selection matches the analysis.** Every module in the
  Sim traces back to a specific need from the research question. Every
  need from the question is served by some module or explicitly noted
  as deferred in `analysis-brief.md`.
- [ ] **No unused module.** Any module in the composition that isn't
  needed by the analysis is either removed or its inclusion is
  justified (e.g. "kept to preserve the default calibration"). Extra
  modules cost run time and complicate interpretation.
- [ ] **Care-seeking parameterized.** If the Sim has any care-based
  intervention (testing, treatment, ANC screening), care-seeking
  intensity is set on the relevant disease modules (`p_symp_care`
  on the SEIS diseases) or on the testing interventions (`rel_test`
  and similar scaling), matching the current STIsim convention.
  Adding the standalone `CareSeeking` module is not the pattern.

## Data pre-processing

- [ ] **Deaths file deduplicated for HIV-burden countries.** For any
  country study with HIV in the composition and any AIDS-era history,
  the deaths CSV must be run through `stisim.data.dedup_deaths` before
  the sim reads it. UN WPP all-cause rates include AIDS deaths; the
  HIV module also kills agents; feeding raw rates double-counts. Skip
  only for very-low-prevalence settings where the AIDS hump is
  negligible.

## Parameters

- [ ] **Beta values set.** For each disease, transmission-relevant
  parameters (`beta_m2f`, `rel_beta_f2m`, etc.) are explicitly set
  rather than left as placeholder zeros. If defaults are intentionally
  used, note the choice in the brief.
- [ ] **Location parameters routed via `dem_pars`.** If the Sim uses
  country-specific demographics (`demographics='<location>'`),
  overrides to migration, pregnancy, or death rates are passed via
  the `dem_pars=` kwarg — not via post-hoc mutation of a demographics
  object that hasn't been instantiated yet.
- [ ] **Per-module pars via `<slot>_pars`.** Disease, network,
  intervention parameters travel via `disease_pars={...}`,
  `network_pars={...}`, `intervention_pars={...}` respectively, keyed
  by module name. No parameters hidden in the top-level `pars` dict
  that should belong to a specific module.

## Wiring and lifecycle

- [ ] **Cross-references in `init_pre(sim)`, not `__init__`.** Any
  custom module (connector, intervention, analyzer) that references
  other modules by name resolves those references in `init_pre(sim)`.
  Cross-refs bound in `__init__` refer to the pre-clone object and
  will not point at the running sim's modules.
- [ ] **Starsim `Arr` API used correctly.** Any custom code reading
  state variables uses `.values`, `.isnan`, `.notnan`, `.notnanvals`.
  No manual `arr[sim.people.auids]` indexing. No `arr.raw` inspection
  for logic (only for debugging).

## Time and stochasticity

- [ ] **Timestep appropriate.** Default is monthly. If the analysis
  depends on a process narrower than one month (specific-week ANC
  screening, short latent phase), the timestep is either reduced or
  the process is re-cast as a hazard rate rather than a scheduled
  event.
- [ ] **Seeds specified.** At least three distinct random seeds are
  set for any stochastic output the analysis will report on. A
  single-seed run is a debugging artifact, not a result.
- [ ] **Event scheduling is integer-safe.** If custom code schedules
  events at `ti + fractional_duration`, the target is floored (or
  otherwise coerced to an integer) before equality-comparing to
  `self.ti`. Floating scheduled times silently miss firing.

## Results and outputs

- [ ] **Analyzers present for every analysis output.** Every quantity
  the analysis needs to compare, plot, or report has a corresponding
  analyzer or a standard `ss.Result` that captures it.
- [ ] **Annual resampling via `.annualize()`.** Any code that produces
  annual summaries from a `ss.Result` uses `.annualize()` or
  `.to_df(resample='year')` — flow vs stock semantics are already
  handled. No manual `groupby('year').mean()` on flow results.
- [ ] **No absolute claims where the model has a structural ceiling.**
  If the analysis will report on a quantity known to have a structural
  ceiling in STIsim (syphilis absolute prevalence is the standing
  example), that limitation is flagged in `analysis-brief.md` under
  **Unresolved decisions** or **Structural limitations**, and the
  analysis frames its claim as a relative contrast rather than an
  absolute number.

## Handoff

- [ ] **`analysis-brief.md` updated.** Module choices and rationales
  are in **Intended model changes**. Any model gaps have a recorded
  disposition (local, upstream, or deferred) in **Unresolved
  decisions** or **Model extensions**.
- [ ] **User instructions given.** The user has been told the target
  file, how to run it, and what to look at first (an expected result,
  a plot, an assertion).
