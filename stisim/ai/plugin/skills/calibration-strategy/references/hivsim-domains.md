# HIVsim / STIsim calibration domains

Reference for the `calibration-strategy` skill. Enumerates the parameter families typically encountered when calibrating an HIVsim / STIsim analysis, with guidance on how each usually classifies against the six-category schema.

**Inspect the current model rather than trusting this list blindly** — HIVsim / STIsim evolves, and parameter names, defaults, and mechanisms change between versions. Cross-reference against [`model-primer/references/calibration-knobs.md`](../../model-primer/references/calibration-knobs.md) for the authoritative index of what STIsim currently exposes.

## HIV prevalence

Usually a **target**, not a knob. Recorded by:

- calendar year (annual or period)
- sex
- age band (5-year usually; sometimes finer for younger ages)
- population group where surveyed separately (FSW / MSM / clients / general population)

Common sources: national household serosurveys (PHIA, SHIMS, DHS AIS), FSW / MSM IBBS studies, UNAIDS Spectrum outputs (with caveat — UNAIDS is itself a model).

## HIV incidence

Usually a **target** where reliable estimates exist. Cohort-based incidence is preferable to biomarker-imputed incidence, which is itself a modeled quantity. Note carefully whether the source is a direct measurement or a Spectrum-style estimate.

## ART

**Observed ART coverage should generally be treated as a strong data constraint**, not freely predicted. HIVsim's `ART` intervention accepts a coverage series and corrects the on-ART population to match each timestep. Feeding the observed series in as a data input is usually correct; calibrating a scalar to reproduce it is usually wrong.

**Diagnosed-pool constraint:** HIVsim can only place *diagnosed* people on ART. If simulated ART coverage falls short of target, the ART coverage series is not the knob — the testing configuration usually is. See main SKILL instruction 5.

## Testing and diagnosis

Often a legitimate calibration domain, because testing rates in the model must be sufficient to feed the ART cascade. Candidate knobs:

- relative testing probability multipliers (by sex, risk group)
- ANC testing rates (for MTCT / pregnant-women pathway)
- symptomatic-testing probability (for STIs with a symptomatic presentation)
- testing frequency by risk group

Do not calibrate testing without checking against observed testing yield or self-reported testing histories if available.

## Sexual network and transmission behavior

The largest genuinely-uncertain calibration domain for HIV / STI models. Candidate knobs typically include (verify current names against the model):

- partnership formation rates
- concurrency (multiple concurrent partners at a given time)
- age mixing preferences
- coital frequency within partnerships
- relationship duration distributions
- population-specific relative-risk multipliers (FSW / MSM / clients)
- per-act transmission probability (rarely — usually a fixed literature parameter, but sometimes opened as a scaling multiplier)
- condom use / effectiveness (usually data-informed if surveys exist)

These are typically the **behavioral / network** category — poorly identifiable from direct evidence, influential on observable epidemiological patterns.

## VMMC

**Coverage is a stock target, not a hazard.** Observed VMMC prevalence by age × year from SHIMS3 / PHIA / DHS should be fed in directly; do not calibrate a scalar to reproduce it. VMMC efficacy (`eff_circ`) is a **fixed literature parameter** anchored to Orange Farm / Rakai / Kisumu RCTs (~0.6 reduction in male acquisition), not a calibration knob.

## PrEP

Historical PrEP coverage where it has been rolled out (relatively recent, limited geographies) is a **data input**. Future PrEP scale-up is an **intervention assumption / scenario input** — never calibrated to historical data.

## Mortality (HIV-specific)

`pars.rel_death` and `pars.rel_death_f` on `sti.HIV` scale CD4-stratified mortality. On-ART mortality has its own age gradient (`art_death_age`) and adherence-specific multipliers. All are candidate knobs, but note the "on-ART never exceeds off-ART" invariant that STIsim enforces at init.

**Also check the background-mortality double-count.** If `ss.Deaths` is being fed all-cause data and `sti.HIV` is also killing agents from AIDS, mortality is double-counted. Use `stisim.data.dedup_deaths` before calibrating any mortality knob — otherwise the calibration will chase a data artifact. This is a **data preprocessing** step, not a calibration step.

## STI-specific parameters

For chlamydia / gonorrhea / trichomoniasis / syphilis / BV, per-disease knobs include:
- transmission probability multipliers
- natural-history durations (acute, latent, recovery)
- symptomatic fraction (usually fixed literature; verify sources)
- treatment efficacy (usually fixed literature)
- congenital / MTCT probabilities by maternal stage (for syphilis)

STI incidence over time is typically the target; treatment / testing coverage feeds in as data.

## Demographics

**Never calibrate demographic inputs.** Fertility, mortality (non-HIV), migration, age structure — all should come from the country's UN or DHS demographic sources via `stisim.data`. If the simulated population trajectory diverges from the observed, the fix is in the data pipeline, not the calibration search space.

## Intervention rollout dates and programmatic quantities

Known historical rollout dates (ART introduction, VMMC scale-up years, PrEP availability) are **data inputs**, never calibrated. If simulated ART cascade differs from observation because a rollout date is wrong in the config, fix the config.

## Behavioral parameters with country-specific evidence

Age at sexual debut, condom use, partnership concurrency reports — where DHS or a national survey has directly measured them for the setting, use the measurement. Do not open these as calibration knobs unless the analysis is *specifically* about uncertainty in the measurement itself.

## Where the current model documentation lives

- **`model-primer/references/calibration-knobs.md`** — authoritative index of what STIsim exposes, updated with the model.
- **`hiv-interventions/references/{art,vmmc,prep,testing}.md`** — per-intervention parameter documentation.
- **STIsim source** — the ultimate ground truth. Any list in these references may lag the code.
