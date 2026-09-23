---
name: hiv-interventions
description: Use to design the intervention set for an HIVsim analysis through a structured conversation with the user. Covers ART, HIV testing, VMMC, and PrEP by eliciting data, presenting HIVsim's options, proposing defensible defaults where data is sparse, and confirming choices before writing configuration.
---

# STIsim HIV interventions — design interview

## When to use

- Composing the intervention slot for a new HIV analysis Sim.
- Reviewing an inherited Sim to check whether its intervention set matches the user's data and question.
- Trigger phrases: "design the interventions", "set up ART and testing", "what interventions should I include", "help me configure PrEP / VMMC / testing".

## When NOT to use

- Theoretical or illustrative studies not tied to a real country program — intervention defaults are fine and detailed elicitation is overhead.
- The analysis is scenario-comparison only against a fixed baseline the user already has — this skill is for authoring the baseline, not for designing scenario contrasts.

## Framing

This skill is an **intervention-design interview**, not a form. Walk through interventions section by section. Explain what HIVsim commonly supports, ask what data the user has, propose defensible defaults where data is thin, and confirm choices before generating configuration. Do not silently populate uncertain values; do not treat STIsim intervention defaults as universally appropriate for a specific country.

For each intervention area, distinguish clearly between:
- **Observed program data** the user should feed in (coverage series, testing rates, VMMC counts).
- **Model parameters** governing behavior, uptake, efficacy, or persistence.
- **Assumptions** the user is making because data are unavailable — these must be recorded in `analysis-brief.md` so they remain visible.

Common failure mode this skill prevents: papering over data gaps with defaults, then treating simulated behavior as if it were the model's prediction rather than a consequence of the assumed input.

## Instructions

Proceed through the four intervention areas in this order. For each, follow the same conversational pattern (details in the section-specific instructions below):

1. Explain what HIVsim represents for that intervention area.
2. Ask what data the user has, and in what shape.
3. Explain what HIVsim can support given typical data shapes.
4. Identify gaps between what's needed and what's available.
5. Propose data sources or explicit assumptions for the gaps.
6. Confirm the resulting design with the user before moving on.

Between sections, restate what the user has committed to so the accumulating design stays visible. After all four intervention areas are complete, summarize the whole intervention design back to the user before writing code or updating `analysis-brief.md`.

### Section 1 — ART

Consult `references/art.md` for parameter names, the coverage-matching mechanism, and the diagnosed-pool constraint.

Open with:

> Let's start with ART. What does your ART data look like? Do you have coverage expressed as the **number of people on ART**, or as the **proportion of people living with HIV who are on ART**? HIVsim can support either — and can support a mix over time (counts in one period, proportions in another) if that's what your data looks like. Is it disaggregated by age and sex?

Explain the two modes explicitly, using the language in `references/art.md`:
- **Match observed coverage** (default when `coverage` is passed): HIVsim actively corrects the on-ART population to match the input each timestep.
- **Emergent from process** (`coverage=None`): coverage emerges from ART initiation and discontinuation parameters; the user tunes those to reproduce observed data.

Ask which mode the user wants. Explain that mode 1 is the usual choice — ART coverage is a data constraint, not something the model should freely predict — but mode 2 is available and matters when the user is deliberately studying the drivers of coverage.

Explain the **diagnosed-pool constraint**: HIVsim can only place diagnosed people on ART. Even in mode 1, if the simulated on-ART count is below the target, it may be because there aren't enough diagnosed people. The fix lives in the testing configuration, not in ART parameters. This will come up in the next section.

Confirm the ART design (data file, shape, mode, target period) before moving on.

### Section 2 — HIV testing

Consult `references/testing.md` for the modality menu, the eligibility-function pattern, and the `rel_test` tuning knob.

Open with:

> What HIV testing data do you have? Rates, counts, proportions tested in a period, diagnoses, or something else? Disaggregated by age, sex, population, or testing modality?

Present the standard modality set from the reference analyses:
- **General-population testing** (usually stratified by age, sex, and calendar time)
- **Key-population testing** (FSW is standard; other key populations if the analysis calls for them)
- **Symptomatic / low-CD4 testing** (people with advanced HIV disease are more likely to test)
- **ANC testing** (pregnant women, first-trimester by convention, tested if not already diagnosed)
- **Infant testing** (HIV-exposed infants after delivery)

Ask which of these are relevant. Do not treat the standard set as mandatory — the analysis may not need all of them.

If the user has thin testing data, walk them through the default construction pattern in `references/testing.md`: near-zero testing in the early epidemic, ramping up as the national response scales, plateauing near current coverage. Ask which interpolation shape they want (linear, step changes, custom) rather than picking silently.

Then raise the **testing-to-ART link** explicitly: if the simulated on-ART count comes in below the ART target, first check whether enough people are being diagnosed by the testing intervention. The knobs to adjust are (a) the testing input series itself or (b) `rel_test` on a specific testing modality (default 1.0). Do not make these adjustments automatically without telling the user.

Confirm the testing design (modalities, per-modality data or defaults, tuning approach) before moving on.

### Section 3 — VMMC

Consult `references/vmmc.md` for the coverage-as-stock semantics, the baseline-vs-program split, and where the protective effect lives.

Open with:

> Is VMMC part of the national HIV response in this country and period? Do you have data on circumcision coverage? Reported as program procedures performed, coverage / proportion circumcised, or both? Disaggregated by age? Do you need to distinguish traditional / non-programmatic circumcision from VMMC?

Explain what HIVsim treats as VMMC coverage: a **stock target** (prevalence of circumcised men), not a flow of procedures. Both reference analyses combine traditional baseline and VMMC scale-up into a single prevalence series rather than modeling them separately. The `traditional_prob` parameter exists for cases where the user wants to model traditional circumcision separately, but it's off by default.

Note where the protective effect lives: `eff_circ` is on the HIV disease module, not on the VMMC intervention. Users who want to tune the protective effect for a sensitivity analysis need to adjust that parameter.

Confirm the VMMC design (data file, age stratification, whether traditional is modeled separately) before moving on.

### Section 4 — PrEP

Consult `references/prep.md` for population targeting, the efficacy / adherence / coverage separation, and the divergence between the reference analyses.

Open with:

> Does PrEP need to be represented in this analysis? PrEP was rolled out in most settings from the mid-to-late 2010s onward. If your analysis covers a period before that, and you're not doing forward projections, PrEP may be out of scope.

If PrEP is in scope, ask:
- Which populations receive PrEP in this setting? (FSW is the STIsim default; adolescent girls and young women, serodiscordant couples, general-population risk groups are all possible.)
- What data does the user have? (Number initiating, number currently using, coverage among an eligible population, program targets, historical rollout dates, age/sex-specific coverage, persistence / discontinuation data.)

Explain the HIVsim parameterization. Efficacy (`prep_eff`), adherence (`prep_adh`), and program coverage (`coverage`) are separate parameters that combine multiplicatively — do not collapse them into a single effective-coverage number if the user has separate observations. Course duration (`prep_dur`) governs how long a course lasts before renewal.

**Flag the divergence** between the reference analyses as a cautionary example: `hivsim_zim` calls `sti.Prep()` with no arguments, which activates a default coverage ramp starting in 2004 — well before evidence-based PrEP scale-up in most settings. `hivsim_eswatini` explicitly disabled PrEP in its production calibration for exactly this reason. Do not accept `sti.Prep()` with no arguments without confirming with the user that its implicit defaults match the intended history.

Confirm the PrEP design (in scope or not, populations, data or defaults, which parameters are being fed from data vs assumed) before moving on.

## Data sources

For each intervention, `references/data-sources.md` lists likely sources (UNAIDS, DHS Statcompiler, PHIA, PEPFAR, national program reports, existing HIVsim country analyses) and, where a programmatic path exists, how to retrieve the data. Prefer APIs and machine-readable downloads over manual transcription when they exist.

## Summarizing before writing

Before generating or editing HIVsim configuration, summarize the intervention design back to the user in a compact table or structured description. Suggested shape:

| Intervention | Data source | Shape | Mode / notes | Assumptions to review |
|---|---|---|---|---|
| ART | e.g. UNAIDS AIDSInfo 2015–2024 | Proportion PLHIV on ART, by year | Match-coverage mode | VLS defaults to 100% (no country data) |
| HIV testing | e.g. DHS 2015, national program 2020 | Per-modality annual probability | General + FSW + low-CD4 (no ANC) | Testing 1990–2010 linearly interpolated from zero to 2010 anchor |
| VMMC | e.g. PEPFAR MER 2013–2023 | Coverage by age band | Coverage-as-stock, single series (traditional folded in) | Ages 15–29 targeted, older bands zero |
| PrEP | e.g. Not in scope | — | Not represented | — |

Ask the user to confirm before writing code. Record every assumption in `analysis-brief.md` under **Unresolved decisions** so it stays visible for later review.

## Handoffs

- `stisim:model-primer` — for the architecture of the interventions module and the parameter catalog in `references/calibration-knobs.md`.
- `stisim:network-data` — for country-specific network-parameter data acquisition (age of debut, partnership durations); intervention design assumes the network is already parameterized.
- `stisim:model-writer` — this skill is typically invoked from step 4 (module selection) of model-writer, after diseases and networks are settled.
- `stisim:extend-model` — if an intervention needed doesn't exist in STIsim as-is (an age-targeted intervention that lacks age targeting, a modality not in the current menu), route to extend-model for local-vs-upstream decision.

## How this skill is evaluated

An intervention set designed with this skill in the loop should leave the analysis brief with, per intervention: (a) the data source used, (b) the mode of use (match-coverage vs emergent, or the equivalent for other interventions), (c) any tuning parameters that differ from STIsim defaults, and (d) an explicit list of assumptions for gaps that couldn't be filled from data. A cold reader of the brief should be able to reconstruct why each intervention was configured the way it was.

## Checks before completion

1. Every intervention area (ART, testing, VMMC, PrEP) explicitly discussed — none silently populated.
2. Every data-driven parameter has a source recorded; every default-based parameter is flagged as an assumption.
3. The testing-to-ART link raised at least once (either in ART discussion or testing discussion).
4. PrEP either scoped out with reason or configured with explicit population and parameter choices — no `sti.Prep()` with no arguments.
5. Summary table presented and confirmed by the user before writing configuration.
