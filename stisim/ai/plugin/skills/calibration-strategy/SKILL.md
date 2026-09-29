---
name: calibration-strategy
description: Use when scoping the calibration for an HIVsim / STIsim analysis — deciding what should be calibrated, what should be populated from data, what should stay fixed from literature, and what belongs to intervention scenarios rather than to historical fit. Owns the model-specific choices (which parameters, which targets, ART/testing dependency, staged sequencing); delegates the machinery (algorithm, optimizer, prior predictive, re-identification, plotting) to the `calib:*` plugin.
metadata:
  version: "0.1"
  versiondate: "2026-09-23"
---

# HIVsim / STIsim calibration strategy

## When to use

- Scoping the calibration for an HIVsim / STIsim analysis before any optimizer is invoked.
- Diagnosing why a calibration cannot fit a target — before adding parameters or widening bounds.
- Reviewing an inherited calibration script that opens a parameter list not tied to the current research question.
- Trigger phrases: "what should I calibrate", "what are the knobs", "how do I set up this calibration", "the calibration won't fit target X", "which parameters should I open up for this HIVsim model", "how many parameters should I calibrate".

## When NOT to use

- The choice of *calibration algorithm* (Optuna vs. ABC vs. HM vs. MCMC) — that's `calib:method-selection`.
- The *functional form of the pseudo-likelihood* — that's `calib:likelihood-design`.
- Generic identifiability, orthogonality, transformations, or the calibration-vs-model parameter distinction — that's `calib:parameter-engineering`. This skill sits on top and adds HIVsim / STIsim specifics.
- Prior predictive checks, re-identification, model-setup sizing, calibration workflow sequencing — the `calib:` plugin owns all of these.
- The user is asking about generic Bayesian calibration theory — that's not this skill.

## Framing

**A parameter being editable does not make it a calibration parameter.** Many HIVsim / STIsim parameters should be supplied directly from data — age at sexual debut, demographic inputs, observed ART coverage, known intervention rollout dates, empirically-measured biological quantities. Do not tune such quantities merely because changing them improves fit. Generic global defaults in HIVsim / STIsim are **placeholders**, not authoritative values — replacing a default with country-specific evidence is not calibration and is not evidence of an STIsim bug.

**Targets are observations the model should reproduce. Knobs are uncertain model parameters that may be adjusted to reproduce them.** These are two different lists, they should be maintained separately, and every entry on each list needs justification.

**Calibration ≠ tuning until it fits.** A good fit does not imply that fitted parameter values are uniquely identified or empirically estimated. Different parameter combinations may produce similar epidemiological outputs. Preserve that uncertainty; do not over-interpret a fitted value as an estimate of the underlying quantity.

**Division of labor with `calib:` skills.** This skill decides *what* to calibrate against *what* and *why*. The `calib:` plugin decides *how* — algorithm, sampler, likelihood, workflow sequencing, prior predictive, re-identification, diagnostics, plotting. Do not reimplement machinery this skill can delegate. When you find a generic capability missing in `calib:`, prefer fixing it upstream in the plugin rather than working around it here (see `extending-stisim` for the upstream-vs-downstream rule).

## Instructions

1. **Classify every candidate parameter into one of six categories.** Do this before any calibration search space is defined. See [`references/hivsim-domains.md`](references/hivsim-domains.md) for the specific HIVsim / STIsim knobs typically found in each.

    | Category | Should it be calibrated? |
    |---|---|
    | **Directly data-informed** — sufficient empirical evidence exists to set the value (age at debut, demographic inputs, observed ART / VMMC coverage series, rollout dates) | No — feed from data, don't include in search space |
    | **Fixed biological / literature** — external evidence transfers across settings (per-act transmission probabilities from RCTs, VMMC efficacy from Orange Farm / Rakai / Kisumu) | Usually fixed; may vary in sensitivity analysis |
    | **Country-specific but poorly observed** — plausibly differs by setting, direct evidence weak or absent (relative FSW risk multipliers, network mixing structure) | Reasonable calibration candidate |
    | **Behavioral / network** — controls a process not directly identifiable from available data but affects observable patterns (partnership formation rates, concurrency, coital frequency) | Often the important knobs |
    | **Implementation / tuning** — exists to help map mechanisms onto observed outcomes | Legitimate but scientifically less interpretable; declare the interpretation explicitly |
    | **Intervention assumption** — hypothetical future scenario (future PrEP uptake, projected ART coverage) | **Never calibrate to historical data.** Scenario input, not a fit knob |

    Produce an explicit list of what is being calibrated and why. Include a list of what is *not* being calibrated but was considered, with the reason for exclusion — future agents inheriting the calibration should be able to see the pruning decisions.

2. **Identify targets separately from knobs, with full metadata.** For every target record: what quantity, population, sex, age, calendar time, units, source, uncertainty if available, and whether it is a point estimate / count / proportion / rate / prevalence / incidence. For every knob record: parameter name, mechanism controlled, allowed range, justification for varying, whether the range is empirical or pragmatic, likely targets affected, constraints or dependencies. Do not create calibration ranges by taking an arbitrary percentage around the default.

3. **Prefer the smallest defensible calibration space.** Start with the smallest set of parameters that are genuinely uncertain, plausibly differ in the modeled setting, influence the selected targets, and cannot be better supplied directly from evidence. Adding a knob should require justification. An excessively large search space weakens identifiability, allows unrelated parameters to compensate for each other, and can hide model misspecification. Hand off to `calib:parameter-engineering` for the generic identifiability / orthogonality / transformation reasoning.

4. **Diagnose before expanding the calibration space.** If the model cannot fit a target, do not immediately add parameters, widen bounds, increase iterations, or blame the algorithm. First check, in order:
    - **Wrong or missing input data.** Was a quantity that should have been populated from country data left at a global placeholder?
    - **Missing mechanism.** Does the model lack a mechanism required to reproduce the observed pattern? If yes, invoke `extending-stisim` — the mechanism may belong upstream, not as a downstream workaround.
    - **Structural constraint.** Is one modeled process limiting another? See instruction 5.
    - **Incompatible targets.** Are two datasets inconsistent with each other?
    - **Incorrect target definition.** Are the model output and the observed data actually measuring the same quantity? Check denominators, age groups, calendar periods, population definitions, units.
    - **Parameter bound problem.** Only after mechanistic inspection should the agent consider whether a legitimate calibration parameter needs a wider range.

5. **ART / testing standard diagnostic.** HIVsim can only initiate ART among *diagnosed* people. If observed ART coverage cannot be reached, the problem may be inadequate testing rather than the ART initiation parameter. Standard diagnostic path:
    1. Check whether enough PLHIV are diagnosed.
    2. Inspect HIV testing rates and modalities.
    3. Inspect testing tuning parameters (relative testing probability, symptomatic testing, ANC testing).
    4. Check ART eligibility / initiation logic.
    5. Only then conclude ART-related calibration parameters need alteration.

    Do not treat ART and testing calibration as independent.

6. **Stage the calibration where the model structure supports it.** Do not calibrate every process simultaneously by default. A typical staged order for an HIV analysis is demographics / population structure → sexual network + epidemic trajectory → testing + diagnosis → treatment → intervention-era quantities. Inspect the actual model dependencies before prescribing an order; the value of staging is exploiting causal structure, not enforcing a fixed sequence.

7. **Prioritize targets — do not silently assign weights.** Not all observations are equally reliable or important. Targets may differ in weight based on measurement quality, sample size, representativeness, scientific importance, and whether the quantity is directly observed or itself modeled by an external source (UNAIDS estimates, for instance, are model outputs, not observations). Document the weighting strategy and surface consequential choices to the researcher for review.

8. **Keep historical fit separate from projection uncertainty.** Calibration asks *what parameterizations are compatible with historical evidence?* Projection asks *under those plausible parameterizations and future assumptions, what happens next?* Do not tune parameters to obtain a desired projection. Future intervention assumptions are scenario inputs, not calibration targets.

9. **Distinguish calibration from sensitivity analysis.** A parameter may be fixed from data, calibrated, varied in sensitivity analysis, or some combination. State which role each consequential parameter plays. Calibration changes uncertain parameters to reproduce historical quantities; sensitivity analysis asks how conclusions change when uncertain assumptions vary. Confusing them produces "sensitivity analyses" over parameters that were themselves calibrated to the same data.

10. **Inspect model-data agreement visually across every relevant dimension after fitting**, not just as an aggregate loss. Time, age, sex, population group, and every distinct target. Systematic discrepancies can be obscured by a single scalar loss. Delegate the plotting to existing calibration tooling — do not build a private plotter here.

11. **Calibration failure is scientific information.** When calibration repeatedly fails, surface a modeling diagnosis rather than continuing to optimize. Possible outcomes to name: calibration bounds inappropriate; a country-specific input is missing; targets are inconsistent; the model lacks a mechanism; the calibration parameters are insufficiently influential; the chosen parameters are not identifiable; testing constrains ART coverage; a target is incorrectly constructed; the data cannot support the intended level of detail. Do not hide these problems by running more iterations.

12. **Integrate with the research question.** Calibration effort should follow from what the analysis is answering. Age-targeted allocation questions need credible age-specific fits; broad national-outcome questions may not warrant fine-grained calibration in dimensions the decision does not use. Spend calibration effort on the aspects of the model that materially affect the research question. Calibration is not an end in itself.

13. **Hand off to `calib:` for machinery.** Once the specification below is written and the researcher has reviewed consequential choices, hand off:
    - `calib:calibration-workflow` — ten-step sequencing (data → visualization → coverage check → likelihood → workflow → edge cases → re-identification → parameter engineering → small real-data calibration → iterate).
    - `calib:parameter-engineering` — deeper identifiability / orthogonality / transformation work if the parameter list needs pruning or restructuring.
    - `calib:coverage-check` — prior predictive check before real-data calibration.
    - `calib:method-selection` — algorithm choice (Optuna / ABC-SMC / MCMC / HM / SBI).
    - `calib:likelihood-design` — functional form of the pseudo-likelihood.
    - `calib:re-identification` — synthetic-data validation before touching real data.
    - `calib:model-setup` — population size, replicates, run-time.

    Do not reimplement these here. If a `calib:` skill is missing a capability, prefer PR'ing the fix upstream rather than working around it in this skill (`extending-stisim`).

14. **Known `calib:` plugin gaps to work around.**
    - **[calib-plugin#1](https://github.com/InstituteforDiseaseModeling/calib-plugin/issues/1)** — HM checkpoints and other calibration artifacts routinely exceed GitHub's 100 MB push limit. Until the plugin scaffold ships a default `.gitignore` for `outputs/hm_*/`, `*.pkl`, `*.msim`, and similar heavy artifacts, add these to `.gitignore` in the calibration spec's Setup section before the first HM wave runs. Save small salvage artifacts (NROY CSVs, `metrics.json`, `run_config.json`, `log.txt`) into a sibling `nroy/{run_label}/` tree that *is* tracked. When a fix lands upstream, remove the local workaround (`extending-stisim` cleanup rule).
    - **[calib-plugin#7](https://github.com/InstituteforDiseaseModeling/calib-plugin/issues/7)** — no dedicated `closing-a-calibration-project` skill exists yet. When the calibration is complete and the analysis is moving to scenarios, follow the archive-tag + curated `calibration/` directory + cherry-pick-to-main pattern described in that issue. The curated dir should contain `README.md`, `calibration_summary.md`, `methodology.md`, `assumptions.md`, `recalibration_guide.md`, and `artifacts/`. If this becomes routine, PR the skill upstream so future users don't reinvent the pattern.

    Both are surfaced here as reminders, not as workarounds this skill maintains permanently. When either is resolved upstream, delete the corresponding paragraph.

## Output — calibration specification

Before invoking `calib:`, produce a concise calibration specification and place it in durable project memory (per `project-memory`). Ask the researcher to review consequential scientific choices before launching a long calibration.

```
**Research question this calibration serves:** <one sentence>

**Calibration targets** (each with: quantity, population, sex, age, year, units, source, uncertainty, weight if non-uniform)

**Direct data inputs that will NOT be calibrated** (with source)

**Fixed literature parameters** (with source)

**Calibration knobs** (each with: parameter name, mechanism, range, justification, targets affected)

**Structural dependencies** (ART requires diagnosis; VMMC coverage requires prevalence-target semantics; etc.)

**Staging strategy** (if applicable: which processes calibrated in which order, and why)

**Objective / loss choice** (delegated to calib:likelihood-design)

**Algorithm choice** (delegated to calib:method-selection)

**Calibration period / target years**

**Known plugin limitations / issue workarounds** (calib-plugin#1, #7 status)

**Diagnostics to check after fitting** (visual by target and dimension; systematic-discrepancy checks)

**Unresolved concerns** (things the researcher and agent flagged as risky or unclear)
```

The goal is not merely low loss. The goal is a model that reproduces relevant historical evidence for scientifically defensible reasons and is fit for the intended analysis.

## Checks before completion

1. Every parameter in the calibration search space has been classified into one of the six categories, with justification for inclusion.
2. Direct data inputs and fixed literature parameters are documented separately from knobs.
3. Targets and knobs are maintained as separate lists with full metadata.
4. The ART / testing structural dependency has been considered explicitly.
5. Target weighting is documented, not silent.
6. Handoff to `calib:` is explicit — this skill has not reimplemented calibration machinery.
7. The calibration specification is written to durable project memory and reviewed by the researcher before any long run is launched.
8. Known `calib:` plugin gaps (#1, #7) are either worked around or explicitly deferred, with the deferral recorded.
