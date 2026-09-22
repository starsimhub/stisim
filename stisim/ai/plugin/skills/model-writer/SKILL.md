---
name: model-writer
description: Use to author a new STIsim Sim in the user's workspace from a scoped research question — select modules, draft the code, self-critique against known traps, and record the intended model changes in the analysis brief.
---

# STIsim model writer

## When to use

- The user has a scoped analysis question and no existing Sim to build on, or wants to compose a Sim from scratch as a new starting point for an analysis.
- Trigger phrases: "write me a Sim", "let's build the model for this analysis", "compose an HIV / syphilis / STI simulation", "author the Sim", "scaffold the model".
- Typically invoked from `stisim:research-workflow` after the analysis-readiness gate (once it ships), or directly by an experienced user who has already scoped their question.

## When NOT to use

- The user wants to modify an existing calibrated Sim they already have in their workspace. Extend or re-parameterise that Sim rather than authoring from scratch.
- The user hasn't yet scoped the research question, comparison, or setting. Send them to `stisim:getting-started` (once it ships) or ask directly before authoring.
- The user's question needs calibration first — hand off to the `calib:` plugin.

## Instructions

1. **Confirm the analysis intent.** Read `analysis-brief.md` at the workspace root if it exists. If it doesn't exist yet, or the intended-changes section is empty, ask the user for:
   - the research question in one paragraph,
   - the setting (country, subregion, or theoretical / illustrative),
   - what comparison or estimand the Sim needs to support.

   Do not proceed to authoring until this is written down. Update the brief with what the user says.

2. **Identify the study type and its data path.** Consult `references/data-guide.md` and classify the analysis as one of:
   - **Country study** — targets a real country. Needs demographic data, initial prevalence, calibration targets, behavioural inputs.
   - **Sub-national study** — targets a subregion or key population. User brings location-specific data in the expected format.
   - **Theoretical / illustrative** — no country-specific data; small population, short horizon.

   For country and sub-national studies, verify that the required data files are present in the workspace (typically under `data/`). If any are missing, walk the user through obtaining them per `data-guide.md` before drafting the Sim. Do not silently substitute defaults.

3. **Load the model primer.** Invoke `stisim:model-primer` and consult `references/architecture.md` to ground module choices in the composition surface `sti.Sim` actually offers. For a specific area (a disease's natural history, a network's partnership dynamics), read the canonical source per `references/canonical-sources.md` before deciding.

4. **Select modules.** From the analysis intent, decide the module composition:
   - **diseases** — which of the available diseases (HIV, syphilis, chlamydia, gonorrhoea, trichomoniasis, BV, GUD) the analysis actually needs.
   - **networks** — default layered structure, or a narrower composition (MSM-only, FSW-focused) if the question calls for it.
   - **demographics** — for country studies, prefer the `demographics='<country>'` string convention so STIsim's location loaders resolve the data files. For sub-national studies, pass a demographics instance pointing at the user's data.
   - **interventions** — baseline standard-of-care plus any the analysis compares.
   - **connectors** — default composition covers common cross-disease coupling; specify explicit connectors (e.g. `sti.hiv_syph`, `sti.hiv_ng`) when overriding parameters.
   - **analyzers** — whatever produces the outputs the analysis needs.
   - **custom modules** — for non-standard behaviour (e.g. fetal-health tracking), pass via the `custom=[...]` kwarg.

   Note each choice, with a one-line rationale, in `analysis-brief.md` under **Intended model changes**.

5. **Handle model gaps.** If the analysis needs a mechanism or option that STIsim does not yet expose (an age-targeted intervention that lacks age-targeting, an outcome that isn't tracked), invoke `stisim:extend-model` to decide whether the fix is local (in the user's workspace) or upstream (a PR to STIsim). Once `extend-model` returns a disposition, record it in the brief and continue authoring around the resolution.

6. **Draft the Sim.** Propose the target file (a `sim.py` script, a notebook cell, a config module — depending on the user's workspace shape) and get approval before writing. Then write the Sim as clearly as possible:
   - use `sti.Sim(...)` for composition — always, not the base `ss.Sim`;
   - route per-module parameters via the `<slot>_pars` kwargs (e.g. `disease_pars={'hiv': {...}}`);
   - route location-specific demographic tuning via `dem_pars=` (e.g. `dem_pars={'rel_migration': 0.5}`), not by post-hoc mutation of a demographics object;
   - parameterise care-seeking on the relevant disease modules (e.g. `p_symp_care` for the SEIS diseases) or on the testing interventions (e.g. `rel_test` scaling) — that is the practical convention in current STIsim projects, rather than adding a separate `CareSeeking` module;
   - set explicit seeds (three or more) for any stochastic claim the analysis will make;
   - if the analysis has cross-module dependencies (e.g. a custom connector or an intervention that references a disease by name), wire cross-references in `init_pre(sim)`, not `__init__`.

7. **Self-critique.** Run through `references/authoring-checklist.md` against the drafted Sim before returning. Report what checks passed, what failed, and what was flagged for the user to decide. Fix mechanical issues; do not silently paper over structural ones.

8. **Update the brief and hand off.** Confirm `analysis-brief.md`'s **Intended model changes** and **Unresolved decisions** reflect what was authored. Tell the user which file was written, what to run to verify it, and what the next step in `stisim:research-workflow` is.

## How this skill is evaluated

A Sim authored with this skill should (a) run without error on the first try in a fresh conda env with the current stisim installed, (b) produce outputs the analysis question actually needs, and (c) leave `analysis-brief.md` current enough that a cold reader can reconstruct what was built and why. Judged manually across the two example analyses as they mature.

## Checks before completion

1. `analysis-brief.md` updated with module choices and rationales.
2. Any model gaps encountered have a recorded disposition (local, upstream, or deferred) — none silently ignored.
3. Authoring checklist run; passed items noted, failed items either fixed or handed back to the user with a clear question.
4. The user has been told which file was written, how to run it, and what to look at first.
