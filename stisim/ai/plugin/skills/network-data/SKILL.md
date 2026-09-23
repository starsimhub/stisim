---
name: network-data
description: Use to fill STIsim's sexual-network parameters (age of debut, partnership durations, risk-group composition, sex-work parameters) with country-specific data rather than accepting the global defaults. Points at the appropriate data source and workflow for each parameter.
---

# STIsim network data

## When to use

- The user is authoring a Sim for a specific country and reaches the network-parameter step.
- The user is reviewing a Sim they've already authored and wants to check whether the network parameters were fed from local data or left at STIsim defaults.
- The user asks how to find data for a network parameter ("what's the age of sexual debut in [country]?", "where do I get partnership duration data?").

## When NOT to use

- Theoretical or illustrative studies. Defaults are fine — the analysis is about mechanism, not a specific setting.
- Sub-national studies where DHS-style national surveys are the wrong resolution — the user brings their own data.
- Parameters that are calibration knobs rather than data-driven (HIV beta, mortality multipliers). See `model-primer/references/calibration-knobs.md` for the distinction.

## Framing

STIsim's network parameter defaults are best-informed placeholders, not values to accept for any specific country. Age of sexual debut, partnership dynamics, risk-group composition, and sex-work parameters all vary substantially across settings. **For a country study, these should be filled with local data — treating them as calibration knobs is fitting on parameters that have observed answers in the literature.**

If a country's data disagrees with the default, that's the normal case, not a defect in STIsim.

## Instructions

1. **Identify which network parameter needs data.** From the Sim's composition and the analysis brief, list the network parameters that are currently at STIsim defaults. Typical candidates:
   - Age of sexual debut (`debut_f`, `debut_m`) → **see `references/dhs-statcompiler.md`**
   - Coital frequency (`acts`) → DHS behavioral indicators (section to be added)
   - Risk-group composition (`prop_f0`, `prop_m0`, `prop_f2`, `prop_m2`) → DHS partnership counts (section to be added)
   - Concurrency (`concurrency_dist`) → DHS or IBBS (section to be added)
   - Partnership durations (`dur_dist`) → DHS or setting-specific studies (section to be added)
   - Sex-work durations (`dur_sw`, `dur_client`) → IBBS surveys (section to be added)

2. **Consult the appropriate data source reference.** Each network parameter has (or will have) a dedicated reference file under `references/` describing where to find the data and how to translate it into STIsim's parameter form. For age of debut, that's `dhs-statcompiler.md`.

3. **Retrieve the data.** Follow the reference's manual or programmatic workflow. Record the source, the survey year, and the raw value in the analysis brief.

4. **Translate to STIsim parameter form.** Data usually comes as a median or a cumulative distribution; STIsim usually wants a distribution shape (`ss.lognorm_ex(mean, std)`). The reference explains the translation for each indicator.

5. **Update the Sim.** Pass the new values as network parameters (either via `network_pars={...}` in `sti.Sim(...)`, or by constructing the network module directly with the values).

6. **Log the change.** In `analysis-brief.md`, under **Intended model changes**, note which network parameters were filled from data and which remain at STIsim defaults. This makes the calibration boundary explicit — everything at a default is either genuinely unknown or an upstream fix waiting to happen.

## Handoffs

- `stisim:model-primer` — for the architecture of the networks module (`references/architecture.md`) and the parameter catalog (`references/calibration-knobs.md`).
- `stisim:model-writer` — this skill is typically invoked from step 4 of model-writer (module selection / parameter defaults) or as a follow-up when reviewing a drafted Sim.
- `stisim:extend-model` — if a network parameter you need doesn't exist yet in the module, that's a case for extend-model rather than a data hunt.

## How this skill is evaluated

A country-scale Sim authored with this skill in the loop should have its network parameters either (a) filled from a cited data source with survey year recorded in the brief, or (b) explicitly noted as "kept at STIsim default because no country-specific data available." No silent defaults.

## Checks before completion

1. Each data-driven network parameter has a source and survey year in the analysis brief.
2. Any parameter still at STIsim default is explicitly flagged as such (not silently accepted).
3. Where data suggests an STIsim default is wrong for the general case (not just this country), the finding is captured — this is an upstream-fix candidate, not a local-only adjustment.
