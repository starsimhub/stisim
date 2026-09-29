# Canonical sources for STIsim

Where to find authoritative information for each part of the model. If a task requires deep semantics, read directly from here rather than paraphrasing the primer.

## Core

| Concept | Source | Docs |
|---|---|---|
| Sim | `stisim/sim.py` | `docs/user_guide/index.md` |
| Parameters | `stisim/parameters.py` | — |
| Utils | `stisim/utils.py` | — |
| Data files | `stisim/data/` | — |

## Diseases

Source lives in `stisim/diseases/`; each disease has its own file.

| Disease | Source | User guide |
|---|---|---|
| Base STI | `stisim/diseases/sti.py` | `docs/user_guide/diseases/` |
| HIV | `stisim/diseases/hiv.py` | `docs/user_guide/diseases/` |
| Syphilis | `stisim/diseases/syphilis.py` | `docs/user_guide/diseases/` |
| Gonorrhea | `stisim/diseases/gonorrhea.py` | `docs/user_guide/diseases/` |
| Chlamydia | `stisim/diseases/chlamydia.py` | `docs/user_guide/diseases/` |
| Trichomoniasis | `stisim/diseases/trichomoniasis.py` | `docs/user_guide/diseases/` |
| BV | `stisim/diseases/bv.py` | `docs/user_guide/diseases/` |
| GUD | `stisim/diseases/gud.py` | `docs/user_guide/diseases/` |

## Networks

Source lives in `stisim/networks/`.

| Network | Source |
|---|---|
| Base | `stisim/networks/base.py` |
| Male-female | `stisim/networks/mf.py` |
| MSM | `stisim/networks/msm.py` |
| FSW | `stisim/networks/fsw.py` |
| Layered | `stisim/networks/layered_networks.py` |
| Matchers | `stisim/networks/matchers.py` |

Docs: `docs/user_guide/networks/`.

## Interventions

Source lives in `stisim/interventions/`.

| Area | Source |
|---|---|
| Base | `stisim/interventions/base_interventions.py` |
| HIV | `stisim/interventions/hiv_interventions.py` |
| Syphilis | `stisim/interventions/syphilis_interventions.py` |
| Gonorrhea | `stisim/interventions/gonorrhea_interventions.py` |
| BV | `stisim/interventions/bv_interventions.py` |
| Pregnancy risk | `stisim/interventions/pregnancy_risk.py` |

Docs: `docs/user_guide/interventions/`.

## Connectors

Source: `stisim/connectors/`.
Docs: `docs/user_guide/connectors.md`.

## Analyzers

Source: `stisim/analyzers.py`.
Docs: `docs/user_guide/analyzers.md`.

## Demographics

Source: `stisim/demographics.py`.

## Care-seeking

Source: `stisim/care_seeking.py`.
Docs: `docs/user_guide/care_seeking.md`.

## Calibration

Source: `stisim/calibration.py`.
Docs: `docs/user_guide/calibration.md`.
Note: for actual calibration workflow, hand off to the `calib:` plugin.

## Logistics

Source: `stisim/logistics/`.

## HIVsim (slim HIV subpackage)

Source: `hivsim/sim.py`.
Docs: `docs/user_guide/hivsim.md`.
Examples: `hivsim_examples/simple/`, `hivsim_examples/zimbabwe/`.

## Tutorials

`docs/tutorials/`:
- `tut_intro.qmd` — introduction
- `tut_interventions.qmd` — building interventions
- `tut_cotransmission.qmd` — coinfection dynamics
- `tut_calibration.qmd` — calibration mechanics
- `tut_results.qmd` — reading outputs

## Worked examples

`docs/examples/`:
- `art_interruptions.qmd`
- `dynamic_debut.qmd`
- `partner_notification.qmd`
- `pmtct_scenario.qmd`
- `pregnancy_risk_modifier.qmd`
- `vmmc_costing.qmd`

## Tests

Source: `tests/`. Consult for expected behavior when semantics are unclear from source alone.

## External

- Starsim itself: the framework STIsim is built on. STIsim's `sti.Sim` extends Starsim's `ss.Sim`. Reference Starsim's docs and source for framework-level questions.
