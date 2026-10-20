# HIV and ART as separate declared modules

Design note for porting STIsim's HIV to the Starsim 4 declared style, and for separating ART from HIV
([stisim#335](https://github.com/starsimhub/stisim/issues/335), "Better logic for ART vs HIV").

## The split

- **`sti.HIV`** is the natural history and transmission only. It has no ART states, parameters or branches.
- **`sti.ART`** is a small declared module of its own: treatment states, durations and parameters, plus its effects on
  HIV, written as `ss.modify()` modifiers keyed on ART's states. HIV's code never asks "is this person on ART?".
- **PrEP and VMMC** are Starsim products (`ss.Vaccine`) delivered to a stock target. Their efficacy is on HIV acquisition (`rel_sus`).
- **Testing** stays a coded intervention (`sti.HIVTest`, plus `ANCTest`, `InfantHIVTest`, `PartnerNotification`).
  It sets HIV's knowledge-of-status arrays (`diagnosed`, `ti_diagnosed`) and refers positives to ART.

### HIV: states and transitions

```
S -> acute                   infection: beta by network (sexual per act; maternal and breastfeeding per month)
acute -> latent              latency, after dur_acute
latent -> falling            progression, after dur_latent
falling -> dead              aids_death, after dur_falling (CD4 reaches ~0)
acute|latent|falling -> dead cd4_death: a hazard by CD4 band (cd4_death_rates), x rel_death, x rel_death_f for women
```

CD4 is an `ss.attr`. `cd4_start` and `cd4_latent` are drawn at infection. `cd4` is derived by stage:

- acute: falls from `cd4_start` to `cd4_latent`;
- latent: `cd4_latent`;
- falling: falls linearly to 1 over the stage.

HIV also owns these modifiers:

- `rel_trans`: `rel_trans_acute` in acute; `rel_trans_falling` where CD4 < 200. Each is drawn once per person at
  infection; STIsim 3 drew it again every step.
- `rel_sus`: `rel_sus_age`.

Seeded infections are back-dated with `ss.spread(init_prev, over=(5 months, 10 years))`, as in STIsim 3's `dist_ti_init_infected`.

### ART: states and transitions

```
naive -> on          start: undriven, moved only by ART deliveries
off   -> on          start (restarting)
on    -> off         stop, after dur_on_art (x rel_dur_on_art, x dur_on_art_trend at the time of starting)
off   -> dead        death, after the post-ART survival time (STIsim 3's ti_zero after stopping)
```

Attributes:

| Attribute | Kind | Set |
|---|---|---|
| `cd4_preart` | recorded from `hiv.cd4` | at start |
| `cd4_potential` | recorded | at the first start only |
| `cd4_stop` | recorded | at stop |
| `suppressed` | drawn, `ss.bernoulli(p_effective_art)` | at each start |
| `care_seeking` | trait | once per person |

ART's effects on HIV are modifiers, where `on` means `art.on`:

| Effect | Modifier |
|---|---|
| Transmission | `ss.modify(hiv.rel_trans, multiply=1 - eff*ramp(on.time), where=on)`; `eff` is `effective_art_efficacy` or `nonsupp_art_efficacy`, by `suppressed` |
| Natural history replaced | `ss.modify(hiv.progression, multiply=0, where=~naive)`, the same for `hiv.aids_death` |
| CD4 on ART | `cd4` is set to logistic recovery from `cd4_preart` to `cd4_potential`, by `on.time` |
| CD4 after stopping | `cd4` is set to a linear decline from `cd4_stop` to 1, by `off.progress` |
| Mortality | `ss.modify(hiv.cd4_death, multiply=rel_art_mortality_*, x art_death_age, x art_death_dur, where=on)` |
| MTCT | `ss.modify(hiv.rel_sus, multiply=1 - pmtct_efficacy, where=ss.mother(on) \| ss.mother(on, born=True))` |

Delivery is one or more `ss.treat(ss.Treatment(moves=art.start), ...)` built from `coverage=`:

- **No coverage** (STIsim 3's default): routine delivery to those referred and ready (after `dur_dx2tx`), with
  `coverage=art_initiation`. Each person is offered ART once, so people who stop do not restart. This is STIsim 3's behaviour.
- **A target** (a number or a share, constant or by year; `n_art` and `p_art`): `target=` (a stock), with
  `priority=-cd4/care_seeking` and `above=art.stop`. People over the target are stopped, highest CD4 first.
- **Stratified by age and sex:** one delivery per stratum, each with its stratum in `where=` and `among=`.

`vls_coverage` corrects `suppressed` among those on ART to a target, by adherence rank. This is coded, and agents only.

### Where each STIsim 3 HIV state and parameter goes

| STIsim 3 (`HIV`) | STIsim 4 |
|---|---|
| `acute`, `latent`, `falling`, `ti_*`, `ti_zero`, `ti_dead` | `hiv.acute`, `hiv.latent`, `hiv.falling` (states); `ti_<state>` comes free; deaths are flows |
| `cd4`, `cd4_start`, `cd4_latent` | `hiv` attributes |
| `cd4_nadir` | dropped (unused) |
| `cd4_preart`, `cd4_potential`, `cd4_postart` | `art` attributes (`cd4_postart` becomes `cd4_stop`) |
| `on_art`, `art_naive`, `art_discontinued`, `ti_art`, `ti_stop_art` | `art.on`, `art.naive`, `art.off`; `ti_on`, `ti_off` |
| `on_effective_art`, `on_nonsuppressive_art` | `art.on` & `art.suppressed` (as properties of the same names on `art`) |
| `diagnosed`, `ti_diagnosed` | stay on `hiv` (knowledge of status, not treatment) |
| `care_seeking`, `baseline_care_seeking`, `maternal_care_scale` | `art` (used only for ART priority) |
| `on_prep`, `prep_*` | the PrEP product (`sti.Prep`); `prep.on_prep` |
| `circumcised`, `circ_traditional`, `ti_circumcised`, `eff_circ`, `eff_circ_traditional` | the VMMC products (`sti.VMMC`); `vmmc.circumcised`, `VMMC(eff_circ=...)` |
| `art_cd4_growth`, `art_cd4_pars`, `effective_art_efficacy`, `nonsupp_art_efficacy`, `time_to_art_efficacy`, `p_effective_art` | `sti.ART` |
| `dur_on_art`, `rel_dur_on_art`, `dur_on_art_trend`, `dur_post_art`, `dur_post_art_scale_factor` | `sti.ART` |
| `rel_art_mortality_effective`, `rel_art_mortality_unsupp_m/f`, `art_death_age`, `art_death_dur` | `sti.ART` (the invariant check is in `ART.init_pre`) |
| `cd4_death_bins/rates`, `rel_death`, `rel_death_f`, `include_aids_deaths`, `rel_trans_*`, `rel_sus_age`, betas, `init_prev*` | `sti.HIV` |
| `init_diagnosed` | `sti.HIV` (sets `diagnosed`, and refers to ART at the start) |
| `ART.pmtct_efficacy`, `art_initiation`, `vls_coverage` | `sti.ART`, as before |
| `HIVTest.dur_dx2tx` | `sti.HIVTest`, as before (a referral delay, recorded on ART as `ti_ready`) |

### Results

- **HIV:**
  - `new_infections`, `n_infected`, `prevalence`, `incidence`, with sex and age;
  - `prevalence_15_49`, FSW and client results;
  - routes: `new_infections_sex`, `_mtct`, `_prenatal`, `_postnatal`;
  - `new_deaths` and `cum_deaths` (HIV's own deaths, plus ART's post-ART deaths);
  - `new_diagnoses`, `n_diagnosed`.
- **ART:**
  - `n_on_art`, `p_on_art`, `new_agents_on_art`, `n_on_art_pregnant`, `n_off`, `death`;
  - `sim.results.hiv.n_on_art` becomes `sim.results.art.n_on_art`.

## Alternatives considered

1. **ART states inside HIV** (`on_art`/`post_art` as HIV states). This is the Starsim test model (`tests/test_attrs.py`) and the design review's port.
   It is clean to declare, but it is exactly what #335 wants to undo: the ART parameters and the post-ART CD4 sit in HIV.
2. **ART as a product only** (`ss.Treatment(efficacy=dict(transmission=, mortality=))`). This is clean for the effects, but:
   - dropout is exponential only;
   - there is no post-ART state or death clock;
   - `mortality=` acts on every outflow, including background death;
   - per-person suppression has no slot.
   The design review reached the same verdict.
3. **ART pauses natural history, which resumes after stopping** (Spectrum's assumption). Epidemiologically this is the
   most natural; the post-ART trajectory would then be HIV's own. But it is blocked in Starsim 4 today:
   - a dwell paused with `multiply=0` loses its remaining time (U3);
   - `X.progress` keeps running while a dwell is paused, so CD4 jumps on resume (U4).
   This is the design to move to once those are fixed. ART would then lose its `death` transition and its post-ART CD4 modifier.
4. **Chosen: ART as its own declared module whose modifiers replace HIV's natural history from the first start**, as
   STIsim 3 does ("ART nullifies all future dates in the natural history"). It keeps STIsim 3's behaviour, and puts every ART state and parameter in ART.

## How PrEP, VMMC and testing attach

- **PrEP:** `ss.Vaccine(efficacy=prep_eff*prep_adh, waning=prep_dur)`, delivered by `ss.vaccinate`.
  - The target is a share of the eligible and uninfected, ranked by a per-person willingness.
  - Holders of another PrEP product are excluded.
- **VMMC:**
  - Program VMMC is `ss.Vaccine(efficacy=eff_circ)`, delivered to a stock target among eligible men, ranked by willingness.
  - Traditional circumcision is a routine delivery at `traditional_age` to men before sexual debut, with `coverage=traditional_prob`.
- **Testing:** `HIVTest` is a coded test that sets `hiv.diagnosed` and refers the positives to ART, with the delay.
  ART's eligibility is diagnosed, infected, ready, and naive or off.

## What is lost or changed compared with STIsim 3

- **Courses and durations:**
  - PrEP courses end after an exponential time with mean `prep_dur`, not exactly `prep_dur`. Renewal is a new delivery to the stock target.
  - A dwell given as a fixed number becomes exponential. Every HIV and ART duration keeps STIsim 3's distribution, written explicitly.
- **Transmission:** `rel_trans_acute` and `rel_trans_falling` are drawn once per infection, not every step. The mean is the same.
- **Removed:** `cd4_nadir`, which nothing read.
- **CD4 at death:** CD4-banded deaths come from `ss.by(cd4=...)` on a stratification by CD4, with STIsim 3's bands and rates.
- **Timing:** stopping ART (`dur_on_art`) is a dwell of the `on` state, drawn at each start, as before. Restarts from `off` re-draw it.
- **Agents only:** the declared model runs as agents only (see below).

## Compartmental runs

The sim does not run compartmentally for free. These block it:

- the STIsim networks, which are agents only;
- `ss.attr` (agents only for now; `test_hiv_cd4_ode` is an xfail in Starsim);
- `ss.mother()` for PMTCT;
- the coded tests and the `vls_coverage` correction;
- per-agent modifier amounts (age and duration multipliers on ART mortality).

The declared HIV and ART, without the STIsim networks and the coded tests, would be a candidate once attributes run as compartments.

## Implementation notes

- ART is a declared module and a cascade (`class ART(ss.Disease, ss.cascade)`). Its deliveries are made in `__init__` as the cascade's steps, so the sim has them before it is built. Its modifiers are added in `init_pre()`.
- Under a coverage target, people who became ready on the step (and take up ART, `art_initiation`) start first, then the rest, lowest CD4 over care seeking first. This is STIsim 3's rule, and it matters: ranking everyone by CD4 alone treats sicker people and averts more deaths (in the Zimbabwe demo, 25% fewer HIV deaths).
- `cd4_preart` and `cd4_potential` are set by a hook (`ss.func`) when people start ART, from HIV's CD4 without ART for a first start, else as at the end of the last step. `cd4_stop` is ART's own on-ART CD4 at the move.

## Starsim friction (upstream changes that would remove the shims)

Worst first. Line numbers are in `/home/robyn/starsim-v4` at a0735405.

1. **A dwell paused by `ss.modify(tr, multiply=0)` never resumes.**
   - Where: `solvers/agents.py:1811` clamps the speed to 1e-9, so the remaining time is scaled by 1e9, and `ti_exit_<state>` is an int16 `StepArr` that saturates at 32767; the remaining time is lost.
   - Repro: `ss.modify(hiv.aids_death, multiply=0, where=art.on)`, then stop ART: the person stays in `falling` forever.
   - Fix: keep the remaining dwell (in steps, as float) when the speed is 0, and resume from it.
   - This blocks the "ART pauses natural history" design (alternative 3).
2. **`X.progress` runs on the wall clock, not the dwell.**
   - Where: `solvers/agents.py:897`, `(ti - start)/(end - start)`; when a modifier rescales the exit, `progress` jumps.
   - Fix: progress as the share of the (rescaled) dwell elapsed.
3. **A target delivery with `moves=` and `where=` is not limited to the transition's sources.**
   - Where: `interventions.py:2058` (`eligible_now`): with `where=`, `default_eligible()` (the sources) is not intersected. People already on ART are chosen, counted as `unneeded`, and use up the gap, so the target is never met.
   - Repro: `ss.treat(ss.Treatment(moves=art.start), target=..., where=hiv.infected)` holds at 60% of a 100% target.
   - Fix: always intersect with the sources of `moves=`.
   - Shim: `where=` includes `naive | off_art` (`ART.ready()`).
4. **Products and routine deliveries made in a module's `init_pre()` get no dimension.**
   - Where: `sim.py:325`. Dimensions (`has_<product>`, `offered_<delivery>`) are made in `_init_solver()`, before `init_pre()`, so a delivery a module adds then fails with `KeyError: not stratified by "offered_..."`.
   - Fix: a way for a module to declare the deliveries and modifiers it brings, collected before the solver (`BaseSTI.add_modifiers` has the same TODO).
   - Shim: ART and the prevention programs are `ss.cascade`s, made in `__init__`.
5. **An attribute of another module, read during a move, sees the mover in its new state.**
   - Where: `solvers/agents.py:1022`. `_leaving` covers only the moving module's own states, so another module's `where=` on the mover is evaluated after the move. `ss.attr(hiv.cd4, on=art.stop)` records CD4 without ART's modifier (438 for someone at 124 on ART).
   - Also: an attribute of another module in the class body loses its module when the module is made (`module=None`).
   - Fix: apply `_leaving` across engines during a move, and keep the module of foreign attributes.
   - Shim: `ART.on_start()` hook, and `cd4_stop` recorded from ART's own expression.
6. **`ss.modify(attr, set=...)` is refused, though the agent engine implements it.**
   - Where: `connectors.py:334` (raises for anything but a parameter); `Agents.modify_attr()` applies set/multiply/add.
   - Fix: allow `set=`/`add=` on attributes.
   - Shim: `sti.set_attr` (in `hiv_interventions.py`).
7. **A modifier's amount cannot read the unmodified value of the quantity it modifies.** Reading `hiv.cd4` inside a modifier of `hiv.cd4` recurses.
   - Fix: e.g. `ss.modify(hiv.cd4, set=lambda base: ...)`, or `hiv.cd4.base`.
   - Shim: `ART.cd4_natural()` temporarily removes ART's modifiers.
8. **There is no integrated attribute (one whose rate of change a modifier sets).** "ART changes CD4's rate of change" cannot be said; CD4 must be closed-form by state.
9. **An amount mixing another module's attribute with `ss.by()` is evaluated for different agents.**
   - Repro: `multiply=art.suppressed*0.25 + ss.by(sex=...)*(1 - art.suppressed)` gives `operands could not be broadcast (3,) (3600,)` during back-filling.
   - Shim: a function of the sim.
10. **`ss.spread()` takes only a scalar share.** A per-person probability (initial prevalence by risk group) needs the undocumented `burnin_spread()` hook, which back-dates over `[0, steps)`, not a window.
11. **Pars of a declared module must be passed in `pars=`**, unlike the kwargs that `BaseSTI` filters. `ss.Disease(art_initiation=1.0)` raises `UnknownParameter`, listing pars without the one given.
12. **An undriven transition cannot branch with `p=`**, so suppression at ART start is an attribute (`suppressed`) rather than two states.
13. **A declared result cannot divide a sum of flows by a state** (`(infection + prenatal + postnatal)/S`: `FlowSum / float` TypeError). Shim: `incidence` in `finalize_results()`.
14. **A distribution par times a par in `dur=` raises a confusing `DuplicateDeclaration`.** Shim: `rel_dur_on_art` is a modifier of `stop`.
15. **A declared module in `interventions=` needs `has_product = False`.**
16. **The automatic result of a `BoolState` defined with `define_states()` on a declared disease stays 0.** Shim: `diagnosed` is a `BoolArr`, with HIV's own `n_diagnosed`.
17. **Run sims with declared HIV and ART need a recursion depth of about 500 to pickle**, which fails in `ss.parallel()` workers under pytest-xdist. Shim: `stisim/__init__.py` raises the recursion limit to 3000.
18. **Per-route infection counts need an infection log.** Shim: transmission is declared as three transitions (sexual, prenatal, postnatal); `sw_stats` keeps an `InfectionLog`.
19. **A beta dict with no networks gives a bare `TypeError`** (`dict * float`), not a coded error.
20. **Results written by a declared module in `update_results()` must be in people** (`ss.results.people_scale()`), and people per agent is private (`sim._current_scale`).
21. **`describe()`/`sim.given`:** `sti.Sim` (which calls `ss.Sim.__init__` without its parameters) left `sim.given` empty, so Starsim turned aging off for every `sti.Sim` without demographics modules (`dimensions.py:789`). This is fixed in STIsim; Starsim could respect a `use_aging` that is already set.
