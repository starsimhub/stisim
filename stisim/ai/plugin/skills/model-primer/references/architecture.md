# STIsim architecture — mental model

The purpose of this document is to give the reader a mental model of how STIsim's components compose. Deep semantics live in the source; this file exists so a reader can hold the shape in working memory without paging in the whole codebase.

## Orientation

STIsim is a Python package for agent-based modeling of sexually transmitted infections and HIV, built on top of the Starsim framework. The user-facing entry point is `sti.Sim` (`stisim/sim.py`), which extends Starsim's `ss.Sim` with STI-specific defaults, module selection, and parameter routing. Anything in Starsim is available; the STIsim additions are conveniences and structure on top.

Two ways to think about STIsim:
- **As a library of ready-to-compose modules.** Diseases, networks, interventions, connectors, analyzers, demographics — pick a set, hand them to `sti.Sim`, run.
- **As an opinionated assembly.** `sti.Sim()` with no arguments produces a reasonable STI simulation with defaults for every module slot. Users override what they care about; everything else uses the default.

## The Sim as assembly point

`sti.Sim` composes modules via keyword arguments: `diseases`, `networks`, `demographics`, `interventions`, `analyzers`, `connectors`. Each may be:
- A string name (`'hiv'`, `'syphilis'`, `'zimbabwe'`) — resolved to the matching STIsim or Starsim class and initialized with defaults.
- An instance of a module already configured by the user.
- A list mixing strings and instances.

Per-module parameters travel via matching `<slot>_pars` kwargs, e.g. `disease_pars={'hiv': {'init_prev': 0.1}}`. This decouples the choice of which modules to include from their configuration.

**Two-phase lifecycle.** Construction (`__init__`) collects and clones the modules. Initialization (`init`, before the run) resolves cross- references between modules, wires up the simulation state, and finalizes parameter derivations. If a module needs to reference another module by name (a connector referring to two diseases, an intervention referring to its target disease), the wiring belongs in `init_pre(sim)` — not in `__init__` — because the object you passed into `sti.Sim` is not the object that runs. See the deep-copy gotcha below.

## Diseases

Disease modules live in `stisim/diseases/`. Two base classes matter:

- **`BaseSTI` / `BaseSTIPars`** — shared foundation for every STI: transmission parameters (`beta_m2f`, `rel_beta_f2m`, `beta_m2c`, `beta_breastfeed`, `beta_m2m`), condom efficacy, timestep, initial prevalence rescaling.
- **`SEIS` / `STIPars`** — extends `BaseSTI` with a specific natural- history shape (susceptible → exposed → infected → susceptible) used by chlamydia, gonorrhea, and trichomoniasis.

HIV, syphilis, BV, and GUD have richer natural histories and extend the base directly rather than reusing the SEIS template. Each disease's own file (`hiv.py`, `syphilis.py`, `bv.py`, `gud.py`) is the source of truth for its mechanics.

**A recurring pattern.** Each disease exposes hooks for testing, diagnosis, symptom expression, and treatment. Interventions attach to these hooks rather than reaching into disease internals — that keeps the disease modules focused on natural history and lets interventions compose across diseases.

## Networks

Sexual networks live in `stisim/networks/`. They describe who is partnered with whom, for how long, and how new partnerships form. Building blocks:

- **`base.py`** — the network base classes and shared machinery.
- **`mf.py`** — male-female networks: stable, casual, and marital layers.
- **`msm.py`** — male-male sexual networks.
- **`fsw.py`** — female sex worker and client dynamics.
- **`layered_networks.py`** — the composition that wraps the above into a single population-level structure that diseases transmit on.
- **`matchers.py`** — the algorithms that choose partners subject to the network's constraints.

The layered structure is what `sti.Sim` uses by default; specialized studies (MSM-only, FSW-focused) can drop or reweight layers.

## Interventions

Intervention modules live in `stisim/interventions/`. They fall into two rough shapes:

- **Screening / testing / diagnosis.** Test-then-treat pipelines, point-of-care versus lab, syndromic management, ANC screening.
- **Prevention / treatment / behavior.** ART, PrEP (including LA-PrEP where represented), condom promotion, partner notification, VMMC, BV treatment.

Base machinery is in `base_interventions.py`; per-disease interventions are in per-disease files (`hiv_interventions.py`, `syphilis_interventions.py`, etc.). Interventions consume disease hooks (e.g. "who is symptomatic and unmet?") rather than reading disease state directly, which keeps them portable across diseases.

## Connectors

Connectors live in `stisim/connectors/` and couple two or more disease modules that would otherwise run independently. Typical use: modifying the transmissibility or severity of disease A in the presence of disease B — HIV amplifying syphilis susceptibility, pregnancy modifying STI-specific risk. Connectors are wired at `init_pre` time and read-only against the diseases they connect: they modify parameters or event rates, not the diseases' internal states.

## Analyzers

Analyzers live in `stisim/analyzers.py`. They observe the simulation without mutating it: per-timestep snapshots by age × sex × disease, running counts, custom summary statistics for downstream analysis. Analyzer output lives on the sim's results structure alongside the disease and network results.

## Demographics

Demographics are in `stisim/demographics.py`. They own the background population dynamics: births, deaths, migration, aging. For location- specific simulations the demographics module loads country-level rates and cohort structure from `stisim/data/`.

**Setup gotcha.** Passing `demographics='<location>'` to `sti.Sim` keeps the value as a string until `init`; parameter overrides targeting migration, pregnancy, or death rates need to be routed via the `dem_pars=` kwarg (not via post-hoc parameter mutation that assumes a module instance is already present).

## Care-seeking

STIsim provides a `CareSeeking` module (`stisim/care_seeking.py`) that represents care-seeking propensity as a per-agent value drawn at birth (lognormal, sex-differential by default) and modifiable by conditions like pregnancy. The intent is that testing and treatment interventions consult this module rather than each disease and intervention re-implementing care-seeking logic.

In current practice, most STIsim projects **parameterize care-seeking on the disease modules and testing interventions directly** rather than composing `CareSeeking` into the Sim. The relevant knobs are `p_symp_care` on the SEIS diseases, `rel_test` and similar scaling factors on testing interventions, and per-project multipliers that sit above those (e.g. a `care_seek_mult` in the analysis-specific code). The `CareSeeking` module is available for analyses that want cross-disease consistency without per-module parameterization.

## Timestep and time discipline

The default timestep is monthly (`dt = 'month'`). Two consequences that trip people up:

1. **Sub-monthly windows.** A month advances gestational age by roughly 4.3 weeks. Any process whose eligibility window is narrower than a month (a specific-week ANC screen, a short latent phase) will fire near-zero times unless it is re-cast in month-compatible units.
2. **Event scheduling.** Scheduling an event at `ti + fractional duration` requires the target timestep to be floored (or otherwise coerced to an integer) before equality-comparing to `self.ti`; floating scheduled times silently miss firing when the sim's integer `ti` walks past them.

If a process really needs a finer temporal resolution, either use a smaller `dt` or model the process as a hazard rate rather than a scheduled event.

## Common idioms and gotchas

*Non-exhaustive; extend as they arise.*

- **`sti.Sim` deep-copies modules on construction.** A module handed to `sti.Sim` is cloned into the sim; cross-references to other modules resolved in `__init__` refer to objects that don't end up running. Resolve cross-refs in `init_pre(sim)` instead.
- **Prefer the Starsim `Arr` API.** Use `.values`, `.isnan`, `.notnan`, `.notnanvals` on `FloatArr` / `BoolArr` state variables — these already filter by the sim's live agent uids (`sim.people.auids`). Indexing `arr[sim.people.auids]` manually or inspecting `arr.raw` is a source of subtle bugs.
- **Annual aggregation.** For annual resampling of a `ss.Result`, use `.annualize()` or `.to_df(resample='year')` — the API already distinguishes flow versus stock semantics. Don't groupby-year and hand-roll mean / sum aggregations.
- **Stochasticity discipline.** A single-seed observation is not evidence of a mechanism. Any claim about model behavior needs three or more seeds, and any proposed mechanism must reconcile with known model idioms — expanded treatment can raise incidence when it selects for AMR, for example, so that pattern is not a bug on its own.
- **Result timing subtleties.** Some per-timestep counters are cleared at the top of `step()` before `update_results()` reads them; if a result reads unexpectedly zero, check the order of operations in the disease module rather than assuming the event didn't happen. Prefer the module's own outcome states (e.g. `syph.cs_outcome`) as the primary source when in doubt.

For canonical implementations, docs, and tests, see `canonical-sources.md`.
