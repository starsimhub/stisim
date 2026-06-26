# What's new

All notable changes to the codebase are documented in this file.

## Version 1.5.8 (2026-06-25)

### Logistics
- Add a `stisim.logistics` package for supply-constrained interventions: `Product` (cost, delivery mode, time-varying efficacy), `Supply`/`Supplies` (quantities and accrued cost), `ProductCategory`/`DeliveryMode` enums, and the abstract `SuppliedIntervention` base that distributes a product across eligible agents subject to available supply. (#492)

### Interventions
- `PartnerNotification`: optional edge-type / partner-sex stratification. `step()` now walks current-channel edges preserving edge-type info and exposes `current_partner_edges` (`{partner_uid: [edge_type_int, ...]}`) to callables on `p_notify_current.p` / `p_attends_current.p`. New `sti.pn_rates(rates)` helper builds such a callable from a dict spec — either `{'stable': 0.20, 'casual': 0.10}` or per-partner-sex `{'stable': {'f': 0.80, 'm': 0.50}, ...}`. Per-edge probabilities for multiply-reached partners are summed and capped at 1. Scalar `p_notify_current` / `p_attends_current` Bernoullis are unchanged (backward compatible). (#505)

### Documentation
- New user-guide pages for analyzers (`docs/user_guide/analyzers.md`) and care-seeking (`docs/user_guide/care_seeking.md`).
- Expand the interventions user guide with syndromic management, antenatal/infant (PMTCT) screening, partner notification, and pregnancy-driven risk reduction.
- List all gallery examples in the sidebar; fix syphilis API drift that broke the docs build.

### Tests
- Add `pn_rates` stratification tests.
- Regenerate `baseline.yaml` / `benchmark.yaml` for starsim 3.4.0 (`n_female`).

## Version 1.5.7 (2026-06-15)

### Diseases
- Syphilis: add treponemal (`trep`) and non-treponemal (`nontrep`) serology states with seroconversion/reversion dynamics (window period, ~80% lifelong trep persistence, non-trep titre decline + post-treatment reversion). New `trep_prevalence`/`nontrep_prevalence` results (overall, by sex, by age×sex) plus `sexually_transmissible`/`symptomatic`/`primary` stage prevalences. Consolidated the duplicate `active` property into `symptomatic` (removed `active_prevalence` and the `*_prevalence_15_64` results). Prevalence denominator age band is now configurable via `age_range`, exposed as a `BaseSTI` constructor argument. (#500)

### Networks
- Add `fsw_mf_conc_mult` on `StructuredSexual` to scale FSW non-sex-work concurrency (values <1 = fewer non-commercial partners). `MFNetwork.set_concurrency` is now negative-binomial-aware: it reparameterizes the target mean correctly instead of silently no-opping when `concurrency_dist` is `ss.nbinom`. (#500)
- Fix the three MSM networks (`AgeMatchedMSM`, `AgeApproxMSM`, `MSMScaleFreeNetwork`): all now honour a participation filter and no longer crash at init. The participation parameter is renamed `msm_share` → `p_msm`. (#502)

### Interventions
- `SyndromicManagement` and `SymptomaticTesting` now take disease **names** (e.g. `diseases=['ng', 'ct']`) rather than module objects, resolved against `sim.diseases` at init. Objects passed at construction were deep-copied into orphans the sim never stepped, giving stale state and dropped result-writes. (#503)

### Analyzers
- Add `sti.PartnershipFormationAnalyzer`. Networks now record the timestep each relationship ends (`None` if still active at sim end), enabling efficient post-run partnership analysis from a shared results store. (#504)

### Results
- Remove the unimplemented `partners_12` / `n_partners_12m` references. (#473, #507)

## Version 1.5.6 (2026-06-02)

### Networks
- **Breaking — default behaviour change.** Split `stisim/networks.py` into a `stisim/networks/` sub-package (`base.py`, `mf.py`, `fsw.py`, `msm.py`, `layered_networks.py`, `matchers.py`). Add a `match_method` string parameter on `MFNetwork`, dispatched through a `MATCHERS` registry of pair-formation algorithms (accepts string or callable). The previous default (`sort_bisect`) collapsed the M-F age gap to ~0 because it didn't honour `age_diff_pars`. Eight matchers ship in the registry. Removes `stisim/pfa_variants.py` and the seven `MFNetwork_*` subclasses. (#472)
- New `closest_age_tapered_seeking` matcher, **now the default**. Orders females by `desired_age = age + normal(target_gap, stddev)` (drawn each step for under-partnered females) and males by age; uses batch-binary-search to find each female's closest available match, with a taper that reduces seeking for older females to bring mean age gaps into data alignment. Realised mean age gap is now within <5% of target (was negatively biased and growing); target stddev preserved. Downstream calibrated models will need recalibration: HIV prevalence shifts ~18% in Zimbabwe and ~59% in Eswatini at default settings. (#477)
- Add `sti.MSMScaleFreeNetwork`: a preferential-attachment MSM sexual network with continuous-time formation/Markovian deletion (Whittles-2019 S2 kernel). Restricted to post-debut males in this first port; subclass hooks (`_get_pool`, `_mix_node_arrays`, `_mix_weights_row`) let you plug in age-/risk-/λ-weighted variants. Documented as not branching-stable under starsim CRN (global weighted-categorical sampling over a pair catalog). Ported from `starsim_x/BespokeNet`, contributed by Stephen Attwood (Oxford). (#488)

### Interventions
- Add `sti.PregnancyRiskReduction`: during pregnancy, optionally clears FSW status (`fsw_redux`), drops high-risk-group membership to `default_risk_group` (`high_risk_redux`), and/or zeros concurrency (`concurrency_redux`). Ported from `syph_dx_zim`. (#482, #484)

### Documentation
- New ART state diagram in the docs illustrating how ART intersects the HIV mortality cascade. (#471)

### Tests
- Regenerate `tests/baseline.yaml` for the new default matcher.
- Add `tests/update_baseline` script (mirroring hpvsim's pattern) to regenerate both `baseline.yaml` and `benchmark.yaml` in one command.

## Version 1.5.5 (2026-05-15)

### Sim / parameter routing
- Add `sti.route_pars` utility (in `stisim/utils.py`) that auto-discovers par categories (sim, sti, nw, dem) from the par-class registry and routes flat kwargs accordingly. Cross-category clashes broadcast to all matches with a printed note; unknown keys raise when `strict=True`. `sti.Sim.separate_pars` now delegates to it and shrinks from ~95 to ~30 lines. (#451, #452)
- Auto-add coinfection connectors when `diseases=[...]` is supplied without an explicit `connectors=`. Lookup tries `sti.<d1>_<d2>` then `sti.<d2>_<d1>`; pairs with no registered class are silently skipped. Replaces the previous `NotImplementedError` stub. (#442, #453)
- Route flat connector pars (e.g. `rel_sus_hiv_syph=3.5`) through to the matching connector instance. (#442, #453)

### Interventions
- Replace `PartnerNotification` with a `PriorPartners`-based implementation. Two channels (current sexual network, optional prior-partner recall network) each with separate notification × attendance probabilities. Tracks notifications and attendance per-channel. The previous transmission-graph-based version is removed. New gallery example `docs/examples/partner_notification.qmd`. (#457)
- Add `dur_dx2tx` parameter to `HIVTest` (default `ss.constant(0)`) for the delay between HIV diagnosis and ART start. `HIVTest` now schedules `hiv.ti_art = ti + delay` on positive results; `ART` becomes passive and starts agents whose `ti_art == ti & ~on_art`. `ANCTest` and initially-diagnosed cases schedule with no delay (backward compatible). (#398, #464)
- Fix `ART` stratified coverage to allocate per-(age, sex) stratum independently. Previously, stratified coverage data was parsed correctly but then collapsed to a single aggregate target before correction, washing out the age/sex differentials in the input. New `compute_stratum_targets` helper in `interventions/utils.py`. (#463)

### HIV
- `new_agents_on_art` result now gates on `on_art` (since `ti_art` may also hold a scheduled future start). (#464)

### Documentation
- New gallery example `docs/examples/vmmc_costing.qmd`: attach unit costs to `sim.results`, discount, and compute an ICER for VMMC vs no-VMMC, framed as a quick STIsim analog of Bansi-Matharu et al., [Lancet Global Health 2023;11(2):e244–e255](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(22)00515-0/fulltext). (#455, #465)

## Version 1.5.4 (2026-05-04)

### Networks
- Split `StructuredSexual` into reusable `MFNetwork` (heterosexual) and `SWNetwork` (sex-work) networks built on a new `BaseNetwork`. `StructuredSexual` is preserved as the combined default. (#307, #445)
- Window sex-work participation by per-agent `age_sw_start` + `dur_sw` instead of lifetime; `fsw` / `client` / `age_sw_stop` are now derived properties. (#307, #445)
- Add `condom_smoothness` parameter forwarded to `sc.smoothinterp` for time-varying condom data. (#445)
- Replace shared `debut` / `debut_pars_f` / `debut_pars_m` with per-sex `debut_f` / `debut_m` distributions; supports scalar shortcuts, full `ss.Dist` instances, and callables. (#397, #420, #423)

### Interventions
- Add `SyndromicManagement` to core stisim, consolidating duplicated implementations from `anc_sti_screening` and `stisim_vddx_zim`. `cervical_diseases` parameter generalises the cervical-symptom check. (#369, #443)
- Add `ANCTest` for single-visit ANC screening (auto-detects HIV + syphilis; per-disease sensitivity; routes to treatment maps). Raises if disease names are mistyped. (#370, #444)
- Add `InfantHIVTest` for newborn/infant testing scheduled by `ANCTest` or by `HIVTest` via the maternal network. (#370, #444)
- Decouple `Prep` eligibility from clinical filters: user-supplied `eligibility=` callables now express targeting only; `~hiv.infected` and `~on_prep` are applied unconditionally. (#330, #441)
- Refactor `Prep` legacy `years=` handling: pop from kwargs before `update_pars` instead of reading back from `self.pars`. (#444)

### HIV
- Update acute HIV defaults to Bellan 2015 central estimates: `rel_trans_acute=N(5.3, 0.5)`, `dur_acute=LogNorm(1.7 mo, 1 mo)`. EHM drops from 15 to ~7.3. (#396, #426)
- Add `rel_sus_age`: list of `(age_lo, age_hi, sex, multiplier)` tuples for age/sex-stratified susceptibility. Default `None` preserves existing behaviour. (#395, #427)
- Remove `BaseSTI.infect()` override; inherits parent `ss.Infection.infect()` for free mixing-pool compatibility. (#412, #421, #422)

### Sim / parameter routing
- `process_networks` / `process_stis` / `process_demographics` now raise on collisions instead of silently winning; `separate_pars` raises on `pars ∩ kwargs` collisions in both `sti.Sim` and `hivsim.Sim`. (#423)
- Improve "ambiguous parameter" error messages to name the offending par. (#423)
- `hivsim.Sim` pulls network-level kwargs from user input into the default `StructuredSexual` instance. (#423)
- Remove `BreastfeedingNet` from `sti.Sim` default networks; only relevant for diseases with postnatal transmission (HIV, syphilis). (#406, #444)

### Bug fixes
- Fix birth data key mismatches in `process_demographics` when `use_pregnancy=False`: `'births'` → `'cbr'`/`'birth'` keys, file path, and metadata columns. (#411, #428)

### Documentation
- Switch docs build from MkDocs to Quarto; convert tutorials from `.ipynb` to `.qmd`. (#439)
- Documentation audit and uplift via IDM standards plugin (folder READMEs, API reference, navigation polish). (#435)
- Add new publication entry: Stuart et al., "Estimating the value of novel tests for active syphilis in Zimbabwe…" forthcoming in *STD*. (#447, #448)

### Tests
- Add `tests/devtests/devtest_sw_networks.py` with 10 integration tests covering MF only, SW only, MF+SW modular, `StructuredSexual`, implicit-SW thresholds, windowed entry/exit, and all `sti.Sim` / `hivsim.Sim` / `hivsim.demo` creation patterns. (#445)
- Add `test_shorter_sw`: shorter SW participation window → fewer FSW-attributable transmissions. (#445)
- Remove fragile `len(sim.diseases/networks/interventions)` count assertions in `test_sim`. (#444)
- Remove fragile `test_doubling_hiv_sexual_beta`. (#445)
- Regenerate `baseline.yaml` to reflect default behaviour changes (windowed SW, BreastfeedingNet removal, Bellan acute pars). (#426, #444, #445)


## Version 1.5.3 (2026-04-17)

### Interventions
- Extract shared coverage-targeting logic into `interventions/utils.py`: `parse_coverage`, `age_sex_mask`, `compute_coverage_target` (#328)
- Replace `_parse_age_bin` with starsim's `ss.parse_age_range` and `ss.apply_age_range` (#326)
- Expose `smoothness` parameter for coverage interpolation (default 0 = linear) (#327)
- Support mixed n/p coverage format: historical absolute numbers transitioning to projected proportions within a single DataFrame or dict (#383)
- Remove deprecated `coverage_data` and `future_coverage` parameters; use `coverage` instead (#383)
- Refactor PrEP to use `parse_coverage` (#327)
- Add `pmtct_efficacy` parameter to ART (default 0.96), replacing hardcoded `rel_sus=0` for prenatal MTCT protection (#404)
- Extend PMTCT protection to breastfed infants via BreastfeedingNet (#407)
- Add ANC HIV testing pattern using `HIVTest` with pregnancy eligibility (#322)

### HIV / MTCT
- Add `new_infections_mtct`, `new_infections_sex`, `new_infections_prenatal`, `new_infections_postnatal` results to BaseSTI (#325)
- Implement breastfeeding HIV transmission via BreastfeedingNet with `beta_breastfeed` parameter (#323)
- Add `p_diagnosed_pregnant` result tracking ANC testing coverage (#322)
- Include BreastfeedingNet in `hivsim.Sim` defaults (conditional on Pregnancy module)

### Bug fixes
- Fix chlamydia `eff_condom` default: updated to 0.4 based on linked citation (#391)
- Fix `hivsim.Sim`: conditionally add BreastfeedingNet only when Pregnancy is present (#413)
- Fix disease count assertion in `test_sim_creation` (#413)

### Documentation
- Add examples gallery with first entry: Modeling ART interruptions (#394)
- Add PMTCT section to HIV user guide (ANC testing, prenatal/postnatal protection)
- Document HIV transmission route results (MTCT, prenatal, postnatal)
- Add docstrings across interventions, diseases, networks, and analyzers (#409)

### Tests
- Add `test_pmtct`: 8-sim combinatorial test for ANC testing, PMTCT efficacy, and breastfeeding duration
- Add `test_mtct`: prenatal + postnatal MTCT consistency checks

## Version 1.5.2 (2026-04-06)

### Bug fixes
- Fix `default_build_fn` to apply `rand_seed` when `reseed=True`, so each calibration trial uses its own seed
- Fix timepar handling in data loaders (remove unnecessary `TimePar` catch)
- Fix stray plot appearing when running pytest (`test_zimbabwe` called `hivsim.demo` with `plot=True`)

### Improvements
- Add `HIV.plot()` with curated 6-panel view (new infections, deaths, prevalence, prevalence 15-49, diagnoses, proportion on ART)
- Pin `starsim>=3.3.1` and `sciris>=3.2.9`

### Tests
- Add baseline regression test and benchmark test using Zimbabwe HIV example
- Refactor test suite with clear file responsibilities: `test_sim.py` (constructor/routing), `test_hiv.py` (HIV scientific validation), `test_stis.py` (non-HIV STI validation), `test_networks.py` (network dynamics)
- Add `hivsim.Sim` constructor tests (defaults, parameter routing, custom modules)
- Add HIV sensitivity tests: `test_par_ranges` checks `beta_m2f`, `init_prev`, `dur_falling`, and `art_efficacy` against both `cum_infections` and `cum_deaths`
- Add `test_prevalence_by_sex` verifying female > male HIV prevalence under defaults
- Add MSM network test
- Rename `test_diseases.py` to `test_stis.py`, `test_hiv_natural_history_verification.py` to `test_hiv.py`
- Merge `test_examples.py` into `test_sim.py`

### Documentation
- Fix calibration tutorial to use disease instances instead of strings
- Add note about `scikit-learn` requirement for `plot_param_importances`

## Version 1.5.1 (2026-03-27)

### Bug fixes
- Fix VMMC eligibility: pass `eligibility` to parent class and use `check_eligibility()` in `step()`, so user-supplied eligibility functions are no longer silently ignored (#353)
- Fix starsim compatibility: convert `ss.years`/`ss.date` to `int` in `get_age_distribution()` for newer starsim `__len__` behavior
- Fix rand_seed placement in beta transmission test (#354)

### Parameter handling
- Remove parameters from `SimPars` that duplicate starsim's `Sim` (`label`, `verbose`, `total_pop`, `pop_scale`, `birth_rate`, `death_rate`, `use_aging`) (#351)
- Add `get_valid_pars()` to detect ambiguous parameters across namespaces (#351)
- Improve error messages for unrecognized parameters with suggestions via `sc.suggest()` (#351)
- Handle `location`/`demographics` conflict in `remap_pars()` (#351)
- `mergepars()` now returns `sc.objdict` for dot-notation access (#351)

### Other changes
- Delete legacy `setup.py` (use `pyproject.toml`)
- Relax `sc.require` to only check starsim, with `die=False`
- Comment out unused `dt` parameter in Migration

### Tests
- Add VMMC tests: reduces male infections, male-only targeting, eligibility targeting (#251, #355)
- Add ART dropout tests: CD4 decline and rel_trans increase after stopping ART (#244)
- Add exponential early-phase epidemic growth test (#260)
- Add increased testing speeds diagnosis test (#253)
- Add no-HIV-without-seeding test (#271)

### Documentation
- Fix `sti.Pregnancy` → `ss.Pregnancy` in Getting Started tutorial (#343)
- Fix tutorials and user guide nav nesting
- Update syphilis docs: exposed stage, treatment/reinfection, key properties
- Fix HIV `dur_on_art` default in docs

## Version 1.5.0 (2026-03-13)

### Interventions
- Refactor ART and VMMC coverage input: accept scalar, dict, or DataFrame formats with flexible column names (#126)
- Add `art_coverage` analyzer for tracking ART coverage by age and sex
- Simplify `art_initiation`: accept plain number, drop `init_prob` backward compat
- Guard pregnancy/maternalnet access in HIV and syphilis interventions so they work without a pregnancy module (#319)
- Document HIV intervention pipeline and ordering (#314)

### HIV
- Add `aids` property to HIV module (cd4 < 200)
- Fix `dur_on_art` being silently scaled by `dur_on_art_trend`: default is now `None` (#336)
- Add time-varying ART duration and wider care-seeking distribution
- Add custom module support to Sim constructor (#318)

### Syphilis dynamics
- Fix `Syphilis.infect()`: use `rel_trans=1` for maternal transmission, stage-specific for sexual
- Fix `NewbornTreatment` to detect MTC-infected babies (susceptible=False, congenital not yet fired)
- Fix congenital syphilis over-counting: clear `ti_*` after events fire, mark babies non-susceptible after MTC
- Add `step_die` to Syphilis to clear states on death
- Add `n_infections` counter and `new_reinfections` result

### Calibration API
- New dot-notation parameter routing: `'hiv.beta_m2f'` automatically finds and sets the right module parameter via `sim.get_module()` (requires Starsim 3.2.0+)
- Support nested parameter format: `dict(hiv=dict(beta_m2f=dict(low=..., high=...)))` alongside flat `{'hiv.beta_m2f': dict(...)}`
- Add `flatten_calib_pars()` to normalize between nested and flat formats
- Add `set_sim_pars(sim, pars)` for setting calibrated parameters on any sim (pre- or post-init)
- Add `Calibration.get_pars(n)` to extract top-N parameter sets as flat dicts
- Add `make_calib_sims()` for creating and running sims from calibrated parameters, with filtering (`check_fn`) and seed replication (`seeds_per_par`)
- Add `Calibration.save()` method for shrink/save workflow (replaces manual shrink + saveobj boilerplate)
- No custom `build_fn` needed -- `default_build_fn` handles all module types automatically

### Documentation
- Add calibration tutorial with ABC philosophy, custom analyzer fitting, production workflow
- Update co-transmission tutorial: HIV-syphilis example with epidemiological explanation of connector effects
- Add custom results section to results tutorial
- Cross-link calibration and results tutorials

### Tests
- Add HIV natural history verification test suite: CD4 decline, transmission doubling, MTCT, AIDS property (#178, #226, #227, #228, #237, #238, #295)
- Add HIV intervention tests: ART coverage formats, duration, effects, parameter sensitivity
- Add `testlib.py` with shared `build_testing_sim()` helper

## Version 1.4.9 (2026-02-24)

- Add `default_build_fn` for calibration: routes parameters by prefix (`hiv_*`, `syph_*`, `nw_*`) to diseases and networks automatically, removing need for custom `build_fn`
- `sti.Calibration` now uses `default_build_fn` when no `build_fn` is provided
- Fix `make_df()` time column to use years from `timevec` instead of integer index
- Fix `parse_study()` to guard `sim_results` reordering when `save_results=False`
- Fix calibration data column names to use dot notation (`hiv.prevalence`)
- Standardize import conventions: `import stisim as sti`, `import hivsim` (no alias)

## Version 1.4.8 (2026-02-23)

- Add `hivsim_examples` package with `simple` and `zimbabwe` pre-configured examples
- Add `hs.demo()` function for quickly running example HIV simulations (e.g. `hs.demo('zimbabwe')`)
- Add `Sim.plot()` override that auto-selects curated HIV result keys when HIV is present
- Add `data` parameter to `sti.Sim` for passing comparison data (e.g. UNAIDS estimates) to starsim's plot overlay
- Fix `process_demographics()` to use `total_pop` from sim_pars when set, and correctly scale age data (values in thousands)
- Fix `separate_pars()` so kwargs override `sim_pars` defaults (e.g. `stop=1995` overrides `sim_pars=dict(stop=2025)`)
- Fix `get_age_distribution()` to handle CSV files without a `year` column

## Version 1.4.7 (2026-02-20)

- Fix bugs in coinfection analyzer: age limit typo, variable name errors for male results
- Remove STIsim Pregnancy module, now superseded by the Starsim Pregnancy module
- Add exposed/incubation period for syphilis (`dur_exposed`), and allow maternal transmission during this stage
- Rename syphilis birth outcome key `active` to `mat_active` to reflect inclusion of exposed stage
- Add `.vscode/` and `*.code-workspace` to `.gitignore`

## Version 1.4.6 (2026-02-20)

- Add sensible `beta_m2f` defaults for all diseases (NG: 0.06, CT: 0.06, TV: 0.1, HIV: 0.05, syphilis: 0.1)
- Add `rel_trans_hiv_ng` parameter to the HIV-NG connector
- Fix connector handling in `Sim.init()` when a single connector is passed
- Set `auto_plot=False` on subpopulation, care-seeking, and detail results so `sim.plot()` shows only high-level results
- Add tutorials: intro (gonorrhea), co-transmission, results/plotting, and interventions
- Add user guide pages for interventions (testing, treatment, ART, VMMC, PrEP)
- Update disease docs with `beta_m2f` defaults and per-act description
- Pre-build font cache in docs CI workflow

## Version 1.4.5 (2026-02-20)

- Fix BV trimester KeyError during initialization
- Add documentation: intro tutorial, disease reference pages, and network guide
- Configure mkdocs to execute notebooks and add user guide structure
- Update README with Python/R install instructions and example repos
- Move devtests to `tests/devtests/` and rename `test_hiv` to `devtest_hiv`

## Version 1.4.3 (2025-12-08)

- Patch to add super calls to init_results for analyzers
- *GitHub info*: PR [160](https://github.com/starsimhub/stisim/pull/160)

## Version 1.4.2 (2025-11-28)

- Patch to ensure that products have different names
- Add `ti_exposed` attribute to HIV module
- *GitHub info*: PR [150](https://github.com/starsimhub/stisim/pull/150)

## Version 1.4 (2025-08-12)

- Add location arg and Sim class
- Update to work with Starsim v3.
- *GitHub info*: PR [148](https://github.com/starsimhub/stisim/pull/148)

## Version 1.3 (2025-06-27)

- Fixes to the pair-matching algorithm within the sexual network to better align partner ages
- Improvements to networks, including analyzers for debut age and partner age differences
- *GitHub info*: PR [143](https://github.com/starsimhub/stisim/pull/143)

## Version 1.2 (2025-06-10)

- Improvements to networks, including analyzers for relationship duration and network degree
- Adds a `PriorPartners` network for recalling past relationships - for use in partner notification
- *GitHub info*: PR [135](https://github.com/starsimhub/stisim/pull/135)

## Version 1.1.2 (2025-06-03)

- Bugfix to calibration class for multisims
- *GitHub info*: PR [138](https://github.com/starsimhub/stisim/pull/138)

## Version 1.1.1 (2025-05-23)

- Bugfixes to calibration class and BV connector
- Replaces the `match_pairs` method of the `StructuredSexual` network with the faster option that was previously in the `FastStructuredSexual` network (now removed).
- *GitHub info*: PR [124](https://github.com/starsimhub/stisim/pull/124)

## Version 1.1.0 (2025-05-13)

- Improvements to the Calibration class: this class now inherits directly from the Starsim calibration class, so users will have access to easier parameter constraints, plotting, flexible fit functions, etc
- Generalization of the coinfection class to handle any two diseases
- Addition of a more detailed BV model
- *GitHub info*: PR [119](https://github.com/starsimhub/stisim/pull/119)

## Version 1.0.5 (2025-05-08)

- Adds results for syphilis transmission by disease stage
- *GitHub info*: PR [112](https://github.com/starsimhub/stisim/pull/112)

## Version 1.0.4 (2025-05-07)

- Adds results for overtreatment among pregnant women
- *GitHub info*: PR [111](https://github.com/starsimhub/stisim/pull/111)

## Version 1.0.3 (2025-04-30)

- Bugfixes for congenital syphilis
- *GitHub info*: PR [85](https://github.com/starsimhub/stisim/pull/85)

## Version 1.0.2 (2025-04-14)

- Bugfixes for syphilis and GUD
- *GitHub info*: PR [83](https://github.com/starsimhub/stisim/pull/83)

## Version 1.0.1 (2025-03-31)

- Track HIV prevalence for 15-49 year olds
- *GitHub info*: PR [79](https://github.com/starsimhub/stisim/pull/79)

## Version 1.0.0 (2024-12-11)

- Updates to work with Starsim v2.2.1
- *GitHub info*: PR [63](https://github.com/starsimhub/stisim/pull/63)

## Version 0.2.0 (2024-11-01)

- Updates to work with Starsim v2.0
- *GitHub info*: PR [62](https://github.com/starsimhub/stisim/pull/62)

## Version 0.1.0 (2024-10-02)

- Collection of updates related to the NG/CT/TV work
- *GitHub info*: PR [59](https://github.com/starsimhub/stisim/pull/59)

## Version 0.0.2 (2024-06-07)

- Initial version of STIsim with structured sexual networks, models of HIV and syphilis, worksflows for model calibration, and interventions for testing and treatment.
- *GitHub info*: PR [22](https://github.com/starsimhub/stisim/pull/22)

## Version 0.0.1 (2024-05-15)

- Pre-release version
