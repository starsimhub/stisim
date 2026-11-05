# HIV

HIV in STIsim is modeled with CD4-based disease progression through acute, latent, and late-stage phases, with ART treatment effects.

**Class:** `sti.HIV` | **Alias:** `'hiv'` | **Base class:** `ss.Disease`

**Declared and coded.** The natural history and care are declared states and transitions: `susceptible`, `acute`, `latent`, `falling`, `on_effective_art`, `on_nonsuppressive_art`, and `art_discontinued`, with the transitions `infection`, `mtct_prenatal`, `mtct_postnatal`, `latency`, `progression`, `aids_death` (when the CD4 count reaches zero), `post_art_death`, `hiv_death` and `art_death` (at a rate by CD4 count), and `art_dropout` (stopping ART early, at `rate_art_dropout`, none by default). ART's reduction of transmission at full efficacy (`effective_art_efficacy`, `nonsupp_art_efficacy`) is part of the declared infection rates, so a modifier (`ss.modify(hiv.pars.effective_art_efficacy, ...)`) changes it. The continuous CD4 count, the rise to ART's full efficacy over its first months, susceptibility (PrEP, circumcision, age), care seeking and diagnosis are code in the same class; starting and stopping ART (`start_art()`, `stop_art()`) move people between the declared states. `infected`, `on_art` and `aids` are conditions over the states. Other modules read and change HIV through these: e.g. MIGHTI's competing risks of death read the hazards of `hiv_death` and `art_death` (`hiv.hiv_death.hazard(uids)`), and its adherence model is a modifier of `rate_art_dropout`.

## States and transitions

```
                        ┌─────────────┐
                        │ Susceptible │
                        └──────┬──────┘
                               │ infection
                               ▼
                        ┌─────────────┐
                        │    Acute    │  CD4 declines from ~800 to ~500
                        │  (~3 mo)    │  Transmissibility: 6x baseline
                        └──────┬──────┘
                               │
                               ▼
                        ┌─────────────┐
                        │   Latent    │  CD4 stable at ~500
                        │  (~10 yr)   │  Transmissibility: 1x baseline
                        └──────┬──────┘
                               │
                    ┌──────────┴──────────┐
                    │                     │
                    ▼                     ▼
             ┌─────────────┐       ┌─────────────┐
             │   Falling   │       │   On ART    │  CD4 reconstitutes
             │   (~3 yr)   │       │   (~3 yr)   │  Transmissibility: 0.04x
             │ 8x baseline │       └──────┬──────┘
             └──────┬──────┘              │ ART dropout
                    │                     ▼
                    │              ┌─────────────┐
                    │              │  Post-ART   │  CD4 declines linearly
                    │              └──────┬──────┘
                    │                     │
                    ▼                     ▼
             ┌─────────────────────────────────┐
             │          AIDS death             │
             │       (CD4 reaches 0)           │
             └─────────────────────────────────┘
```

## States and transitions with Treatment

```
                       ┌─────────────────┐
                       │   Susceptible   │
                       └────────┬────────┘
                                │ infection (beta_m2f × rel_beta_f2m)
                                │ condom use, VMMC reduce transmission here
                                ▼
                       ┌─────────────────┐
                       │      Acute      │  CD4 declines from ~800 → ~500 (linear decline)
                       │    (~3 mo)      │  rel_trans: 6× baseline
                       └────────┬────────┘
                                │                         ╔══════════════════════════╗
                                │    ┌────────────────────║   ART initiation can     ║
                                │    │  ART (testing      ║   occur from Acute,      ║
                                │    │   + diagnosis      ║   Latent, or Falling.    ║
                                │    │   required)        ║   Nullifies all pending  ║
                                ▼    │                    ║   stage transitions.     ║
                       ┌─────────────────┐                ╚══════════════════════════╝     
                       │     Latent      │
                       │    (~10 yr)     │
                       │  CD4 stable     │
                       │  at ~500        │
                       │  rel_trans: 1×  │
                       └────────┬────────┘
                                │
              ┌─────────────────┤
              │  no treatment   │  diagnosis + ART initiation
              ▼                 ▼
     ┌─────────────────┐   ┌─────────────────┐
     │    Falling      │   │    On ART       │  CD4 reconstitutes (logistic)
     │    (~3 yr)      │──▶│    (~3 yr ×     │  toward cd4_potential ceiling
     │  CD4: ~500 → 0  │   │  rel_dur_on_art)│  rel_trans: 0.04× baseline
     │  rel_trans: 8×  │   │                 │  (96% reduction, ramps over 6mo)
     └────────┬────────┘   └────────┬────────┘
              │                     │ ART dropout
              │                     │ (stochastic, dur_on_art ~ LogNormal)
              │                     ▼
              │            ┌─────────────────┐
              │            │    Post-ART     │  CD4 declines linearly → 0
              │            │                 │  Rate depends on CD4 at dropout
              │            │                 │  and cd4_potential
              │            └────────┬────────┘
              │                     │
              ▼                     ▼
     ┌─────────────────────────────────────────┐
     │                AIDS death               │
     │   Deterministic: CD4 reaches 0 (ti_zero)│
     │   Stochastic: p_death(CD4) per timestep │
     │     CD4 > 350: 0.3%/yr                  │
     │     200–350:   0.5%/yr                  │
     │     50–200:    1.0%/yr                  │
     │     0–50:      5.0%/yr                  │
     │     ~0:       30.0%/yr                  │
     │   (untreated agents; on-ART agents      │
     │    also die here via a separate,        │
     │    smaller age/sex/adherence/CD4-       │
     │    stratified rate, see below)          │
     └─────────────────────────────────────────┘
```

## On-ART mortality

Agents on ART may still die from HIV-related causes. Mortality while on ART depends on age, sex, CD4 count, and ART adherence status (effective/virally-suppressive vs. non-suppressive ART), in the following directions:

- **Lower CD4 → higher mortality.** Even while on ART, lower CD4 count leads to higher mortality. On-ART mortality is anchored to the same off-ART CD4-based hazard used below (`cd4_death_bins`/`cd4_death_rates`), so it rises sharply as CD4 falls, exactly like the off-ART hazard does. By default, mortality is divided into six CD4 count bins: 0-50, 50-200, 200-350, 350-500, 500-1000, 1000+.
- **Older age → higher mortality.** Older individuals have higher HIV-related mortality than younger individuals. `art_death_age` scales mortality up with age, from 1.0x (under 25) to 1.32x (45+). By default, ART mortality is divided into four age bins: 0-25, 25-35, 35-45, 45+.
- **Non-suppressive ART → higher mortality than effective ART.** Failing to achieve viral suppression carries a mortality penalty on top of the CD4-based risk.
- **Men on non-suppressive ART → higher mortality than women on non-suppressive ART** Failing treatment carries a proportionally worse mortality penalty for men than for women. The non-suppressive/effective mortality ratio is roughly 2x higher for men than for women by default (2.8x vs. 1.4x).

`on_art_death_rate()` computes the rate of `art_death`, **bounded from above by the off-ART CD4-based rate** for the same age/sex/cd4 count demographic subgroup:

```
rate = off_art_rate(cd4) * rel_art_mortality[effective? / sex] * age_mult(age) * rel_death_f?
```
On-ART mortality is bounded by `off_art_rate` such that someone who is on ART will never have mortality hazard higher than if they were off ART.

`rel_art_mortality_effective` (default 0.25) is the fraction of off-ART mortality rate while on effective ART - a 75% mortality reduction relative to being off ART for both male and female. For individuals on non-suppressive ART, there is an additional parameter – `rel_art_mortality_unsupp_m` (default 0.7) for males or `rel_art_mortality_unsupp_f` (default 0.35) for females - which parameterizes how nonsuppressed males are at higher mortality risk than nonsuppressed females.

## On-ART transmission

ART adherence status also affects transmission. Individuals who achieve viral load suppression on effective ART have lower transmission rates than individuals on non-suppressive ART. (That is to say, "undetectable = untransmissible.") 

In `HIV.update_transmission()`, `rel_trans` for an on-ART agent ramps linearly (over `time_to_art_efficacy`, default 6 months, starting from `ti_art`) down to `1 - efficacy`, where `efficacy` is `effective_art_efficacy` (default 0.99, i.e. a 99% reduction — consistent with "U=U": PARTNER/PARTNER2/Opposites Attract found zero linked transmissions from virally suppressed partners across thousands of couple-years) for virally-suppressed agents, or `nonsupp_art_efficacy` (default 0.35, i.e. only a 35% reduction) for non-suppressive agents. In other words: a non-suppressive agent remains substantially infectious despite nominally being "on ART," at roughly `0.65/0.01 = 65x` the transmissibility of an effective-ART agent once both are fully ramped. This mirrors the mortality story above — both mortality and transmission risk are much closer to untreated levels for agents who fail to achieve viral suppression, even though they're recorded as being on ART.

**Where does `nonsupp_art_efficacy = 0.35` come from?** It's anchored to a viral-load-based estimate rather than measured directly (there's no equivalent trial data for "non-suppressive ART" the way PARTNER/PARTNER2 give us "suppressive ART"). Using the Quinn (Rakai 2000) hazard ratio of ~2.45x transmission risk per +1 log10 viral load, and assuming a non-suppressed viral load distribution of log10(VL) ~ N(3.75, 0.8²) relative to untreated chronic infection at ~4.5 log10:

```
2.45^(3.75 - 4.5) * exp(ln(2.45)² * 0.8² / 2) ≈ 0.51 * 1.29 ≈ 0.66
```

giving a relative infectiousness of ~0.66, i.e. efficacy ≈ 1 - 0.66 = 0.34 — close to the shipped 0.35. The `exp(...)` term is a Jensen correction (the hazard is exponential in log10(VL), so integrating over the assumed distribution differs materially from evaluating at its mean). This estimate is sensitive to the assumed VL distribution — a mean of 4.0 gives ~0.16, a mean of 3.5 gives ~0.47 — so `nonsupp_art_efficacy` should be treated as anchored-but-revisitable rather than a hard empirical estimate, and is a reasonable target for recalibration or sensitivity analysis.

### Other Notes on ART

- ART can be initiated from **any** infected stage (Acute, Latent, or Falling)
- On initiation: all pending stage transitions (`ti_latent`, `ti_falling`, `ti_zero`) are cleared for future events
- ART duration is drawn per-agent: `dur_on_art ~ LogNormal(3yr, 1.5yr) × rel_dur_on_art`
  - `rel_dur_on_art` is a calibrated scalar (range 1–20) that stretches mean duration
- After dropout, `ti_zero` is redrawn based on CD4 at dropout and `cd4_potential`
- **No age or sex dependence** in natural history or ART dropout timing. Mortality and transmission *while on ART* are age/sex/CD4/adherence-dependent (see "On-ART mortality"/"On-ART transmission" above); other stages are not.

### What is not modeled

- Viral load / viral suppression (ART = binary effective/non-suppressive, not a continuous measure)
- Age-varying progression rates (natural history is age-independent; only on-ART mortality varies by age)
- Sex-varying natural history (only on-ART mortality varies by sex)
- Re-infection or superinfection
- Duration-since-ART-initiation effects on mortality beyond what's mediated through CD4 (see "On-ART mortality" above)
- EMOD's "terminal failing window" (a temporary loss of transmission suppression in the months before an ART-mortality-table-scheduled death)

## Parameters

Pars of `sti.HIV`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd), which combines these tables with those of the other modules.

### Natural history

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `cd4_start` | normal(800, 50) | CD4 count at infection | Not recorded |
| `cd4_latent` | normal(500, 50) | CD4 count during the latent phase | Not recorded |
| `dur_acute` | lognorm(1.7 mo, 1 mo) | Duration of acute infection | [Bellan 2015](https://doi.org/10.1371/journal.pmed.1001801) |
| `dur_latent` | lognorm(10 yr, 3 yr) | Duration of latent infection (untreated) | Not recorded |
| `dur_falling` | lognorm(3 yr, 1 yr) | Duration of late-stage CD4 decline | Not recorded |
| `p_aids_death` | 1.0 | Probability of dying of AIDS when the CD4 count reaches zero (0 turns AIDS deaths off) | Not applicable |
| `cd4_death_bins` | [1000, 500, 350, 200, 50, 0] | Lower edges of the off-ART CD4 mortality bins (descending) | Not recorded |
| `cd4_death_rates` | [0.003, 0.003, 0.005, 0.01, 0.05, 0.30] | Off-ART annual mortality rate in each CD4 bin | Not recorded |
| `rel_death` | 1.0 | Multiplier on all HIV death probabilities, off and on ART | Not applicable |

### Transmission

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `beta_m2f` | 0.05 | Per-act male-to-female transmission probability; usually calibrated | Not recorded |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not recorded |
| `beta_m2c` | 0.025 | Prenatal mother-to-child transmission rate, per month | Not recorded |
| `beta_breastfeed` | 0.005 | Postnatal (breastfeeding) transmission rate, per month | Not recorded |
| `beta_m2m` | 0 | Per-act male-to-male transmission probability (used on an `msm` network) | Not applicable |
| `rel_trans_acute` | normal(5.3, 0.5) | Relative transmissibility during the acute phase | [Bellan 2015](https://doi.org/10.1371/journal.pmed.1001801) |
| `rel_trans_falling` | normal(8, 0.5) | Relative transmissibility during late-stage CD4 decline | Not recorded |
| `eff_condom` | 0.9 | Condom efficacy against transmission | Not recorded |
| `eff_circ` | 0.6 | Reduction in male susceptibility from program (medical) circumcision | Not recorded |
| `eff_circ_traditional` | 0.6 | Reduction in male susceptibility from traditional circumcision | Not recorded |
| `rel_sus_age` | None | Age- and sex-stratified susceptibility multipliers | Not applicable |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.05 | Initial prevalence | Not recorded |
| `rel_init_prev` | 1 | Multiplier on initial prevalence | Not applicable |
| `init_diagnosed` | 0.0 | Proportion of initial infections already diagnosed | Not recorded |
| `dist_ti_init_infected` | uniform(-120, -5) | Time of initial infection relative to the start (timesteps; months at the default `dt`) | Not recorded |

### Treatment (ART)

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `effective_art_efficacy` | 0.99 | Reduction in transmission on effective (virally suppressive) ART | [Rodger 2016](https://doi.org/10.1001/jama.2016.5148); [Rodger 2019](https://doi.org/10.1016/S0140-6736(19)30418-0); [Bavinton 2018](https://doi.org/10.1016/S2352-3018(18)30132-2) |
| `nonsupp_art_efficacy` | 0.35 | Reduction in transmission on non-suppressive ART (see "On-ART transmission" above) | Derived from [Quinn 2000](https://doi.org/10.1056/NEJM200003303421303) |
| `time_to_art_efficacy` | 6 months | Time to full ART efficacy (linear ramp) | Not recorded |
| `p_effective_art` | bernoulli(1.0) | Probability that ART initiation achieves viral suppression | Not recorded |
| `art_cd4_growth` | 0.1 | Logistic growth rate of CD4 reconstitution on ART | Fitted to the [EMOD](https://docs.idmod.org/projects/emod-hiv/en/latest/hiv-model-healthcare-systems.html) reconstitution curve; [Seoane 2008](https://doi.org/10.1186/1471-2334-8-20); [Darraj 2018](https://doi.org/10.1016/j.jiph.2017.08.001) |
| `art_cd4_pars` | cd4_max=1000, cd4_healthy=500 | Upper bound and healthy reference CD4 count for reconstitution | Not recorded |
| `dur_on_art` | lognorm(3 yr, 1.5 yr) | Duration on ART before discontinuation | Not recorded |
| `rel_dur_on_art` | 1.0 | Multiplier on `dur_on_art`; usually calibrated | Not applicable |
| `dur_on_art_trend` | None | Time-varying multiplier on `dur_on_art` (years and values) | Not applicable |
| `rate_art_dropout` | 0/yr | Rate of stopping ART early (e.g. from poor adherence), besides at the end of `dur_on_art` | Not applicable |
| `dur_post_art_scale_factor` | 0.1 | Standard deviation of post-ART CD4 decline duration, as a fraction of its mean | Not recorded |
| `rel_art_mortality_effective` | 0.25 | Fraction of off-ART CD4-based mortality retained on effective ART, both sexes (see "On-ART mortality" above) | Not recorded |
| `rel_art_mortality_unsupp_m` | 0.7 | Fraction of off-ART CD4-based mortality retained on non-suppressive ART, males | Not recorded |
| `rel_art_mortality_unsupp_f` | 0.35 | Fraction of off-ART CD4-based mortality retained on non-suppressive ART, females; gives a non-suppressive/effective mortality ratio of 2.8 (males) vs 1.4 (females) | Not recorded |
| `rel_death_f` | 0.74 | Mortality multiplier for females, on and off ART | Not recorded |
| `art_death_age` | (0, 25, 1.0), (25, 35, 1.10), (35, 45, 1.21), (45, 125, 1.32) | On-ART mortality multipliers by age, as `(age_lo, age_hi, mult)` | Not recorded |
| `art_death_dur` | None | Reserved for on-ART mortality multipliers by time since initiation; currently unused | Not applicable |

### Care seeking

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `care_seeking` | normal(1, 0.5) | Relative care-seeking propensity (per agent) | Not recorded |
| `maternal_care_scale` | 2 | Multiplier on care seeking during pregnancy | Not recorded |

## Results

HIV produces STIsim's standard results (`new_infections`, `prevalence`, `incidence`, etc.), one for each state and transition (e.g. `acute`, `art_dropout`), plus HIV-specific results:

### HIV-specific results

| Result | Description |
|--------|-------------|
| `new_deaths` | HIV/AIDS deaths this timestep |
| `cum_deaths` | Cumulative HIV/AIDS deaths |
| `new_diagnoses` | Newly diagnosed this timestep |
| `cum_diagnoses` | Cumulative diagnoses |
| `new_agents_on_art` | Agents starting ART this timestep |
| `p_on_art` | Proportion of infected agents on ART |
| `prevalence_15_49` | HIV prevalence among 15-49 year olds |
| `n_on_art_pregnant` | Number of pregnant women on ART (only when pregnancy is in the sim) |
| `p_diagnosed_pregnant` | Proportion of HIV+ pregnant women who are diagnosed (only when pregnancy is in the sim) |

### Transmission route results

All STI diseases (including HIV) track infections by transmission route:

| Result | Description |
|--------|-------------|
| `new_infections` | Total new infections (sexual + MTCT) |
| `new_infections_sex` | New infections via sexual transmission |
| `new_infections_mtct` | New infections via mother-to-child transmission |

These are always consistent: `new_infections_sex + new_infections_mtct == new_infections`.

When pregnancy is modeled, prenatal and postnatal MTCT are tracked separately:

| Result | Description |
|--------|-------------|
| `new_infections_prenatal` | New infections via prenatal (in utero) transmission |
| `new_infections_postnatal` | New infections via postnatal (breastfeeding) transmission |

These satisfy: `new_infections_prenatal + new_infections_postnatal == new_infections_mtct`.

**Accessing MTCT results:**

```python
sim.run()

# Total MTCT infections over the simulation
total_mtct = sim.results.hiv.new_infections_mtct.sum()

# Time series
plt.plot(sim.t.yearvec, sim.results.hiv.new_infections_mtct)
```

## PMTCT (prevention of mother-to-child transmission)

STIsim models PMTCT through three mechanisms:

### 1. ANC testing

ANC (antenatal care) testing identifies HIV-positive pregnant women so they can start ART. This is implemented as an `HIVTest` with pregnancy-based eligibility, the same pattern used for FSW-targeted testing:

```python
import stisim as sti

# Test undiagnosed pregnant women in first trimester
anc_test = sti.HIVTest(
    test_prob_data=0.9,
    dt_scale=False,
    name='anc_test',
    eligibility=lambda sim: (sim.demographics.pregnancy.P & (sim.demographics.pregnancy.P.time < ss.weeks(13)) # First trimester
                             & ~sim.diseases.hiv.diagnosed),
)
```

ANC testing is not included in `hivsim.Sim` defaults (which use a single general-population `HIVTest`). Add it explicitly when modeling targeted testing pathways. The `p_diagnosed_pregnant` result tracks what proportion of HIV+ pregnant women have been diagnosed, which measures the effectiveness of the ANC testing pathway.

### 2. ART retention during pregnancy

Pregnant women on ART are less likely to drop out thanks to the `maternal_care_scale` parameter (default 2), which doubles care-seeking behavior during pregnancy. This makes it much less likely that pregnant women stop ART, keeping them on treatment through delivery and breastfeeding.

### 3. Prenatal protection (MaternalNet)

When a pregnant woman is on ART, her unborn infant's susceptibility is reduced by the `pmtct_efficacy` parameter on the ART intervention (default 0.96). This is applied each timestep via the MaternalNet.

### 4. Postnatal protection (BreastfeedingNet)

When a breastfeeding mother is on ART, her infant's susceptibility is similarly reduced by `pmtct_efficacy` via the BreastfeedingNet.

For both prenatal and postnatal transmission, total protection compounds two effects: the infant's reduced susceptibility (`rel_sus`, from `pmtct_efficacy`) and the mother's reduced transmissibility (`rel_trans`, from `effective_art_efficacy`/`nonsupp_art_efficacy`).

### Configuring PMTCT

```python
import stisim as sti

# Custom PMTCT efficacy
art = sti.ART(pmtct_efficacy=0.98)

# Complete protection (previous default behavior)
art = sti.ART(pmtct_efficacy=1.0)
```

## hivsim.Sim parameter routing

`hivsim.Sim` is a thin wrapper around `sti.Sim` that supplies HIV-appropriate defaults. It accepts parameters in several forms, all routed by `sti.Sim.separate_pars`:

| Kwarg | What it contains | Example |
|---|---|---|
| `pars` | Any flat parameter dict — entries are routed to the right module automatically | `pars=dict(n_agents=500, beta_m2f=0.03)` |
| `sim_pars` | Parameters specifically for the base `sti.Sim` (start/stop/n_agents/dt/…) | `sim_pars=dict(start=2000, n_agents=500)` |
| `hiv_pars` | Explicit HIV module pars (legacy; prefer flat kwargs or `hiv=dict(...)`) | `hiv_pars=dict(beta_m2f=0.03)` |
| flat kwargs | Any recognised par name, routed automatically | `hivsim.Sim(beta_m2f=0.03, debut_f=18)` |
| `hiv=dict(...)` | Disease-keyed dict — applies only to the HIV module | `hiv=dict(rel_trans_acute=10)` |

In most cases you only need flat kwargs:

```python
import hivsim

# Override HIV pars
sim = hivsim.Sim(beta_m2f=0.03, rel_trans_acute=10)

# Override network pars
sim = hivsim.Sim(debut_f=18, fsw_shares=0.03)

# Mix
sim = hivsim.Sim(beta_m2f=0.03, n_agents=5000, start=1990)
```

The same patterns work via `hivsim.demo`:

```python
sim = hivsim.demo('simple', run=False, beta_m2f=0.03)
sim = hivsim.demo('zimbabwe', run=False, hiv=dict(rel_trans_acute=10))
```

**Rule**: pass a module instance (e.g. `diseases=sti.HIV(...)`) **or** pars for the default — not both. `sti.Sim` raises an error if you do both for the same module slot.
