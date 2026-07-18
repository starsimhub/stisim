# HIV

HIV in STIsim is modeled with CD4-based disease progression through acute, latent, and late-stage phases, with ART treatment effects.

**Class:** `sti.HIV` | **Alias:** `'hiv'` | **Base class:** `BaseSTI`

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
     │   (untreated agents only --             │
     │    on-ART agents can also die here if   │
     │    use_art_mortality_table=True, see    │
     │    "ART mortality: two selectable       │
     │    modes" below)                        │
     └─────────────────────────────────────────┘
```

## ART mortality: two selectable modes

This section zooms into a single box from the diagram above: **On ART**. In that diagram, "On ART" has exactly one way out — ART dropout, into Post-ART — because by default agents cannot die from HIV while `on_art = True`; only the off-ART CD4-based hazard (bottom box, "AIDS death") applies, and only once they've dropped out. `HIVPars.use_art_mortality_table` (default `False`) controls whether that's still true, i.e. whether "On ART" gains a **second** way out, straight to AIDS death, alongside the existing dropout path:

```
                        ┌─────────────────┐
                        │     On ART      │  (from the diagram above)
                        └───┬─────────┬───┘
             always: ART    │         │  only if use_art_mortality_table=True:
             dropout        │         │  ARTMortalityTable-style hazard, every timestep
                            ▼         ▼
                  ┌─────────────┐ ┌─────────────┐
                  │  Post-ART   │ │ AIDS death  │
                  └─────────────┘ └─────────────┘
```

Everything else in the diagram above (Falling, Post-ART, the AIDS-death CD4 hazard table) is unaffected by this toggle either way — it only changes whether "On ART" itself is a possible source of death, and if so, how that risk is computed:

```
                     ┌─────────────────────────┐
                     │  On ART (on_art = True) │
                     └───────────┬─────────────┘
                                 │
                      use_art_mortality_table?
              ┌──────────────────┴──────────────────┐
              │                                     │
            False (default)                         True
              │                                     │
              ▼                                     ▼
┌─────────────────────────────┐       ┌────────────────────────────────────┐
│ "Simple" mode               │       │ ARTMortalityTable mode             │
│                             │       │                                    │
│ Mortality = 0 while on ART. │       │ get_art_mortality_hazard():        │
│ Only the off-ART CD4-based  │       │ look up annual hazard by:          │
│ hazard (see diagram above)  │       │   - age          (4 bins)          │
│ ever applies.               │       │   - sex          (m / f)           │
└─────────────────────────────┘       │   - adherence    (effective /      │
                                      │                   non-suppressive) │
                                      │   - duration since ti_art          │
                                      │                   (5 bins)         │
                                      │   - CURRENT cd4  (7 bins)          │
                                      └────────────────────────────────────┘
```

The ARTMortalityTable-mode lookup uses the agent's **current** CD4 count (not frozen at ART initiation) — since STIsim already recomputes CD4 every timestep for the reconstitution curve, there's no reason to freeze it the way EMOD does (see `art_implementation_notes.md` section 9 for the rationale). The default table's numbers (both effective- and non-suppressive-ART variants, by sex) ship from the EMOD-HIV campaign example in `HIVSim_notes/art.md`; pass your own `art_mortality_table` dict of the same shape to use different data. EMOD's "terminal failing window" (a temporary loss of transmission suppression in the months before an ART-mortality-table-scheduled death) is not modeled in either mode.

### Notes on ART

- ART can be initiated from **any** infected stage (Acute, Latent, or Falling)
- On initiation: all pending stage transitions (`ti_latent`, `ti_falling`, `ti_zero`) are cleared for future events
- ART duration is drawn per-agent: `dur_on_art ~ LogNormal(3yr, 1.5yr) × rel_dur_on_art`
  - `rel_dur_on_art` is a calibrated scalar (range 1–20) that stretches mean duration
- After dropout, `ti_zero` is redrawn based on CD4 at dropout and `cd4_potential`
- **No age or sex dependence** in natural history or ART dropout timing. Mortality *while on ART* is age/sex-dependent only if `use_art_mortality_table=True` (see above); it's otherwise zero regardless of age or sex.

### What is not modeled

- Viral load / viral suppression (ART = fully suppressed, binary)
- Age-varying progression rates (natural history is age-independent regardless of ART-mortality mode)
- Age/sex-varying mortality while on ART, unless `use_art_mortality_table=True` (see above)
- Sex-varying natural history
- Re-infection or superinfection
- EMOD's "terminal failing window" (see above)

## Parameters

### Natural history

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cd4_start` | normal(800, 50) | Initial CD4 count at infection |
| `cd4_latent` | normal(500, 50) | CD4 count during latent phase |
| `dur_acute` | lognorm(3 mo, 1 mo) | Duration of acute infection |
| `dur_latent` | lognorm(10 yr, 3 yr) | Duration of latent infection (untreated) |
| `dur_falling` | lognorm(3 yr, 1 yr) | Duration of late-stage CD4 decline |
| `include_aids_deaths` | True | Whether to include AIDS mortality |

### Transmission

| Parameter | Default | Description |
|-----------|---------|-------------|
| `beta_m2f` | 0.05 | Per-act male-to-female transmission probability |
| `rel_beta_f2m` | 0.5 | Female-to-male transmission relative to male-to-female |
| `beta_m2c` | 0.025/mo | Prenatal mother-to-child transmission probability |
| `beta_breastfeed` | 0.005/mo | Postnatal (breastfeeding) transmission probability |
| `rel_trans_acute` | normal(6, 0.5) | Relative transmissibility during acute phase |
| `rel_trans_falling` | normal(8, 0.5) | Relative transmissibility during late stage |
| `eff_condom` | 0.9 | Condom efficacy for reducing transmission |

### Initialization

| Parameter | Default | Description |
|-----------|---------|-------------|
| `init_prev` | 0.05 | Initial prevalence |
| `init_diagnosed` | 0.0 | Proportion initially diagnosed |
| `dist_ti_init_infected` | uniform(-120, -5) | Time of initial infection (months before start) |

### Treatment (ART)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `art_efficacy` | 0.96 | ART efficacy at reducing transmission |
| `time_to_art_efficacy` | 6 months | Time to reach full ART efficacy (linear ramp) |
| `art_cd4_growth` | 0.1 | Logistic growth rate for CD4 reconstitution on ART |
| `dur_on_art` | lognorm(3 yr, 1.5 yr) | Duration on ART before dropout |
| `use_art_mortality_table` | False | If True, apply the ARTMortalityTable-style age/sex/adherence/duration-stratified hazard to on-ART agents instead of zero mortality (see "ART mortality: two selectable modes" above) |
| `art_mortality_table` | EMOD example table | Age/sex/adherence/duration/CD4-stratified annual hazard table used when `use_art_mortality_table=True`; pass your own dict of the same shape to override |

### Care seeking

| Parameter | Default | Description |
|-----------|---------|-------------|
| `care_seeking` | normal(1, 0.5) | Relative care-seeking behavior (per agent) |
| `maternal_care_scale` | 2 | Multiplicative increase in care seeking during pregnancy |

## Results

HIV produces the standard BaseSTI results (`new_infections`, `prevalence`, `incidence`, etc.) plus HIV-specific results:

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
    eligibility=lambda sim: sim.demographics.pregnancy.tri1_uids[
        ~sim.diseases.hiv.diagnosed[sim.demographics.pregnancy.tri1_uids]
    ],
)
```

ANC testing is not included in `hivsim.Sim` defaults (which use a single general-population `HIVTest`). Add it explicitly when modeling targeted testing pathways. The `p_diagnosed_pregnant` result tracks what proportion of HIV+ pregnant women have been diagnosed, which measures the effectiveness of the ANC testing pathway.

### 2. ART retention during pregnancy

Pregnant women on ART are less likely to drop out thanks to the `maternal_care_scale` parameter (default 2), which doubles care-seeking behavior during pregnancy. This makes it much less likely that pregnant women stop ART, keeping them on treatment through delivery and breastfeeding.

### 3. Prenatal protection (MaternalNet)

When a pregnant woman is on ART, her unborn infant's susceptibility is reduced by the `pmtct_efficacy` parameter on the ART intervention (default 0.96). This is applied each timestep via the MaternalNet.

### 4. Postnatal protection (BreastfeedingNet)

When a breastfeeding mother is on ART, her infant's susceptibility is similarly reduced by `pmtct_efficacy` via the BreastfeedingNet.

For both prenatal and postnatal transmission, total protection compounds two effects: the infant's reduced susceptibility (`rel_sus`, from `pmtct_efficacy`) and the mother's reduced transmissibility (`rel_trans`, from `art_efficacy`).

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
