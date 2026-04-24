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
                        │  (~3 mo)   │  Transmissibility: 6x baseline
                        └──────┬──────┘
                               │
                               ▼
                        ┌─────────────┐
                        │   Latent    │  CD4 stable at ~500
                        │  (~10 yr)  │  Transmissibility: 1x baseline
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
