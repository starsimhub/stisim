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
             │   (~3 yr)   │       │  (~18 yr)   │  Transmissibility: 0.04x
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
| `beta_m2f` | 0.05 | Male-to-female transmission probability |
| `rel_beta_f2m` | 0.5 | Female-to-male transmission relative to male-to-female |
| `beta_m2c` | 0.025/mo | Mother-to-child transmission probability |
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
| `dur_on_art` | lognorm(18 yr, 5 yr) | Duration on ART before dropout |

### Care seeking

| Parameter | Default | Description |
|-----------|---------|-------------|
| `care_seeking` | normal(1, 0.1) | Relative care-seeking behavior (per agent) |
| `maternal_care_scale` | 2 | Multiplicative increase in care seeking during pregnancy |
