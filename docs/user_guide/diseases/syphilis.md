# Syphilis

Syphilis in STIsim progresses through primary, secondary, latent, and tertiary stages, with stage-dependent transmissibility and congenital syphilis outcomes.

**Class:** `sti.Syphilis` | **Alias:** `'syphilis'` / `'syph'` | **Base class:** `BaseSTI`

## States and transitions

```
                         ┌─────────────┐
                         │ Susceptible │
                         └──────┬──────┘
                                │ infection
                                ▼
                         ┌─────────────┐
                         │   Primary   │  Chancre stage
                         │  (~6 wk)   │  Transmissibility: 1x
                         └──────┬──────┘
                                │
                                ▼
                         ┌─────────────┐
                         │  Secondary  │  Systemic symptoms
                         │ (~3.6 mo)  │  Transmissibility: 1x
                         └──────┬──────┘
                                │
                                ▼
                  ┌─────────────────────────────┐
                  │        Latent               │
                  │   Early (~12-14 mo)         │  Transmissibility: decays
                  │   Late (indefinite)         │  with half-life of 1 yr
                  └──────────┬──────────────────┘
                             │
                 ┌───────────┼───────────┐
                 │           │           │
                 ▼           ▼           ▼
          ┌──────────┐ ┌──────────┐ ┌───────────┐
          │Reactivate│ │ Tertiary │ │ Remains   │
          │   (35%)  │ │  (35%)   │ │  latent   │
          │→Secondary│ │ (~20 yr) │ │           │
          └──────────┘ └────┬─────┘ └───────────┘
                            │
                            ▼
                     ┌─────────────┐
                     │ Death (5%)  │
                     └─────────────┘
```

## Parameters

### Natural history

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dur_primary` | normal(6 wk, 1 wk) | Duration of primary stage |
| `dur_secondary` | lognorm(3.6 mo, 1.5 mo) | Duration of secondary stage |
| `dur_early` | uniform(12 mo, 14 mo) | Duration of early latent stage |
| `p_reactivate` | 0.35 | Probability of reactivation from latent to secondary |
| `time_to_reactivate` | lognorm(1 yr, 1 yr) | Time to reactivation |
| `p_tertiary` | 0.35 | Probability of progressing to tertiary |
| `time_to_tertiary` | normal(20 yr, 2 yr) | Time to tertiary stage |
| `p_death` | 0.05 | Probability of death from tertiary syphilis |
| `time_to_death` | lognorm(5 yr, 5 yr) | Time to death after tertiary |

### Transmission

| Parameter | Default | Description |
|-----------|---------|-------------|
| `beta_m2f` | 0.1 | Male-to-female transmission probability |
| `eff_condom` | 0.0 | Condom efficacy (syphilis transmits via skin contact) |
| `rel_trans_primary` | 1 | Relative transmissibility during primary stage |
| `rel_trans_secondary` | 1 | Relative transmissibility during secondary stage |
| `rel_trans_latent` | 1 | Baseline latent transmissibility (decays exponentially) |
| `rel_trans_tertiary` | 0.0 | Relative transmissibility during tertiary (non-infectious) |
| `rel_trans_latent_half_life` | 1 yr | Half-life of latent-stage transmissibility decay |

### Congenital syphilis

Birth outcomes depend on maternal stage of infection (probabilities sum to 1 across 5 outcomes):

| Maternal stage | Miscarriage | Neonatal death | Stillborn | Congenital syphilis | Normal birth |
|---------------|-------------|----------------|-----------|--------------------|----|
| Active (primary/secondary) | 0.00 | 0.10 | 0.20 | 0.45 | 0.25 |
| Early latent | 0.00 | 0.05 | 0.10 | 0.40 | 0.45 |
| Late latent | 0.00 | 0.00 | 0.10 | 0.10 | 0.80 |

### Initialization

| Parameter | Default | Description |
|-----------|---------|-------------|
| `init_prev` | 0.0 | Initial active syphilis prevalence |
| `init_latent_prev` | 0.0 | Initial latent syphilis prevalence |
