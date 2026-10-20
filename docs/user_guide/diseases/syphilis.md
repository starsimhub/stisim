# Syphilis

Syphilis in STIsim progresses through exposed, primary, secondary, latent, and tertiary stages, with stage-dependent transmissibility and congenital syphilis outcomes. Treatment at any stage clears the infection but leaves the agent seropositive and susceptible to reinfection.

**Class:** `sti.Syphilis` | **Alias:** `'syphilis'` / `'syph'` | **Base class:** `BaseSTI`

## States and transitions

```
                         ┌─────────────┐
                         │ Susceptible │◄──────────────────────────────┐
                         └──────┬──────┘                               │
                                │ infection                            │
                                ▼                                      │
                         ┌─────────────┐                               │
                         │   Exposed   │  Incubation                   │
                         │  (~50 d)    │  Not yet infectious           │
                         └──────┬──────┘                               │
                                │                                      │
                                ▼                                      │
                         ┌─────────────┐                               │
                         │   Primary   │  Chancre stage        treatment at
                         │  (~6 wk)    │  Transmissibility: 1x   any stage
                         └──────┬──────┘                         clears infection
                                │                                (sets susceptible
                                ▼                                 + ever_exposed)
                         ┌─────────────┐                               │
                         │  Secondary  │  Systemic symptoms            │
                         │ (~3.6 mo)   │  Transmissibility: 1x        │
                         └──────┬──────┘                               │
                                │                                      │
                                ▼                                      │
                  ┌─────────────────────────────┐                      │
                  │        Latent               │  Transmissibility:   │
                  │   Early (~12-14 mo)         │  decays with         │
                  │   Late (indefinite)         │  half-life of 1 yr   │
                  └──────────┬──────────────────┘                      │
                             │                                         │
                 ┌───────────┼───────────┐                             │
                 │           │           │                             │
                 ▼           ▼           ▼                             │
          ┌──────────┐ ┌──────────┐ ┌───────────┐                     │
          │Reactivate│ │ Tertiary │ │ Remains   │                     │
          │   (35%)  │ │  (35%)   │ │  latent   │                     │
          │→Secondary│ │ (~20 yr) │ │  (30%)    │                     │
          └──────────┘ └────┬─────┘ └───────────┘                     │
                            │                                         │
                            ▼
                     ┌─────────────┐
                     │ Death (5%)  │
                     └─────────────┘
```

**Treatment and reinfection:** Treatment (via `SyphTx`) clears all stage states and returns the agent to susceptible. The `ever_exposed` flag persists -- these agents are seropositive (detectable by treponemal tests) but can be reinfected. This reinfection cycle is a key driver of syphilis epidemics: symptomatic agents seek care, get treated, become susceptible again, and are reinfected by their partners.

**Key properties:**

- `active` = primary | secondary (used by connectors for HIV susceptibility effects)
- `infectious` = active | latent
- `naive` = never exposed (seronegative)
- `sus_not_naive` = susceptible & ever_exposed (seropositive but susceptible to reinfection)

## Parameters

Pars of `sti.Syphilis`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Natural history

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur_exposed` | normal(50 d, 10 d) | Duration of incubation period | [Garnett 1997](https://doi.org/10.1097/00007435-199704000-00002) |
| `dur_primary` | normal(6 wk, 1 wk) | Duration of primary stage | [Garnett 1997](https://doi.org/10.1097/00007435-199704000-00002) |
| `dur_secondary` | lognorm(3.6 mo, 1.5 mo) | Duration of secondary stage | [Garnett 1997](https://doi.org/10.1097/00007435-199704000-00002) |
| `dur_early` | uniform(22 mo, 24 mo) | Duration of the early latent stage | WHO Europe definition of early latent syphilis (≤2 years after infection) |
| `p_reactivate` | 0.35 | Probability of reactivation from latent to secondary | [O'Byrne 2019](https://doi.org/10.1136/bmj.l4159) |
| `time_to_reactivate` | lognorm(1 yr, 1 yr) | Time from latent onset to reactivation | [O'Byrne 2019](https://doi.org/10.1136/bmj.l4159) |
| `p_tertiary` | 0.35 | Probability of progressing from latent to tertiary | [Wolujewicz 2016](https://doi.org/10.3399/bjgp16X686065); [O'Byrne 2019](https://doi.org/10.1136/bmj.l4159) |
| `time_to_tertiary` | normal(20 yr, 2 yr) | Time from latent onset to tertiary | [O'Byrne 2019](https://doi.org/10.1136/bmj.l4159) |
| `p_death` | 0.05 | Probability of death from tertiary syphilis | Assumption |
| `time_to_death` | lognorm(5 yr, 5 yr) | Time from tertiary onset to death | Assumption |

### Symptoms and serology

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_symp_primary` | [0.3, 0.5] | Probability that a primary chancre is noticed [F, M] | [O'Byrne 2019](https://doi.org/10.1136/bmj.l4159) |
| `p_symp_secondary` | [0.85, 0.85] | Probability that secondary symptoms are noticed [F, M] | Not recorded |
| `time_to_undetectable` | 15 yr | Mean time in late latent until the non-treponemal titre falls below detection | Not recorded |
| `p_trep_persists` | 0.8 | Probability that a positive treponemal test persists for life | Not recorded |
| `dur_nontrep_revert` | 9 mo | Mean time after early treatment until the non-treponemal (RPR) titre falls below detection | Not recorded |
| `time_to_trep_revert` | normal(10 yr, 5 yr) | Time to treponemal sero-reversion, for those who revert | Not recorded |

### Transmission

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `beta_m2f` | 0.1 | Per-act male-to-female transmission probability; usually calibrated | Not recorded |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not recorded |
| `beta_m2c` | 0 | Monthly rate of in-utero infection of a fetus whose mother is in primary, secondary, or latent syphilis | Not applicable |
| `beta_m2m` | 0 | Per-act male-to-male transmission probability (used on an `msm` network) | Not applicable |
| `eff_condom` | 0.0 | Condom efficacy against transmission | Not recorded |
| `rel_trans_primary` | 1 | Relative transmissibility during the primary stage | Assumption |
| `rel_trans_secondary` | 1 | Relative transmissibility during the secondary stage (reference) | Assumption |
| `rel_trans_tertiary` | 0.0 | Relative transmissibility during the tertiary stage | Assumption |
| `rel_trans_latent_half_life` | 1 yr | Half-life of latent-stage transmissibility, which decays from the secondary-stage value | Assumption |

### Congenital syphilis

Birth outcomes depend on the maternal stage of infection. Probabilities sum to 1 across the five outcomes:

| Maternal stage | Miscarriage | Neonatal death | Stillborn | Congenital syphilis | Normal birth |
|---------------|-------------|----------------|-----------|--------------------|----|
| Active (primary/secondary) | 0.00 | 0.25 | 0.25 | 0.45 | 0.05 |
| Early latent | 0.00 | 0.05 | 0.10 | 0.40 | 0.45 |
| Late latent | 0.00 | 0.00 | 0.10 | 0.10 | 0.80 |

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_miscarriage_active` | 0.00 | Probability of miscarriage, for an active (primary or secondary) maternal infection | [Gomez 2013](https://doi.org/10.2471/BLT.12.107623); [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_nnd_active` | 0.25 | Probability of neonatal death, for an active (primary or secondary) maternal infection | [Gomez 2013](https://doi.org/10.2471/BLT.12.107623); [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_stillborn_active` | 0.25 | Probability of stillbirth, for an active (primary or secondary) maternal infection | [Gomez 2013](https://doi.org/10.2471/BLT.12.107623); [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_congenital_active` | 0.45 | Probability of congenital syphilis, for an active (primary or secondary) maternal infection | [Gomez 2013](https://doi.org/10.2471/BLT.12.107623); [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_miscarriage_early` | 0.00 | Probability of miscarriage, for an early latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_nnd_early` | 0.05 | Probability of neonatal death, for an early latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_stillborn_early` | 0.10 | Probability of stillbirth, for an early latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_congenital_early` | 0.40 | Probability of congenital syphilis, for an early latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_miscarriage_late` | 0.00 | Probability of miscarriage, for a late latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_nnd_late` | 0.00 | Probability of neonatal death, for a late latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_stillborn_late` | 0.10 | Probability of stillbirth, for a late latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `p_congenital_late` | 0.10 | Probability of congenital syphilis, for a late latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.0 | Initial active syphilis prevalence | Not applicable |
| `init_latent_prev` | 0.0 | Initial latent syphilis prevalence | Not applicable |
