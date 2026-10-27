# Syphilis

Syphilis in STIsim progresses through exposed, primary, secondary, latent, and tertiary stages, with stage-dependent transmissibility and congenital syphilis outcomes. Treatment at any stage clears the infection but leaves the agent seropositive and susceptible to reinfection.

**Class:** `sti.Syphilis` | **Alias:** `'syphilis'` / `'syph'` | **Base class:** `ss.Infection`

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

The stages are states (`exposed`, `primary`, `secondary`, `early`, `late`, `tertiary`). The time each one ends is drawn when it starts: at infection, the incubation and primary stage; on entering the secondary stage, its length; and on entering the latent stage, the time to late latent, and a relapse to secondary or tertiary syphilis (and death) to come, timed from the start of the latent stage, so moving from early to late latent does not restart them. A fetus is infected at a rate `beta_m2c` per month while its mother is infectious, into a state by the outcome to come, which is decided by the mother's stage: `fetus_miscarriage` and `fetus_stillborn` die at the delivery (in the step before it, while unborn), `fetus_nnd` dies once born, `fetus_congenital` is born with congenital syphilis (to `congenital`), and `fetus_normal` is born susceptible.

**Treatment and reinfection:** Treatment (via `SyphTx`) moves the agent back to `susceptible` from any stage. The `ever_exposed` flag persists -- these agents are seropositive (detectable by treponemal tests) but can be reinfected. This reinfection cycle is a key driver of syphilis epidemics: symptomatic agents seek care, get treated, become susceptible again, and are reinfected by their partners. Treating a pregnant woman also makes an infected fetus's outcome `fetus_normal` (`syph.cure_fetuses()`).

**Not stages:** serology (`trep`, `nontrep`, and their timers), whether the chancre and rash are visible (`chancre_visible`, `rash_visible`), and the decay of infectiousness in the latent stage are arrays of each person's, updated by code at the events of the stages.

**Key properties:**

- `symptomatic` = primary | secondary (used by connectors for HIV susceptibility effects)
- `latent` = early | late
- `infectious` = symptomatic | latent
- `infected` = exposed | infectious | tertiary
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
| `time_to_undetectable` | lognorm(15 yr, 5 yr) | Time from late-latent onset until the non-treponemal titre falls below detection | Not recorded |
| `dur_to_trep` | normal(3 wk, 1 wk) | Time from infection to a positive treponemal test (window period) | Not recorded |
| `p_trep_persists` | 0.8 | Probability that a positive treponemal test persists for life | Not recorded |
| `time_to_trep_revert` | normal(10 yr, 5 yr) | Time to treponemal sero-reversion, for those who revert | Not recorded |
| `anc_detection` | 0.8 | Proportion of prevalent infections in pregnant women counted in the `detected_pregnant_prevalence` result | Not recorded |

### Transmission

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `beta_m2f` | 0.1 | Per-act male-to-female transmission probability; usually calibrated | Not recorded |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not recorded |
| `beta_m2c` | 0 | Mother-to-fetus transmission rate while the mother is infectious, per month | Not applicable |
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
| `birth_outcomes` (mat_active) | [0.00, 0.25, 0.25, 0.45, 0.05] | Outcome probabilities for an active (exposed, primary, or secondary) maternal infection | [Gomez 2013](https://doi.org/10.2471/BLT.12.107623); [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `birth_outcomes` (early) | [0.00, 0.05, 0.10, 0.40, 0.45] | Outcome probabilities for an early latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |
| `birth_outcomes` (late) | [0.00, 0.00, 0.10, 0.10, 0.80] | Outcome probabilities for a late latent maternal infection | [Stoltey 2015](https://doi.org/10.1071/SH14174); [Arnold 2000](https://doi.org/10.1093/pch/5.8.463) |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.0 | Initial prevalence, seeded into `exposed` | Not applicable |
| `init_latent_prev` | 0.0 | Initial prevalence of early latent syphilis | Not applicable |
