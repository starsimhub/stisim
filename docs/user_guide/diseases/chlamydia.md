# Chlamydia

Chlamydia (*Chlamydia trachomatis*) follows the SEIS (susceptible-exposed-infectious-susceptible) pattern with sex-stratified symptom probabilities and natural clearance.

**Class:** `sti.Chlamydia` | **Alias:** `'ct'` | **Base class:** `SEIS`

## States and transitions

```
     ┌─────────────┐
     │ Susceptible │◄──────────────────────────────────────┐
     └──────┬──────┘                                       │
            │ infection                                     │
            ▼                                               │
     ┌─────────────┐                                       │
     │   Exposed   │  Latent period                        │
     │   (~1 wk)   │  (not yet infectious)                 │
     └──────┬──────┘                                       │
            │                                               │
            ├────────────────────┐                          │
            ▼                    ▼                          │
     ┌─────────────┐     ┌──────────────┐                  │
     │Asymptomatic │     │ Symptomatic  │                  │
     │  F: 80%     │     │  F: 20%      │                  │
     │  M: 46%     │     │  M: 54%      │                  │
     └──────┬──────┘     └──────┬───────┘                  │
            │                   │                           │
            │                   ├───────────┐               │
            │                   ▼           ▼               │
            │            ┌───────────┐  ┌───────┐          │
            │            │ Seeks care│  │  PID  │ (F only) │
            │            └───────────┘  └───────┘          │
            │                                               │
            └───────── natural clearance ──────────────────┘
                    F: ~18 mo   M: ~12 mo
```

## Parameters

Pars of `sti.Chlamydia`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Natural history

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur_exp` | 1 wk | Duration of the exposed (non-infectious) period | Not recorded |
| `p_symp` | [0.20, 0.54] | Probability of developing symptoms [F, M] | [Althaus 2010](https://doi.org/10.1016/j.epidem.2010.04.002) |
| `dur_presymp` | F: [1 wk, 10 wk]; M: [0.25 wk, 1 wk] | Time from infection to symptom onset, lognormal [mean, sd] | Not recorded |
| `dur_asymp2clear` | F: [18 mo, 1 mo]; M: [12 mo, 1 mo] | Duration of untreated asymptomatic infection, lognormal [mean, sd] | [Althaus 2010](https://doi.org/10.1016/j.epidem.2010.04.002) |
| `dur_symp2clear` | F: [18 mo, 1 mo]; M: [12 mo, 1 mo] | Duration of untreated symptomatic infection, lognormal [mean, sd] | Assumption |
| `p_clear` | 1.0 | Probability that an asymptomatic infection in a woman clears spontaneously; otherwise it persists | Not applicable |
| `dur_persist` | 100 yr | Duration of a persistent infection | Not applicable |

### Care seeking

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_symp_care` | [0.42, 0.83] | Probability of seeking care if symptomatic [F, M] | Unpublished evidence review |
| `dur_symp2care` | F: [2 mo, 1 mo]; M: [1 wk, 2 wk] | Time from symptom onset to care seeking, lognormal [mean, sd] | Not recorded |

### Complications

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_pid` | 0.0 | Probability of pelvic inflammatory disease (females) | Assumption, following [Gray 2009](https://doi.org/10.1086/598983) and [Wiesenfeld 2002](https://doi.org/10.1016/s0029-7844(02)02118-x) |
| `dur_prepid` | lognorm(1.5 mo, 3 mo) | Time from infection to PID onset | Not recorded |
| `p_pid_care` | 0.1 | Probability of seeking care for PID | Not recorded |
| `dur_pid2care` | lognorm(2 wk, 4 wk) | Time from PID onset to care seeking | Not recorded |
| `dur_pid2clear` | lognorm(52 wk, 5 wk) | Duration of untreated PID | Not recorded |

### Transmission

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `beta_m2f` | 0.06 | Per-act male-to-female transmission probability; usually calibrated | Not recorded |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not recorded |
| `beta_m2c` | None | Mother-to-child transmission probability per timestep | Not applicable |
| `beta_m2m` | None | Per-act male-to-male transmission probability (used on an `msm` network) | Not applicable |
| `eff_condom` | 0.4 | Condom efficacy against transmission | [Paz-Bailey 2005](https://doi.org/10.1001/archpedi.159.6.536) |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.01 | Initial prevalence | Not recorded |
