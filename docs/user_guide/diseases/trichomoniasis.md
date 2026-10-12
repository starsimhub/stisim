# Trichomoniasis

Trichomoniasis (*Trichomonas vaginalis*) follows the SEIS pattern with a key feature: most women do not spontaneously clear the infection, leading to very long durations.

**Class:** `sti.Trichomoniasis` | **Alias:** `'tv'` | **Base class:** `SEIS`

## States and transitions

```
     ┌─────────────┐
     │ Susceptible │◄──────────────────────────────────────┐
     └──────┬──────┘                                       │
            │ infection (immediately infectious)            │
            │                                               │
            ├────────────────────┐                          │
            ▼                    ▼                          │
     ┌─────────────┐     ┌──────────────┐                  │
     │Asymptomatic │     │ Symptomatic  │                  │
     │  F: 60%     │     │  F: 40%      │                  │
     │  M: 50%     │     │  M: 50%      │                  │
     └──────┬──────┘     └──────┬───────┘                  │
            │                   │                           │
            │                   ▼                           │
            │            ┌───────────┐                      │
            │            │ Seeks care│                      │
            │            └───────────┘                      │
            │                                               │
            ├── 10% of F: natural clearance ───────────────┤
            │   (F: ~48 mo, M: ~26 wk)                     │
            │                                               │
            └── 90% of F: persistent infection              │
                M: natural clearance ──────────────────────┘
```

## Parameters

Pars of `sti.Trichomoniasis`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Natural history

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur_exp` | 0 | Duration of the exposed (non-infectious) period | Not recorded |
| `p_symp` | [0.40, 0.50] | Probability of developing symptoms [F, M] | [Bowden 2000](https://doi.org/10.1136/sti.76.4.248) |
| `dur_presymp` | F: [1 wk, 12 wk]; M: [0.25 wk, 1 wk] | Time from infection to symptom onset, lognormal [mean, sd] | Not recorded |
| `dur_asymp2clear` | F: [48 mo, 6 mo]; M: [26 wk, 4 wk] | Duration of untreated asymptomatic infection, lognormal [mean, sd] | [Bowden 2000](https://doi.org/10.1136/sti.76.4.248) |
| `dur_symp2clear` | F: [20 wk, 4 wk]; M: [18 wk, 4 wk] | Duration of untreated symptomatic infection, lognormal [mean, sd] | Assumption |
| `p_symp_clear` | [0.0, 0.0] | Probability that symptoms resolve before clearance [F, M] | Not applicable |
| `dur_symp` | lognorm(1 wk, 26 wk) [F, M] | Duration of symptoms, for those whose symptoms resolve | Not recorded |
| `p_clear` | 0.10 | Probability that an asymptomatic infection in a woman clears spontaneously; otherwise it persists | [Bowden 2000](https://doi.org/10.1136/sti.76.4.248) |
| `dur_persist` | 100 yr | Duration of persistent infection (effectively lifelong) | Not recorded |

### Care seeking

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_symp_care` | [0.39, 0.27] | Probability of seeking care if symptomatic [F, M] | Not recorded |
| `dur_symp2care` | F: [2 mo, 1 mo]; M: [1 wk, 2 wk] | Time from symptom onset to care seeking, lognormal [mean, sd] | Not recorded |

### Complications

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_pid` | 0.0 | Probability of pelvic inflammatory disease (females) | Not recorded |
| `dur_prepid` | lognorm(6 wk, 4 wk) | Time from infection to PID onset | Not recorded |
| `p_pid_care` | 0.1 | Probability of seeking care for PID | Not recorded |
| `dur_pid2care` | lognorm(2 wk, 4 wk) | Time from PID onset to care seeking | Not recorded |
| `dur_pid2clear` | lognorm(52 wk, 5 wk) | Duration of untreated PID | Not recorded |

### Transmission

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `beta_m2f` | 0.1 | Per-act male-to-female transmission probability; usually calibrated | Not recorded |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not recorded |
| `beta_m2c` | None | Mother-to-child transmission probability per timestep | Not applicable |
| `beta_breastfeed` | None | Breastfeeding transmission probability per timestep | Not applicable |
| `beta_m2m` | None | Per-act male-to-male transmission probability (used on an `msm` network) | Not applicable |
| `eff_condom` | 0.0 | Condom efficacy against transmission | Not recorded |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.01 | Initial prevalence | Not recorded |
| `rel_init_prev` | 1 | Multiplier on initial prevalence | Not applicable |
