# Gonorrhea

Gonorrhea (*Neisseria gonorrhoeae*) follows the SEIS pattern with immediate infectiousness (no latent period) and tracks drug resistance.

**Class:** `sti.Gonorrhea` | **Alias:** `'ng'` | **Base class:** `SEIS`

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
     │  F: 65%     │     │  F: 35%      │                  │
     │  M: 35%     │     │  M: 65%      │                  │
     └──────┬──────┘     └──────┬───────┘                  │
            │                   │                           │
            │                   ├───────────┐               │
            │                   ▼           ▼               │
            │            ┌───────────┐  ┌───────┐          │
            │            │ Seeks care│  │  PID  │ (F only) │
            │            └───────────┘  └───────┘          │
            │                                               │
            └───────── natural clearance ──────────────────┘
                    F: ~8 mo   M: ~6 mo
```

## Parameters

Pars of `sti.Gonorrhea`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Natural history

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur_exp` | 0 | Duration of the exposed (non-infectious) period | Not recorded |
| `p_symp` | [0.35, 0.65] | Probability of developing symptoms [F, M] | [Newman 2015](https://doi.org/10.1371/journal.pone.0143304) |
| `dur_presymp` | F: [1 wk, 12 wk]; M: [0.25 wk, 1 wk] | Time from infection to symptom onset, lognormal [mean, sd] | Not recorded |
| `dur_asymp2clear` | F: [8 mo, 2 mo]; M: [6 mo, 3 mo] | Duration of untreated asymptomatic infection, lognormal [mean, sd] | Not recorded |
| `dur_symp2clear` | F: [9 mo, 2 mo]; M: [6 mo, 3 mo] | Duration of untreated symptomatic infection, lognormal [mean, sd] | Assumption |
| `p_symp_clear` | [0.0, 0.0] | Probability that symptoms resolve before clearance [F, M] | Not applicable |
| `dur_symp` | lognorm(1 wk, 26 wk) [F, M] | Duration of symptoms, for those whose symptoms resolve | Not recorded |

### Care seeking

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_symp_care` | [0.66, 0.83] | Probability of seeking care if symptomatic [F, M] | Unpublished evidence review |
| `dur_symp2care` | F: [2 mo, 1 mo]; M: [1 wk, 2 wk] | Time from symptom onset to care seeking, lognormal [mean, sd] | Not recorded |

### Complications

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_pid` | 0.0 | Probability of pelvic inflammatory disease (females) | Not recorded |
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
| `beta_breastfeed` | None | Breastfeeding transmission probability per timestep | Not applicable |
| `beta_m2m` | None | Per-act male-to-male transmission probability (used on an `msm` network) | Not applicable |
| `eff_condom` | 0.9 | Condom efficacy against transmission | Not recorded |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.01 | Initial prevalence | Not recorded |
| `rel_init_prev` | 1 | Multiplier on initial prevalence | Not applicable |
