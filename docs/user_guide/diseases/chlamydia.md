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

### Natural history

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dur_exp` | 1 wk | Duration of exposed (latent) period |
| `p_symp` | [0.20, 0.54] | Probability of symptoms [F, M] |
| `dur_presymp` | F: [1 wk, 10 wk]; M: [0.25 wk, 1 wk] | Time to symptom onset [mean, std] |
| `dur_asymp2clear` | F: [18 mo, 1 mo]; M: [12 mo, 1 mo] | Duration of asymptomatic infection [mean, std] |
| `dur_symp2clear` | F: [18 mo, 1 mo]; M: [12 mo, 1 mo] | Duration of symptomatic infection [mean, std] |

### Care seeking

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_symp_care` | [0.42, 0.83] | Probability of seeking care if symptomatic [F, M] |
| `dur_symp2care` | F: [2 mo, 1 mo]; M: [1 wk, 2 wk] | Time from symptoms to care seeking [mean, std] |

### Complications

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_pid` | 0.0 | Probability of pelvic inflammatory disease (females) |
| `dur_prepid` | lognorm(1.5 mo, 3 mo) | Time to PID onset |

### Transmission

| Parameter | Default | Description |
|-----------|---------|-------------|
| `beta_m2f` | None | Male-to-female transmission probability (set by network) |
| `rel_beta_f2m` | 0.5 | Female-to-male transmission relative to male-to-female |
| `beta_m2c` | None | Mother-to-child transmission probability |
| `beta_m2m` | None | Male-to-male transmission probability |
| `eff_condom` | 0.0 | Condom efficacy |
| `init_prev` | 0.01 | Initial prevalence |
