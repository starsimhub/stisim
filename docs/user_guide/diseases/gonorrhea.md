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

### Natural history

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dur_exp` | 0 | No latent period (immediately infectious) |
| `p_symp` | [0.35, 0.65] | Probability of symptoms [F, M] |
| `dur_presymp` | F: [1 wk, 12 wk]; M: [0.25 wk, 1 wk] | Time to symptom onset [mean, std] |
| `dur_asymp2clear` | F: [8 mo, 2 mo]; M: [6 mo, 3 mo] | Duration of asymptomatic infection [mean, std] |
| `dur_symp2clear` | F: [9 mo, 2 mo]; M: [6 mo, 3 mo] | Duration of symptomatic infection [mean, std] |

### Care seeking

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_symp_care` | [0.66, 0.83] | Probability of seeking care if symptomatic [F, M] |
| `dur_symp2care` | F: [2 mo, 1 mo]; M: [1 wk, 2 wk] | Time from symptoms to care seeking [mean, std] |

### Complications

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_pid` | 0.0 | Probability of pelvic inflammatory disease (females) |
| `dur_prepid` | lognorm(1.5 mo, 3 mo) | Time to PID onset |

### Transmission

| Parameter | Default | Description |
|-----------|---------|-------------|
| `beta_m2f` | 0.06 | Per-act male-to-female transmission probability |
| `rel_beta_f2m` | 0.5 | Female-to-male transmission relative to male-to-female |
| `beta_m2c` | None | Mother-to-child transmission probability |
| `beta_m2m` | None | Male-to-male transmission probability |
| `eff_condom` | 0.9 | Condom efficacy |
| `init_prev` | 0.01 | Initial prevalence |
