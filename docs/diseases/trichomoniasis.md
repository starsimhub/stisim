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

### Natural history

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dur_exp` | 0 | No latent period |
| `p_symp` | [0.40, 0.50] | Probability of symptoms [F, M] |
| `dur_presymp` | F: [1 wk, 12 wk]; M: [0.25 wk, 1 wk] | Time to symptom onset [mean, std] |
| `dur_asymp2clear` | F: [48 mo, 6 mo]; M: [26 wk, 4 wk] | Duration of asymptomatic infection [mean, std] |
| `dur_symp2clear` | F: [20 wk, 4 wk]; M: [18 wk, 4 wk] | Duration of symptomatic infection [mean, std] |
| `p_clear` | 0.10 | Probability of spontaneous clearance in women |
| `dur_persist` | 100 yr | Duration of persistent infection (effectively lifelong) |

### Care seeking

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_symp_care` | [0.39, 0.27] | Probability of seeking care if symptomatic [F, M] |
| `dur_symp2care` | F: [2 mo, 1 mo]; M: [1 wk, 2 wk] | Time from symptoms to care seeking [mean, std] |

### Transmission

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eff_condom` | 0.0 | Condom efficacy |
| `init_prev` | 0.01 | Initial prevalence |
