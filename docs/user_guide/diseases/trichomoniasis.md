# Trichomoniasis

Trichomoniasis (*Trichomonas vaginalis*) follows the SEIS pattern with a key feature: most women do not spontaneously clear the infection, leading to very long durations.

**Class:** `sti.Trichomoniasis` | **Alias:** `'tv'` | **Base class:** `sti.SEIS` (declared: its states and transitions)

## States and transitions

```
susceptible ──infection──▶ exposed ──onset──▶ presymptomatic ──symptoms──▶ symptomatic ──clear_symp──▶ susceptible
             (latent period: none)  │    (F 40% / M 50%)
                                     └──▶ asymptomatic ──clear_asymp (or, for 90% of women, persist)──▶ susceptible

symptomatic ──seek_care──▶ symptomatic                   (a share seek care, and stay symptomatic)
presymptomatic, asymptomatic, symptomatic ──pid_onset──▶ pid ──clear_pid──▶ susceptible   (women)
```

Women may also develop pelvic inflammatory disease (`pid`, with probability `p_pid`), which clears after `dur_pid2clear`; a share of them seek care (`seek_care_pid`). `seek_care` and `seek_care_pid` are events that interventions can be triggered by (e.g. `who=ng.seek_care` for an `ss.treat`, or `eligibility=lambda sim: sim.diseases.ng.seek_care.uids` for an STIsim treatment); `infected` and `infectious` are conditions over the states.

Values that differ by sex are given with `sti.by_sex(female, male)`, and lognormal durations that do with `sti.lognorm_by_sex([mean, sd], [mean, sd])`.

## Parameters

Pars of `sti.Trichomoniasis`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Natural history

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur_exp` | constant(0) | Duration of the exposed (non-infectious) period: none | Not recorded |
| `p_symp` | by_sex(0.40, 0.50) | Probability of developing symptoms | [Bowden 2000](https://doi.org/10.1136/sti.76.4.248) |
| `dur_presymp` | F: [1 wk, 12 wk]; M: [0.25 wk, 1 wk] | Time from becoming infectious to symptom onset, lognormal [mean, sd] | Not recorded |
| `dur_asymp2clear` | F: [48 mo, 6 mo]; M: [26 wk, 4 wk] | Duration of untreated asymptomatic infection, lognormal [mean, sd] | [Bowden 2000](https://doi.org/10.1136/sti.76.4.248) |
| `dur_symp2clear` | F: [20 wk, 4 wk]; M: [18 wk, 4 wk] | Duration of untreated symptomatic infection from symptom onset, lognormal [mean, sd] | Assumption |
| `p_clear` | 0.10 | Probability that an asymptomatic infection in a woman clears spontaneously; otherwise it persists | [Bowden 2000](https://doi.org/10.1136/sti.76.4.248) |
| `dur_persist` | constant(100 yr) | Duration of a persistent infection (effectively lifelong) | Not recorded |

### Care seeking

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_symp_care` | by_sex(0.39, 0.27) | Probability of seeking care if symptomatic | Not recorded |
| `dur_symp2care` | F: [2 mo, 1 mo]; M: [1 wk, 2 wk] | Time from symptom onset to care seeking, lognormal [mean, sd] | Not recorded |

### Complications

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_pid` | 0.0 | Probability of pelvic inflammatory disease (females) | Not recorded |
| `dur_prepid` | lognorm(6 wk, 4 wk) | Time from becoming infectious to PID onset | Not recorded |
| `p_pid_care` | 0.1 | Probability of seeking care for PID | Not recorded |
| `dur_pid2care` | lognorm(2 wk, 4 wk) | Time from PID onset to care seeking | Not recorded |
| `dur_pid2clear` | lognorm(52 wk, 5 wk) | Duration of untreated PID | Not recorded |

### Transmission

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `beta_m2f` | 0.1 | Per-act male-to-female transmission probability; usually calibrated | Not recorded |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not recorded |
| `beta_m2c` | 0 | Mother-to-child transmission rate in pregnancy, per month | Not applicable |
| `beta_m2m` | 0 | Per-act male-to-male transmission probability (used on an `msm` network) | Not applicable |
| `eff_condom` | 0.0 | Condom efficacy against transmission | Not recorded |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.01 | Initial prevalence, seeded into `exposed` | Not recorded |
