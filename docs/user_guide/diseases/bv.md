# Bacterial vaginosis (BV)

BV is modeled using community state types (CSTs) representing vaginal microbiome composition. Unlike the other STIs, BV can arise spontaneously (without sexual transmission) and is driven by behavioral risk factors.

**Class:** `sti.BV` | **Alias:** `'bv'` | **Base class:** `BaseSTI`

## States and transitions

```
     ┌─────────────────────────────────────────────┐
     │            Female agents only                │
     │                                              │
     │   ┌───────────┐     ┌───────────┐           │
     │   │   CST 1   │────▶│   CST 3   │           │
     │   │ L. crisp. │◄────│ L. iners  │           │
     │   │  (10%)    │     │  (30%)    │           │
     │   └───────────┘     └─────┬─────┘           │
     │                           │                  │
     │                           ▼                  │
     │                     ┌───────────┐            │
     │                     │   CST 4   │            │
     │                     │   (BV)    │            │
     │                     │  (60%)    │            │
     │                     └─────┬─────┘            │
     │                           │                  │
     │                    ┌──────┴──────┐           │
     │                    ▼             ▼           │
     │             ┌───────────┐ ┌───────────┐     │
     │             │Symptomatic│ │Asymptomatic│     │
     │             └───────────┘ └───────────┘     │
     │                                              │
     └──────────────────────────────────────────────┘

     Risk factors modify CST transition probabilities:
       - Douching (RR 1.21)
       - Poor menstrual hygiene (RR 4.1)
       - Concurrent partners (RR 1.28)
       - Uncircumcised partner (RR 4.0)
```

CST transitions are spontaneous and probabilistic. CST 4 (polymicrobial/BV) is the "infected" state. The stable CST represents each woman's baseline microbiome, which she tends to return to after clearance.

## Parameters

Pars of `sti.BV`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Microbiome states

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `stable_cst_distribution` | CST 1: 0.10, CST 3: 0.30, CST 4: 0.60 | Distribution of baseline (stable) CST | [Munoz 2021](https://doi.org/10.1186/s40168-021-01096-9) |
| `p_cst_change_cst1` | 0.10 per month | Probability of a change from CST 1 to CST 3 | Not recorded |
| `p_cst_change_cst3` | 0.05 per month | Probability of a change from CST 3 to CST 4 (BV) | Not recorded |
| `rr_stable_cst1` | 0.1 | Relative risk of CST transition for women with stable CST 1 | Not recorded |
| `rr_stable_cst3` | 1.0 | Relative risk of CST transition for women with stable CST 3 (reference) | Not applicable |
| `rr_stable_cst4` | 2.0 | Relative risk of CST transition for women with stable CST 4 | Not recorded |

### Clearance

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur2clear_cst3` | uniform(1 wk, 18 wk) | Duration of a transient CST 3, in women with stable CST 1 | Not recorded |
| `dur2clear_cst4` | uniform(5 wk, 50 wk) | Duration of untreated BV, after which women return to CST 1 (stable CST 1) or CST 3 | [Klebanoff 2004](https://doi.org/10.1016/j.ajog.2003.08.020); [Tamarelle 2022](https://doi.org/10.1097/OLQ.0000000000001662) |
| `dur_post_tx` | constant(0) | Time after treatment before return to the stable CST (0: at once) | Assumption |
| `recur_rate` | 0.06 per month | Rate of recurrence just after treatment | Assumption |
| `dur_recur_decay` | 3 mo | Time over which the rate of recurrence falls by a factor of e | Assumption |

### Symptoms

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_symp` (stable_cst1) | 0.80 | Probability of symptoms, stable CST 1 | Not recorded |
| `p_symp` (stable_cst3) | 0.70 | Probability of symptoms, stable CST 3 | Not recorded |
| `p_symp` (stable_cst4) | 0.10 | Probability of symptoms, stable CST 4 | Not recorded |
| `dur_presymp` | uniform(1 wk, 2 wk) | Time to symptom onset | Not recorded |

### Behavioral risk factors

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `p_douching` | 0.64 | Proportion of women who douche | [Bassey 2021](https://doi.org/10.1016/j.jpag.2020.06.003); [Ziba 2019](https://doi.org/10.1186/s12905-019-0720-7) |
| `rr_douching` | 1.21 | Relative risk of BV from douching | [Brotman 2008](https://doi.org/10.1093/aje/kwn103) |
| `p_poor_menstrual_hygiene` | 0.55 | Proportion of women with poor menstrual hygiene | [Anbesu 2023](https://doi.org/10.1186/s12889-022-14942-8) |
| `rr_poor_menstrual_hygiene` | 4.1 | Relative risk of BV from poor menstrual hygiene | Calibrated (Kenya) |
| `n_partners` | poisson(0.05) | Number of current partners, used when `count_partners` is False | Not recorded |
| `count_partners` | False | Whether concurrency is counted from the sexual network rather than drawn from `n_partners` | Not applicable |
| `rr_concurrency` | 1.28 | Relative risk of BV from concurrent partners | [Kenyon 2016](https://doi.org/10.1136/sextrans-2016-052652) |
| `p_circumcised` | 0.40 | Proportion of male partners circumcised | Not recorded |
| `rr_uncircumcised` | 4.0 | Relative risk of BV from an uncircumcised partner | Calibrated (Kenya) |

### Transmission

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `beta_m2f` | 0 | Not used: BV is not transmitted | Not applicable |
| `rel_beta_f2m` | 0.5 | Not used: BV is not transmitted | Not applicable |
| `beta_m2c` | 0 | Not used: BV is not transmitted | Not applicable |
| `beta_m2m` | 0 | Not used: BV is not transmitted | Not applicable |
| `eff_condom` | 0.45 | Condom effectiveness in reducing BV risk | [Hutchinson 2007](https://doi.org/10.1097/EDE.0b013e3181567eaa) |
