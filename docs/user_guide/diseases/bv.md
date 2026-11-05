# Bacterial vaginosis (BV)

BV is modeled using community state types (CSTs) representing vaginal microbiome composition. Unlike the other STIs, BV can arise spontaneously (without sexual transmission) and is driven by behavioral risk factors.

**Class:** `sti.BV` | **Alias:** `'bv'` | **Base class:** `ss.Disease` (declared: its CSTs and transitions, with the risk factors as code)

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

CST transitions are spontaneous, with a probability per month raised by each woman's risk factors. CST 4 (polymicrobial/BV) is the "infected" state, as the declared states `presymptomatic`, `asymptomatic`, and `symptomatic` (`cst4` is their union, for women over 15); `cst1` and `cst3` are states too. The stable CST represents each woman's baseline microbiome, drawn at birth, which she returns to after clearance (a woman whose stable CST is 4 returns to CST 3, and starts there as a girl). The risk factors, and the probabilities of a worse CST they give, are code (as agents).

## Parameters

Pars of `sti.BV`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Microbiome states

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `stable_cst_distribution` | CST 1: 0.10, CST 3: 0.30, CST 4: 0.60 | Distribution of baseline (stable) CST | [Munoz 2021](https://doi.org/10.1186/s40168-021-01096-9) |
| `p_cst_change` | cst1: 0.1, cst3: 0.05 | Probability of transition to a worse CST (CST 1 to 3; CST 3 to 4) | Not recorded |
| `rr_stable_cst1` | 0.1 | Relative risk of CST transition for women with stable CST 1 | Not recorded |
| `rr_stable_cst3` | 1.0 | Relative risk of CST transition for women with stable CST 3 (reference) | Not applicable |
| `rr_stable_cst4` | 2.0 | Relative risk of CST transition for women with stable CST 4 | Not recorded |

### Clearance

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur_cst3` | uniform(1 wk, 18 wk) | Time in CST 3 before a woman whose stable CST is 1 returns to it | Not recorded |
| `dur_cst4` | uniform(5 wk, 50 wk) | Time in CST 4 before it clears (from symptom onset, for those with symptoms) | [Klebanoff 2004](https://doi.org/10.1016/j.ajog.2003.08.020) |
| `dur_male` | uniform(1 wk, 18 wk) | Duration of male colonization | Not recorded |
| `spontaneous_clearance` | cst1→1, cst3→3, cst4→3 | CST reached after spontaneous clearance, by stable CST | [Tamarelle 2022](https://doi.org/10.1097/OLQ.0000000000001662) |

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
| `beta_m2f` | 0 | Per-act male-to-female transmission probability of BV-associated bacteria (none by default) | Not applicable |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not applicable |
| `eff_condom` | 0.45 | Condom effectiveness in reducing BV risk | [Hutchinson 2007](https://doi.org/10.1097/EDE.0b013e3181567eaa) |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0 | Not used: each woman's CST is set by `assign` (from `stable_cst_distribution`), at the start and at birth | Not applicable |

### Pregnancy outcomes

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `or_ptb` (trimester 1) | 3 | Odds ratio for preterm birth, BV in 1st trimester | [Juliana 2020](https://doi.org/10.3389/fpubh.2020.567885); [Dunlop 2021](https://doi.org/10.3389/fcimb.2021.641005) |
| `or_ptb` (trimester 2) | 2 | Odds ratio for preterm birth, BV in 2nd trimester | [Juliana 2020](https://doi.org/10.3389/fpubh.2020.567885) |
| `or_ptb` (trimester 3) | 1 | Odds ratio for preterm birth, BV in 3rd trimester | [Juliana 2020](https://doi.org/10.3389/fpubh.2020.567885) |
