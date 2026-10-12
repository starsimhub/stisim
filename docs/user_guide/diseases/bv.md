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
| `p_cst_change` | cst1: 0.1, cst3: 0.05 | Probability of transition to a worse CST (CST 1 to 3; CST 3 to 4) | Not recorded |
| `rr_stable_cst1` | 0.1 | Relative risk of CST transition for women with stable CST 1 | Not recorded |
| `rr_stable_cst3` | 1.0 | Relative risk of CST transition for women with stable CST 3 (reference) | Not applicable |
| `rr_stable_cst4` | 2.0 | Relative risk of CST transition for women with stable CST 4 | Not recorded |

### Clearance

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur2clear` (cst3) | uniform(1 wk, 18 wk) | Time to return to the stable CST from CST 3 | Not recorded |
| `dur2clear` (cst4) | uniform(5 wk, 50 wk) | Time to return to the stable CST from CST 4 | [Klebanoff 2004](https://doi.org/10.1016/j.ajog.2003.08.020) |
| `dur2clear` (male) | uniform(1 wk, 18 wk) | Duration of male colonization | Not recorded |
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
| `sexual_transmission` | False | Whether BV-associated bacteria are sexually transmitted | Not applicable |
| `beta_m2f` | None | Per-act male-to-female transmission probability (used if `sexual_transmission` is True) | Not applicable |
| `rel_beta_f2m` | 0.5 | Female-to-male per-act probability relative to `beta_m2f` | Not applicable |
| `beta_m2c` | None | Mother-to-child transmission probability per timestep | Not applicable |
| `beta_breastfeed` | None | Breastfeeding transmission probability per timestep | Not applicable |
| `beta_m2m` | None | Per-act male-to-male transmission probability | Not applicable |
| `eff_condom` | 0.45 | Condom effectiveness in reducing BV risk | [Hutchinson 2007](https://doi.org/10.1097/EDE.0b013e3181567eaa) |

### Initialization

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `init_prev` | 0.23 | Initial BV prevalence | [WHO BV fact sheet](https://www.who.int/news-room/fact-sheets/detail/bacterial-vaginosis) |
| `rel_init_prev` | 1 | Multiplier on initial prevalence | Not applicable |

### Pregnancy outcomes

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `or_ptb` (trimester 1) | 3 | Odds ratio for preterm birth, BV in 1st trimester | [Juliana 2020](https://doi.org/10.3389/fpubh.2020.567885); [Dunlop 2021](https://doi.org/10.3389/fcimb.2021.641005) |
| `or_ptb` (trimester 2) | 2 | Odds ratio for preterm birth, BV in 2nd trimester | [Juliana 2020](https://doi.org/10.3389/fpubh.2020.567885) |
| `or_ptb` (trimester 3) | 1 | Odds ratio for preterm birth, BV in 3rd trimester | [Juliana 2020](https://doi.org/10.3389/fpubh.2020.567885) |
