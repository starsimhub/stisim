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

### Microbiome states

| Parameter | Default | Description |
|-----------|---------|-------------|
| `stable_cst_distribution` | [1: 10%, 3: 30%, 4: 60%] | Distribution of baseline CST types |
| `p_cst_change` | cst1: 0.1, cst3: 0.05 | Probability of transitioning to worse CST |
| `rr_stable_cst1` | 0.1 | Relative risk of CST transition for CST 1 women |
| `rr_stable_cst3` | 1.0 | Relative risk for CST 3 (reference) |
| `rr_stable_cst4` | 2.0 | Relative risk for CST 4 women |

### Clearance

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dur2clear` (cst3) | uniform(1 wk, 18 wk) | Time to return to stable CST from CST 3 |
| `dur2clear` (cst4) | uniform(5 wk, 50 wk) | Time to return to stable CST from CST 4 |
| `spontaneous_clearance` | cst1→1, cst3→3, cst4→3 | CST to transition to after clearance |

### Symptoms

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_symp` (stable_cst1) | 0.80 | Probability of symptoms if stable CST 1 |
| `p_symp` (stable_cst3) | 0.70 | Probability of symptoms if stable CST 3 |
| `p_symp` (stable_cst4) | 0.10 | Probability of symptoms if stable CST 4 |
| `dur_presymp` | uniform(1 wk, 2 wk) | Time to symptom onset |
| `init_prev` | 0.23 | Initial BV prevalence |

### Behavioral risk factors

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_douching` | 0.64 | Proportion of women who douche |
| `rr_douching` | 1.21 | Relative risk of BV from douching |
| `p_poor_menstrual_hygiene` | 0.55 | Proportion with poor menstrual hygiene |
| `rr_poor_menstrual_hygiene` | 4.1 | Relative risk from poor menstrual hygiene |
| `rr_concurrency` | 1.28 | Relative risk from concurrent partners |
| `p_circumcised` | 0.40 | Proportion of men circumcised |
| `rr_uncircumcised` | 4.0 | Relative risk from uncircumcised partner |

### Pregnancy outcomes

| Parameter | Default | Description |
|-----------|---------|-------------|
| `or_ptb` (trimester 1) | 3 | Odds ratio for preterm birth, 1st trimester |
| `or_ptb` (trimester 2) | 2 | Odds ratio for preterm birth, 2nd trimester |
| `or_ptb` (trimester 3) | 1 | Odds ratio for preterm birth, 3rd trimester |
