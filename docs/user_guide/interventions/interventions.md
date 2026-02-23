# Interventions

STIsim separates the **product** (what is being delivered) from the **intervention** (how it is delivered). A diagnostic product stores the characteristics of a test -- its sensitivity, specificity, and what states it can detect. An intervention stores the delivery information -- who is eligible, how often they are reached, and when the program starts and stops. This separation means you can swap out a diagnostic product without changing the delivery strategy, or target the same product to different populations.

## Architecture

```
Intervention (delivery)          Intervention (delivery)
  - eligibility                    - eligibility
  - test probability               - treatment efficacy
  - start/stop years               - start/stop years
        │                                │
        ▼                                ▼
  Diagnostic product               Disease module
  (product characteristics)        (clears infection)
  - sensitivity by state
  - specificity
```

Testing and treatment are separate interventions. A testing intervention uses a diagnostic product to determine outcomes, and can optionally trigger a treatment intervention for positive results. Alternatively, treatment can be applied directly (e.g., syndromic management without a lab test).

## Testing

### STITest

The base testing class. Controls who gets tested, how often, and with what diagnostic.

```python
test = sti.STITest(
    product=my_diagnostic,      # Diagnostic product (required)
    test_prob_data=0.1,         # Per-timestep probability of testing
    eligibility=lambda sim: sim.people.female,  # Who is eligible
    start=2015,                 # When testing begins
)
```

Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `product` | (required) | Diagnostic product (e.g., `STIDx`, `HIVDx`) |
| `test_prob_data` | 1.0 | Testing probability -- scalar, array, or time-varying |
| `eligibility` | all agents | Function `f(sim) -> UIDs` defining who can be tested |
| `start` / `stop` | sim start/end | Active period for the intervention |
| `dt_scale` | True | Scale probability by timestep length |

### HIVTest

HIV-specific testing. Defaults to testing only undiagnosed agents and using a perfect diagnostic.

```python
hiv_test = sti.HIVTest(test_prob_data=0.1, start=2000)
```

On a positive result, sets `sim.diseases.hiv.diagnosed = True`.

### SymptomaticTesting

Syndromic management: test symptomatic care-seekers and route them to treatment based on results. Used for bacterial STIs where lab confirmation may not be available.

```python
syndromic = sti.SymptomaticTesting(
    diseases=[ng, ct, tv],
    treatments=[ng_tx, ct_tx, metro],
    disease_treatment_map={'ng': ng_tx, 'ct': ct_tx, 'tv': metro},
)
```

## Treatment

### STITreatment

The base treatment class. Clears infection for eligible agents with a given efficacy.

```python
ct_tx = sti.STITreatment(
    diseases='ct',         # Which disease(s) to treat
    treat_eff=0.95,        # 95% efficacy
    name='ct_tx',
)
```

Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `diseases` | (required) | Disease name or list of names to treat |
| `treat_prob` | 1.0 | Probability of accepting treatment |
| `treat_eff` | 0.9 | Treatment efficacy |
| `eligibility` | all agents | Function defining who can be treated |

Treatment outcomes are tracked as successful (infection cleared), failed (infection persists), or unnecessary (treated a susceptible agent).

### Disease-specific treatments

| Class | Disease | Notes |
|-------|---------|-------|
| `sti.STITreatment` | Any STI | Generic; specify disease(s) by name |
| `sti.GonorrheaTreatment` | Gonorrhea | Tracks antimicrobial resistance via `rel_treat` |
| `sti.SyphTx` | Syphilis | Treats pregnant mothers and fetuses; stage-specific |
| `sti.treat_BV` | BV | CST-based treatment with durable post-treatment effects |

## HIV-specific interventions

### ART

Antiretroviral therapy with CD4-based prioritization and coverage targets.

```python
art = sti.ART(coverage_data=art_df)  # DataFrame with year index and n_art or p_art column
```

ART reduces transmissibility by 96% (default) and reconstitutes CD4 counts. Coverage can be specified as absolute numbers or proportions.

### VMMC

Voluntary medical male circumcision. Reduces male susceptibility to HIV by 60% (default).

```python
vmmc = sti.VMMC(coverage_data=vmmc_df)
```

### PrEP

Pre-exposure prophylaxis targeting FSWs. Reduces susceptibility by 80% (default).

```python
prep = sti.Prep(coverage=[0, 0.5], years=[2020, 2025])
```

## Combining interventions

Interventions are passed to the sim as a list:

```python
sim = sti.Sim(
    diseases='ng',
    interventions=[test, treatment],
)
```

For a complete worked example, see the [Interventions tutorial](../../tutorials/tut_interventions.ipynb).
