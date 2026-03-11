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
    test_prob_data=0.1,         # Annual testing rate (10% per year)
    eligibility=lambda sim: sim.people.female,  # Who is eligible
    start=2015,                 # When testing begins
)
```

Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `product` | (required for STITest; HIVTest auto-creates a perfect HIVDx) | Diagnostic product (e.g., `STIDx`, `HIVDx`) |
| `test_prob_data` | 1.0 | Annual testing probability (scalar or array over years). Converted to per-timestep probability via `ss.probperyear` when `dt_scale=True`. |
| `eligibility` | all agents (HIVTest: undiagnosed only) | Function `f(sim) -> BoolArr` defining who can be tested, e.g. `lambda sim: sim.people.female & (sim.people.age >= 15)` |
| `start` / `stop` | first/last sim year | Calendar year (inclusive) when the intervention activates/deactivates |
| `dt_scale` | True | If True (default), interpret `test_prob_data` as an annual probability. Set to False to interpret as a per-timestep probability. |

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

## HIV treatment pipeline

HIV treatment in STIsim follows a specific pipeline. Understanding this flow is essential for setting up HIV interventions correctly:

```
HIVTest                          ART
  │                               │
  │ test eligible agents          │ check newly diagnosed
  │ (per-year probability)        │ (ti_diagnosed == this step)
  │                               │
  ▼                               ▼
hiv.diagnosed = True    ──►   art_initiation filter (90% default)
hiv.ti_diagnosed = ti         │
                              ▼
                         If coverage specified:
                           prioritize by CD4, match target
                         If no coverage:
                           treat all who initiate
                              │
                              ▼
                         hiv.on_art = True
                         (reduces transmissibility, reconstitutes CD4)
```

**Key points:**

- **Intervention order matters**: interventions run in the order they appear in the list. HIVTest must come before ART so that agents diagnosed this timestep can initiate ART in the same step.
- Agents must be **diagnosed** before they can go on ART. If you add ART without HIVTest, a warning is raised and no agents will be treated.
- `test_prob_data` is an **annual testing probability** by default (`dt_scale=True`). It is converted to a per-timestep probability via `ss.probperyear`. With monthly timesteps, `test_prob_data=0.1` means ~0.88% per month, correctly recovering ~10% per year. To use a per-timestep probability directly, set `dt_scale=False`.
- ART `art_initiation` (default 0.9) controls what fraction of newly diagnosed agents are willing to start treatment.
- If `coverage` is provided, ART force-fits the number on treatment to match the target each timestep. Agents are added or removed immediately (not gradually) to hit the target, prioritized by CD4 count and care-seeking propensity.
- If no `coverage` is provided (the default), ART simply treats everyone who initiates, with no capacity constraint.

### HIVTest

```python
# Test 20% of undiagnosed agents per year
hiv_test = sti.HIVTest(test_prob_data=0.2, start=2000, name='hiv_test')

# Higher rate for FSWs
fsw_test = sti.HIVTest(
    test_prob_data=0.5,
    name='fsw_test',
    eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed,
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `test_prob_data` | 1.0 | Annual testing probability (converted via `ss.probperyear` when `dt_scale=True`) |
| `eligibility` | undiagnosed | Function `f(sim) -> BoolArr` defining who can be tested |
| `start` | sim start | Year testing begins |

### ART

ART accepts coverage in multiple formats:

```python
# No coverage target — 90% of newly diagnosed initiate ART (default)
art = sti.ART()

# Treat ALL diagnosed, no coverage constraint
art = sti.ART(art_initiation=1)

# Constant proportion of infected
art = sti.ART(coverage=0.8)

# Time-varying proportion
art = sti.ART(coverage={'year': [2000, 2010, 2025], 'value': [0, 0.5, 0.9]})

# From a DataFrame (absolute numbers)
art = sti.ART(coverage=pd.read_csv('art_data.csv').set_index('year'))

# Age/sex stratified (columns: Year, Gender/Sex, AgeBin, + value)
art = sti.ART(coverage=stratified_df)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `coverage` | None (no target; treat all who initiate) | Coverage target: scalar, dict, DataFrame, or stratified DataFrame. Dict values are linearly interpolated between years. DataFrame requires index=years and a column named `n_art` (absolute numbers) or `p_art` (proportion of infected). Stratified DataFrames need columns Year, Gender/Sex, AgeBin `[lo,hi)`, and a numeric value column. Gender accepts: 0/f/female (female), 1/m/male (male). |
| `art_initiation` | bernoulli(0.9) | Probability a newly diagnosed person initiates ART. Set to 1 to treat all diagnosed. |

### VMMC

Voluntary medical male circumcision. Targets all males by default; reduces susceptibility to HIV acquisition by 60%.

```python
vmmc = sti.VMMC(coverage=0.3)
vmmc = sti.VMMC(coverage={'year': [2010, 2025], 'value': [0, 0.4]})
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `coverage` | None (does nothing without coverage data) | Coverage target (same formats as ART). Required for VMMC to have any effect. |
| `eff_circ` | 0.6 | Efficacy (60% reduction in HIV acquisition risk) |
| `eligibility` | all males | Optional function to further restrict who is eligible |

### PrEP

Pre-exposure prophylaxis. By default targets HIV-negative FSWs (female sex workers) who are not already on PrEP. Use the `eligibility` parameter to target a different population.

```python
prep = sti.Prep(coverage=[0, 0.5], years=[2020, 2025])
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `coverage` | [0, 0.01, 0.5, 0.8] | Coverage values at each year (linearly interpolated between years) |
| `years` | [2004, 2005, 2015, 2025] | Calendar years corresponding to coverage values |
| `eff_prep` | 0.8 | Efficacy (80% reduction in HIV acquisition risk) |
| `eligibility` | HIV-negative FSWs not on PrEP | Optional function to override default targeting |

## Combining interventions

Interventions are passed to the sim as a list:

```python
sim = sti.Sim(
    diseases='ng',
    interventions=[test, treatment],
)
```

For a complete worked example, see the [Interventions tutorial](../../tutorials/tut_interventions.ipynb).
