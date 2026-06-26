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
    diseases=['ng', 'ct', 'tv'],
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

# Historical absolute numbers then projected proportions (mixed n/p format)
df = pd.read_csv('n_art.csv').set_index('year')
df['p_art'] = np.nan
df.loc[2023:, 'p_art'] = 0.90
art = sti.ART(coverage=df)

# Age/sex stratified (columns: Year, Gender/Sex, AgeBin, + value)
art = sti.ART(coverage=stratified_df)

# Smoother interpolation (default 0 = linear)
art = sti.ART(coverage={'year': [2000, 2025], 'value': [0, 0.9]}, smoothness=5)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `coverage` | None (no target; treat all who initiate) | Coverage target: scalar, dict, DataFrame, or stratified DataFrame. Dict values are interpolated between years (controlled by `smoothness`). DataFrame requires index=years and a column named `n_art` (absolute numbers) or `p_art` (proportion of infected). Dual-column DataFrames with both `n_art` and `p_art` support mixed format (use `format_priority` to resolve). Stratified DataFrames need columns Year, Sex, AgeBin `[lo,hi)`, and a numeric value column. Sex accepts: 0/f/female, 1/m/male. |
| `art_initiation` | bernoulli(0.9) | Probability a newly diagnosed person initiates ART. Set to 1 to treat all diagnosed. |
| `pmtct_efficacy` | 0.96 | Efficacy of maternal ART in reducing infant susceptibility to HIV, applied to both prenatal (MaternalNet) and postnatal (BreastfeedingNet) transmission. |
| `smoothness` | 0 | Interpolation smoothness for coverage data (0 = linear, higher = smoother S-curves). |
| `format_priority` | 'n' | When both `n_art` and `p_art` are non-NaN, prefer this format (`'n'` or `'p'`). |

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

## Syndromic management

`SyndromicManagement` treats vaginal/urethral discharge syndromes presumptively,
without a confirmatory lab test: symptomatic care-seekers are routed to treatment
based on the syndrome rather than a diagnosis. Use it for bacterial STIs in
settings where point-of-care diagnostics are unavailable.

```python
sm = sti.SyndromicManagement(
    diseases=['ng', 'ct', 'tv'],      # diseases by NAME (resolved against sim.diseases)
    treatments=[ng_tx, ct_tx, metro],
    outcome_tx_map={'discharge': [ng_tx, ct_tx]},
)
```

| Parameter | Description |
|-----------|-------------|
| `diseases` | Disease names presenting as the syndrome (e.g. `['ng','ct']`). Passed as names, not module objects. |
| `cervical_diseases` | Subset checked via the cervical-symptom path. |
| `treatments` / `outcome_tx_map` | Treatment interventions and the syndrome → treatment routing. |
| `treat_prob_data` | Probability a symptomatic care-seeker is managed. |

> **Stub** — expand with the symptom-detection logic and a worked discharge example.
> See [`SymptomaticTesting`](#symptomatictesting) for the test-and-treat variant, and
> the API reference for [`interventions.base_interventions`](../../api/interventions.base_interventions.qmd).

## Antenatal and infant screening (PMTCT)

STIsim supports single-visit antenatal care (ANC) screening that auto-schedules
newborn/infant follow-up, used for preventing mother-to-child transmission (PMTCT)
of HIV and congenital syphilis.

```python
anc = sti.ANCTest(
    disease_names=['hiv', 'syphilis'],   # auto-detects HIV + syphilis
    visit_prob=0.9,                      # ANC attendance probability
    disease_treatment_map={'hiv': art, 'syphilis': syph_tx},
    newborn_tests=[sti.InfantHIVTest(test_prob=0.8)],
)
```

| Class | Role |
|-------|------|
| `sti.ANCTest` | Multi-disease ANC visit testing, scheduled once per pregnancy; per-disease sensitivity; routes positives to treatment. |
| `sti.InfantHIVTest` | HIV test for infants born to mothers diagnosed in pregnancy (scheduled by `ANCTest` or `HIVTest` via the maternal network). |
| `sti.ANCSyphTest` | Syphilis-specific ANC screening with a diagnostic product; can schedule a newborn test. |
| `sti.NewbornSyphTest` | Syphilis test for newborns of mothers diagnosed in pregnancy. |
| `sti.NewbornTreatment` | Treatment for congenital syphilis in newborns. |

> **Stub** — expand with the ANC → treatment → newborn-test cascade, sensitivity
> defaults, and a PMTCT worked example. Requires a `Pregnancy` module (and
> `MaternalNet`/`BreastfeedingNet` for infant scheduling). See
> [`interventions.hiv_interventions`](../../api/interventions.hiv_interventions.qmd) and
> [`interventions.syphilis_interventions`](../../api/interventions.syphilis_interventions.qmd).

## Partner notification

`PartnerNotification` reaches the sexual partners of newly diagnosed index cases
and offers them follow-up testing. It works over two channels — the current sexual
network and an optional prior-partner recall network — each with separate
notification and attendance probabilities.

```python
pn = sti.PartnerNotification(
    eligibility=lambda sim: sim.diseases.hiv.ti_diagnosed == sim.ti,  # index cases
    test=sti.HIVTest(test_prob_data=1.0),                            # offered to attendees
    pars=dict(p_notify_current=0.6, p_attends_current=0.5),
)
```

| Parameter | Description |
|-----------|-------------|
| `eligibility` | `f(sim) -> uids` returning index cases (e.g. those diagnosed this step). |
| `test` | Testing intervention scheduled for partners who attend follow-up. |
| `p_notify_*` / `p_attends_*` | Per-channel (current/prior) notification × attendance probabilities; may be callables for edge-type stratification (see `sti.pn_rates`). |

> **Stub** — expand with the two-channel cascade and per-channel result tracking.
> See the gallery example [Partner notification](../../examples/partner_notification.qmd).

## Pregnancy-driven risk reduction

`PregnancyRiskReduction` lowers sexual-risk behaviour during pregnancy and restores
it afterwards — useful when behaviour change during pregnancy materially affects
transmission.

```python
prr = sti.PregnancyRiskReduction(pars=dict(
    fsw_redux=True,          # clear FSW status while pregnant
    high_risk_redux=True,    # drop high-risk-group membership
    concurrency_redux=True,  # zero concurrency
))
```

During pregnancy the module optionally clears FSW status, drops high-risk-group
membership to `default_risk_group`, and/or zeros concurrency; each agent's prior
state is restored once the pregnancy ends.

> **Stub** — expand with parameter semantics and the gallery example
> [Pregnancy risk modifier](../../examples/pregnancy_risk_modifier.qmd).

## Combining interventions

Interventions are passed to the sim as a list:

```python
sim = sti.Sim(
    diseases='ng',
    interventions=[test, treatment],
)
```

For a complete worked example, see the [Interventions tutorial](../../tutorials/tut_interventions.qmd).
