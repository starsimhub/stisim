# STIsim HIV Module – GitHub Issues

> **Convention:**
> - **Milestones** = Phases below
> - **Features** = GitHub Issues (`label: HIV`, assigned to milestone)
> - **Tasks** = Sub-issues nested under each Feature issue

---

## Milestone: Phase 1 – Structure, Defaults & Templates

**Success criteria:**
- New users can create and run a sim, plot/export outputs, and change parameters
- Python and R users are supported
- AI support integration is available
- New devs understand how and where to contribute (test format, contribution guide, CoC)

---

### Feature: Repo Structure
`label: HIV` `milestone: Phase 1`

Set up and finalize the repository structure for the HIV module.

- [ ] Decide on and implement repo structure
- [ ] Make project board / roadmap public
- [ ] Add links from README to project board
- [ ] Consider adding a `ROADMAP.md` in root folder

---

### Feature: API for Creating an HIV Sim
`label: HIV` `milestone: Phase 1`

Implement the core API so users can instantiate and run an HIV sim with a single line.

- [ ] Ensure `sim = hiv.Sim(location='zambia')` works end-to-end

---

### Feature: Data Loading & Sim Creation
`label: HIV` `milestone: Phase 1`

Create methods for loading, caching, and overwriting simulation data, with defaults for countries with and without data.

- [ ] Create methods for data loading, caching, and overwriting (assumes data exists)
- [ ] Implement default parameter sets for ~1–2 countries with data
- [ ] Implement generic parameter sets for countries without data
- [ ] Implement calibration target formatting system

---

### Feature: Default Sim Plotting
`label: HIV` `milestone: Phase 1`

Ensure the sim produces useful default plots out of the box.

- [ ] Implement default sim plotting

---

### Feature: User Onboarding Materials
`label: HIV` `milestone: Phase 1`

Documentation and resources to help new users get started quickly.

- [ ] Write detailed documentation of model defaults and assumptions
- [ ] Write / include AI support integration
- [ ] Update README with install instructions for both Python and R
- [ ] Point to 1–2 example repos

---

### Feature: Dev Onboarding Materials
`label: HIV` `milestone: Phase 1`

Resources to help new contributors understand how and where to contribute.

- [ ] Document preferred format for tests
- [ ] Document anything incremental over the Starsim user guide
- [ ] Write Code of Conduct / contribution guide

---

## Milestone: Phase 2-RES – Scientific Validation

**Success criteria:**
- Model landscaping complete
- Parameters validated
- Structure validated
- Model comparison ready

> *(No feature-level tasks defined yet — to be scoped in planning.)*

---

## Milestone: Phase 2-SW – Test-Driven Development & Validation

---

### Feature: Test Epi / Model Dynamics
`label: HIV` `milestone: Phase 2-SW`

Validate natural history and disease progression without treatment.

- [ ] Test CD4 counts decline over time without treatment
- [ ] Test median time from infection to AIDS
- [ ] Test untreated HIV mortality follows expected patterns
- [ ] Test perinatally infected children progress faster

---

### Feature: Test Transmission Dynamics
`label: HIV` `milestone: Phase 2-SW`

Validate transmission probabilities, sensitivity, and network effects.

- [ ] Test transmission probability increases with viral load
- [ ] Test acute phase has elevated transmission
- [ ] Test late-stage HIV has higher transmission than chronic phase
- [ ] Test MTCT rate without any intervention
- [ ] Test MTCT risk increases with maternal VL
- [ ] Test postnatal transmission through breastfeeding
- [ ] Test KP prevalence is higher than general population
- [ ] Test sensitivity to transmission probability (beta)
- [ ] Test sensitivity to initial prevalence
- [ ] Test sensitivity to network structure parameters
- [ ] Test sensitivity to CD4 decline rate
- [ ] Test sensitivity to ART efficacy assumptions
- [ ] Test that doubling beta approximately doubles incidence (at low prevalence)
- [ ] Test that removing sexual networks eliminates transmission

---

### Feature: Test Treatment Effects
`label: HIV` `milestone: Phase 2-SW`

Validate ART effects and treatment interruption dynamics.

- [ ] Test U=U (undetectable = untransmittable)
- [ ] Test that ART suppresses viral load
- [ ] Test that ART reduces HIV mortality
- [ ] Test that CD4 counts increase on ART
- [ ] Test that poor adherence leads to virologic failure
- [ ] Test viral rebound after treatment interruption
- [ ] Test MTCT rate with maternal ART (Option B+)

---

### Feature: Test Network & Sexual Behavior
`label: HIV` `milestone: Phase 2-SW`

Validate partnership formation, concurrency, and degree distribution.

- [ ] Test that partnership formation matches demographic expectations
- [ ] Test age-disparate relationships
- [ ] Test concurrent partnership prevalence
- [ ] Test partnership duration by type
- [ ] Test sexual partner degree distribution

---

### Feature: Test Other Interventions
`label: HIV` `milestone: Phase 2-SW`

Validate prevention interventions and counterfactual scenarios.

- [ ] Test voluntary medical male circumcision effectiveness
- [ ] Test PrEP effectiveness with adherence
- [ ] Test that increased testing leads to earlier diagnosis
- [ ] Test "no intervention" counterfactual produces higher mortality
- [ ] Test that removing ART returns to pre-ART mortality levels
- [ ] Test that ART interruption scenario shows mortality rebound
- [ ] Test "early ART" scenario (e.g., 2000 vs 2004) reduces cumulative deaths
- [ ] Test that achieving 95-95-95 by 2025 reduces incidence to threshold levels
- [ ] Test pessimistic vs optimistic scenarios produce expected outcome ranges

---

### Feature: Test Epidemic Patterns
`label: HIV` `milestone: Phase 2-SW`

Validate macro-level epidemic dynamics and demographic stratifications.

- [ ] Test exponential growth in early epidemic
- [ ] Test age-specific prevalence patterns
- [ ] Test HIV prevalence by sex

---

### Feature: Validate Outputs
`label: HIV` `milestone: Phase 2-SW`

Validate that all standard outputs are calculated correctly and accessible.

- [ ] Test all standard HIV metrics (prevalence, incidence, PLHIV, deaths, ART coverage)
- [ ] Test age-disaggregated outputs are available and stored correctly
- [ ] Test sex-disaggregated outputs are available
- [ ] Test outputs can be stratified by key populations
- [ ] Test cumulative metrics (e.g., total infections averted)
- [ ] Test results can be exported in multiple formats (CSV, JSON)
- [ ] Test that results include timestamped metadata

---

### Feature: Intervention Scale-Up / Coverage
`label: HIV` `milestone: Phase 2-SW`

Validate that intervention scale-up dynamics behave correctly.

- [ ] Test that intervention effects match meta-analysis estimates from clinical trials

---

### Feature: Test Edge Cases & Robustness
`label: HIV` `milestone: Phase 2-SW`

Ensure the model handles degenerate and extreme inputs gracefully.

- [ ] Test zero infections scenario
- [ ] Test 100% prevalence scenario
- [ ] Test very small and very large populations
- [ ] Test unusual/extreme parameter values
- [ ] Test missing data handling

---

## Milestone: Phase 3-I – Model Comparison: Eswatini

**Success criteria:**
- Serves as a fit-for-purpose replacement for Adam's Eswatini HIV analyses (currently using EMOD)

> *(Feature-level tasks to be scoped in planning.)*

---

## Milestone: Phase 3-II – Model Comparison: Kenya

**Success criteria:**
- Replicates or is close to Monisha's EMOD model

---

### Feature: Kenya Model Comparison
`label: HIV` `milestone: Phase 3-II`

- [ ] Test that key outputs are comparable to EMOD
- [ ] Test that MTCT rates match WHO/UNAIDS estimates by intervention coverage

---

## Milestone: Phase 3-III – Model Comparison: Zambia

**Success criteria:**
- Replicates or is close to the NYU EMOD Zambia model

---

### Feature: Zambia Model Comparison
`label: HIV` `milestone: Phase 3-III`

- [ ] Test that key outputs are comparable to EMOD
- [ ] Test that MTCT rates match WHO/UNAIDS estimates by intervention coverage
