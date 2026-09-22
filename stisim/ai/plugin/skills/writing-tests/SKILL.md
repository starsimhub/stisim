---
name: writing-tests
description: Use when writing or reviewing tests for HIVsim / STIsim development or downstream analysis code. Optimises for a small number of scientifically meaningful tests over comprehensive enumeration of trivial ones. Every proposed test should have a one-sentence answer to "what meaningful bug would this catch?"; if the answer is "it confirms a value we assigned still has that value", do not add it.
---

# Writing tests for HIVsim / STIsim code

*This skill is intentionally compact. Refine as we gain more experience with which testing patterns hold up on agent-generated scientific-model code.*

## When to use

- Adding or reviewing tests for HIVsim / STIsim itself (upstream package tests).
- Adding tests to a downstream analysis repo — custom interventions, analyzers, disease modules, eligibility / targeting logic, data preprocessing whose errors would change scientific results.
- Diagnosing a bug — writing the regression test that captures it.
- Trigger phrases: "add a test", "write tests for X", "cover this with a test", "why doesn't this have a test", "make sure this is tested".

## When NOT to use

- Reviewing code where no test change is being proposed — that's just code review.
- The user is running an existing test suite, not writing new ones.

## The one rule

**Every new test should have a one-sentence answer to: "what meaningful bug would this test catch?"**

- If the answer is *"it confirms that a value we assigned still has the value we assigned"* — do not add the test.
- If the answer is *"it catches a vaccination intervention that selects the right agents but fails to actually reduce their susceptibility"* — that is a useful test.
- If the answer is *"it catches a disease intervention whose internal state changes but produces no effect on transmission because the disease module never reads that state"* — that is an especially valuable test.

Optimise for a small set of high-value tests. Coverage metrics may be informative but must not drive the creation of meaningless simulation tests.

## Framing

### Upstream package tests vs downstream analysis tests

Two places tests may live, with different rules:

- **Upstream (HIVsim / STIsim itself)**: TDD-style — reproduce the failure in a minimal test, watch it fail, implement the fix, verify it passes, run the surrounding suite for regressions. The test captures a **general invariant of the package**, not the details of the downstream analysis that revealed it. Regression tests get names describing the invariant (`test_art_target_can_change_from_count_to_proportion`), not the story (`test_exp_5_art_bug`). This is the `comment-hygiene` rule applied to test names.
- **Downstream (analysis repository)**: more freedom to reflect the specific analysis, but still understandable and scientifically meaningful. Strongly encourage tests for custom interventions, analyzers, disease modules, eligibility logic, targeting logic, scenario utilities, and data-processing functions whose errors would affect scientific results.

If work on a downstream analysis uncovers a generic upstream bug, the fix belongs upstream *with* a general-invariant test — see `extending-stisim`.

### Scientific validity is the priority

Do not stop at "the intervention initialises without error." For a new mechanism, at minimum consider tests at two levels:

1. **Mechanism** — does the intervention actually change the intended per-agent state? (Reduced `rel_sus`, updated infection state, correct eligibility mask, correct targeting.)
2. **Population outcome** — in a *controlled scenario* designed to make the effect large and detectable, does the intervention alter the relevant epidemiological outcome (incidence / prevalence / infections / person-time infected) versus a paired no-intervention control?

Both levels catch different failure modes. A mechanism test alone won't catch a state that the disease module never reads; a population-outcome test alone won't isolate whether the failure was mechanism, targeting, or transmission wiring.

### Do not encode simplistic monotonic expectations

Simulation outcomes can be counterintuitive. Do not automatically assert *"more diagnosis must reduce incidence"* — for some STIs, faster diagnosis and treatment shortens infection duration and increases the susceptible pool, and under some circumstances can raise incidence via turnover.

Before writing an outcome-direction test, ask:

- Is this direction theoretically guaranteed by the mechanism?
- Or is it merely what we usually expect?
- Could a competing mechanism reverse it?
- Is it robust to stochastic variation at the population size the test uses?

If the sign is not guaranteed, test the direct mechanism instead, or design a simplified controlled scenario where the direction is unambiguous.

## Instructions

1. **Identify where the test belongs.** Upstream package or downstream analysis repo? If a downstream discovery reveals a generic bug, the test belongs upstream and captures the general invariant, not the analysis story.

2. **Inspect existing tests and reuse fixtures.** For upstream work look at `stisim/tests/testlib.py` (the `build_testing_sim` helper), `stisim/tests/test_hiv.py`, and neighbouring `test_*.py` files. Match naming, organisation, and assertion conventions. Prefer minimal existing demo (`ss.demo(...)`, `sti.Sim(...)`, `build_testing_sim(...)`) + small modification over hand-building a full simulation stack for one test. Do not introduce a new test framework or elaborate fixture architecture unless there is a clear need.

3. **Design the smallest controlled scenario that exercises the behaviour.** Rather than running a realistic country model with dozens of competing processes, use small populations, short durations, simplified networks, minimal modules, and make the intervention effect large enough to be detectable. `stisim/tests/test_hiv.py` conventions: `tiny_pop=10`, `small_pop=100`, `medium_pop=1000`, `large_pop=4000` — sized so assertions fail from stochastic noise <5% of the time. Reuse those sizes where possible.

4. **Prefer invariants over exact trajectories.** Avoid fragile assertions on exact incidence values or entire time series unless the process is intentionally deterministic. Robust assertions include: at least one eligible person was vaccinated; vaccinated people received the expected protection state; no ineligible people were vaccinated; protection values lie within valid bounds; the intervention arm differs materially from the control in a controlled setup; a quantity that must be non-negative is; mutually exclusive states never coexist; counts reconcile when required. Fix the seed if the package convention supports it.

5. **Handle stochasticity explicitly, not by inflating tolerances.** A test that fails randomly is harmful. Options: control the seed; simplify the scenario; increase the intervention effect; test the direct mechanism rather than the noisy aggregate endpoint; use a robust inequality (arm A > arm B by at least X) rather than an exact equality; use replicates only when the mechanism genuinely needs them. Do not "solve" flake by making tolerances arbitrarily enormous.

6. **Test public behaviour, not implementation plumbing.** If an intervention promises *vaccinate eligible individuals and reduce susceptibility*, test those behaviours. Do not couple to private helper names, temporary internal arrays, or exact implementation order unless those are themselves part of an intentional contract.

7. **Keep upstream tests fast.** Any test added to HIVsim / STIsim must run quickly — the whole suite should not become a 30-minute analysis. If a test needs a large population or many replicates just to observe the expected effect, the test design usually can be simplified — isolate the mechanism more directly.

8. **Custom analyzers, disease modules, and analysis-level smoke tests get specific attention:**
    - **Analyzers**: test that the recorded quantity, denominator / population, and time indexing agree with directly inspectable simulation state in a small controlled case. Use the framework's analyzer architecture; do not mutate core sim state to collect results (starsim explicitly treats analyzers as observers).
    - **Custom disease modules**: test initialisation, susceptibility / infection state transitions, transmission behaviour, recovery / treatment, immunity or altered susceptibility, result counting, and interaction with relevant interventions. `sim.run()` completes is *not* sufficient for a custom disease.
    - **Analysis repositories**: at least one inexpensive smoke test that constructs the core sim, runs a short simulation, and confirms essential outputs exist and are finite. Do not run the full calibrated national analysis in routine CI — the smoke test is there to catch broken configuration and API changes cheaply.

9. **Do not proliferate trivial tests.** Do not generate a separate test for every setter, constructor argument, or obvious line of code. Return to the one rule: what meaningful bug would this catch?

## Suggested workflow

1. Identify upstream vs downstream.
2. Inspect existing tests and fixture helpers.
3. Identify the specific behaviour or scientific invariant to protect.
4. Reuse existing demo / fixture infrastructure.
5. Construct the smallest controlled scenario that exercises the behaviour.
6. Write the test; where practical, verify it fails without the fix / behaviour.
7. Implement or verify the code.
8. Run the specific test, then the surrounding test set for regressions.
9. Check runtime; simplify if unnecessarily expensive.

## Checks before completion

1. Every new test has an answered "what meaningful bug would this catch?" question.
2. Upstream regression tests are named for the invariant, not the downstream story.
3. Existing demo / fixture infrastructure was reused rather than reinvented.
4. Controlled scenario is as small as possible while still exercising the mechanism.
5. Assertions are robust to stochastic variation without inflated tolerances.
6. Upstream tests are fast enough not to slow the suite meaningfully.
7. No trivial "value I assigned still equals value I assigned" tests were added.
