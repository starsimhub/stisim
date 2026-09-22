---
name: analysis-selector
description: Use very early in research intake — before any tool is assumed — to classify what class of analytical method a research question actually calls for. Prevents routing a statistical, ML, causal-inference, or descriptive question into HIVsim / STIsim (or any transmission model) simply because the researcher entered through that workflow. Sequence is research question → analytical objective → need for dynamics → need for mechanistic transmission → required representation of heterogeneity → data sufficiency → specific tool.
---

# Analysis selector

## When to use

- Very early in research intake — after the provisional research question is stated clearly enough to classify, but before any tool-specific configuration begins.
- Called from `analysis-intake` once the first round or two has settled the provisional question, populations, and decision purpose.
- The researcher entered through an HIVsim / STIsim workflow but the question does not obviously require dynamic transmission modelling.
- The researcher asks "is HIVsim right for this?" or "should I use HIVsim or something else?".
- Trigger phrases: "should I use HIVsim for this", "is this an HIVsim question", "what's the right method", "help me decide the approach", "scope the method before I start".

## When NOT to use

- The research question has already been triaged through this skill in an earlier session and the decision is recorded in the analysis brief.
- The user is inside a downstream skill (`calibration-workflow`, `model-writer`, etc.) and the method choice is already settled.
- The question is a follow-up refinement on an already-scoped analysis, not a new scoping question.

## Framing

**Start from the question, not from the tool.** A researcher arriving through an HIVsim workflow may nevertheless have a question best answered with regression, machine learning, clustering, descriptive statistics, causal inference, a compartmental model, or a time-series forecast. Do not contort a question into HIVsim simply because the researcher entered through the HIVsim door.

Two failure modes this skill exists to prevent:

1. **Tool-first framing.** *"I have HIVsim, so this must be an HIVsim question."* A question about *predicting individual acquisition risk from a cohort* is a statistical / ML question even though the outcome is HIV. A question about *whether a program increased testing* is a causal-inference question even though the substrate is HIV.

2. **Complexity inflation.** *"The system is complicated, so I need an agent-based model."* The correct combination is not `complicated → ABM`. It is `question complexity + required heterogeneity/interactions + sufficient information to parameterize those mechanisms`. An elaborate simulation with poorly informed inputs is not automatically preferable to a well-scoped simpler analysis.

The correct sequence:

> research question → analytical objective → need for dynamics → need for mechanistic transmission → required representation of heterogeneity / interactions → data / evidence sufficiency → specific tool

Redirection away from HIVsim is a **successful** outcome of intake, not a failure. Correctly recognising that HIVsim is unnecessary saves the researcher months.

## Instructions

1. **Classify the analytical objective first.** Ask the researcher what they are ultimately trying to learn, and place it in one of these categories:

    | Objective | Examples | Likely method family |
    |---|---|---|
    | **Description / estimation** | Current HIV prevalence; ART coverage over time; PrEP use by age × sex | Descriptive statistics; survey analysis |
    | **Association / prediction** | What features predict HIV acquisition? Predict ART discontinuation | Statistics / ML — not HIVsim, even though the outcome is HIV |
    | **Structure / clustering / dimensionality** | Latent behavioural sub-populations; reduce measures to fewer dimensions | Clustering / latent-variable / dimensionality reduction |
    | **Causal estimation from observed data** | Did policy X increase testing? Effect of intervention X on outcome Y | Causal inference / study-design methods before any transmission model |
    | **Forecasting** | Infections in the next 4 weeks; prevalence next year | Time-series if an empirical forecast is sufficient; dynamic transmission only if mechanism matters |
    | **Dynamic mechanistic** | Long-run intervention impact; indirect population effects of a scale-up | Dynamic transmission model — proceed to step 2 |

    Do not skip this step just because HIVsim is the tool the researcher walked in with.

2. **Only if the objective is dynamic-mechanistic, choose the representation.** Compartmental and agent-based transmission models both propagate epidemiological processes through time — the choice between them is not "do we need dynamics?" but "what representation of the population and its interactions does the question require?"

    **Compartmental may be sufficient when:** population- or group-average dynamics are adequate; heterogeneity can be represented with a manageable number of strata; individual histories are not important; explicit contact networks are not important; interventions operate at population / group level; the question does not depend strongly on within-individual correlations; computational simplicity is valuable. Modern compartmental models can be highly sophisticated (age, sex, risk groups, stages, intervention states, time-varying parameters, stochasticity, spatial structure) — do not recommend an ABM just because the system is complicated.

    **Agent-based / individual-based adds value when the answer depends on:**
    - **Individual heterogeneity** attached to the same person across time (age, sex, risk, intervention history, adherence, partnership behaviour)
    - **Correlated characteristics** — combinations that would require an unwieldy number of compartments (eligibility depending simultaneously on age, sex, prior PrEP use, partnership behaviour, testing history)
    - **Individual histories / path dependence** — what happens next depends on what happened previously *to that specific person*
    - **Explicit networks or interactions** — who interacts with whom rather than average mixing (especially for STIs, where partnership structure is often central)
    - **Targeted interventions** defined by combinations of individual characteristics, histories, or relationships
    - **Emergent population effects** from many heterogeneous individual interactions

3. **If an agent-based transmission model is appropriate, check whether HIVsim / STIsim specifically is appropriate.** Inspect actual capabilities, do not assume the mechanism exists. Relevant strengths: HIV/STI transmission over time; sexual partnership / network dynamics; individual demographic and behavioural heterogeneity; HIV testing and diagnosis; ART; PrEP; VMMC; disease progression; targeted intervention strategies; downstream transmission effects. If a required mechanism is missing, determine which of these applies:

    - (a) the mechanism can reasonably be added (route to `extending-stisim` for the upstream-vs-downstream decision);
    - (b) another existing model is more appropriate — name it and say so; or
    - (c) a different analytical method is required — return to step 1 with the refined understanding.

4. **Check data / evidence sufficiency.** Ask: do we have enough information to parameterise the aspects of heterogeneity that are essential to this question? If not, options are: simplify the question; simplify the model; gather or review additional evidence; explicitly represent uncertainty; or conclude that the proposed analysis is not currently well supported. An elaborate simulation with poorly informed inputs is not automatically preferable to a simpler analysis that respects what the evidence supports.

5. **Recognise compound questions.** A single project may require multiple methods. Do not force everything into one bucket. For example:

    > "Which characteristics predict PrEP uptake, and what would transmission look like if future uptake followed those patterns?"

    naturally decomposes into (1) statistical / ML analysis to estimate uptake relationships, (2) translation of those estimates into intervention assumptions, (3) HIVsim to estimate downstream transmission consequences. Produce a multi-method plan when the question warrants it.

6. **Recognise the redirect examples explicitly.** Use these as pattern-matching anchors for the researcher's phrasing:

    - *"What features of this dataset are most predictive of HIV risk?"* → statistical / ML prediction, not HIVsim.
    - *"Can you cluster these behavioural variables into groups?"* → clustering / unsupervised learning, not HIVsim.
    - *"Did this program increase HIV testing?"* → causal inference if observational or experimental data support the estimand; a simulation may complement but should not substitute for identification from data.
    - *"If this testing program were scaled nationally, how would earlier diagnoses and ART initiation affect incidence over 20 years?"* → potentially HIVsim; transmission dynamics and indirect effects now matter.
    - *"If we have 100,000 long-acting PrEP courses annually, how should they be distributed across age / sex / risk groups to maximise population impact?"* → potentially HIVsim; heterogeneous individuals, intervention targeting, sexual networks, and dynamic transmission effects likely required.

7. **Write the methodological assessment into the analysis brief.** Use the output format below. The assessment explains *reasoning*, not just tool name — the researcher should be able to push back on the classification, not just on the recommendation.

## Output format

Return a short methodological assessment suitable for the `Methodological assessment` field of the analysis brief:

```
**Analytical objective:** <one-sentence classification>
**Dynamic transmission required:** <yes / no> — <reasoning>
**Individual-level representation potentially valuable:** <yes / no / n/a> — <reasoning>
**HIVsim / STIsim suitability:** <plausibly strong / not required / unsuitable / requires missing mechanism>
**Evidence gaps:** <what needs to be gathered or reviewed>
**Next step:** <route to a specific downstream skill or workflow>
```

Two example completions:

**Example — plausibly HIVsim:**

> **Analytical objective:** Estimate the population impact of alternative allocations of long-acting PrEP.
> **Dynamic transmission required:** Yes — the outcome includes infections indirectly prevented via changes in onward transmission.
> **Individual-level representation potentially valuable:** Yes — allocation is defined by combinations of individual characteristics and impacts propagate through sexual partnerships.
> **HIVsim suitability:** Plausibly strong; next step is to verify the required PrEP targeting, uptake, persistence, and relevant population structure can be represented.
> **Evidence gaps:** Population-specific uptake and persistence assumptions need review.
> **Next step:** Proceed to evidence discovery and `hiv-interventions` intake.

**Example — not HIVsim:**

> **Analytical objective:** Identify variables predictive of HIV acquisition in an existing cohort.
> **Dynamic transmission required:** No.
> **HIVsim suitability:** This question does not require HIVsim.
> **Likely approach:** Statistical / ML prediction; method selection depends on intended use of the prediction, outcome structure, sample size, and validation requirements.
> **Next step:** Route to an appropriate statistical / ML analysis workflow. HIVsim intake ends here.

## Checks before completion

1. The analytical objective has been classified (one of the six categories), with reasoning.
2. If dynamic-mechanistic, the compartmental-vs-ABM decision is explicit with reasoning.
3. If ABM, HIVsim / STIsim suitability has been assessed against actual capabilities, not assumed.
4. Data / evidence sufficiency has been considered.
5. Compound-question decomposition has been considered if the question spans method families.
6. The methodological assessment is written into the analysis brief with the reasoning, not just the tool name.
7. If the assessment routes the researcher away from HIVsim, that redirection is stated positively — a successful intake outcome, not a failure.
