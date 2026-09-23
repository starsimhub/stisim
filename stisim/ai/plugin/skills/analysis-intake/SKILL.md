---
name: analysis-intake
description: Use at the start of any new HIVsim analysis to convert an initially vague research idea into a provisional analysis specification through an iterative, grilling-style design-tree conversation. Investigates facts the agent can discover (literature, data, HIVsim capabilities); asks the researcher only for decisions that require their scientific judgment. Establishes working style, constraints, and durable memory before code is written.
---

# HIVsim analysis intake

## When to use

- User is starting a new HIVsim analysis and no analysis brief exists yet.
- User is entering a project via any of three paths: **question-first** ("how should I distribute long-acting PrEP…"), **data-first** ("we have new efficacy data — can we do something with it?"), or **tool-first** ("I'd like to use HIVsim to look at X").
- Trigger phrases: "start a new HIVsim analysis", "help me set up a project", "we have new data on X, what should we do with it?", "I'd like to use HIVsim", "where do I start with this analysis", "scope this study with me".

## When NOT to use

- Analysis brief already exists — hand off to `model-writer`, the calib plugin, or the appropriate downstream skill.
- User is resuming an established project — that is session resumption, not intake. See `session-close` for the handoff conventions.
- The question is a genuinely one-off calculation with no follow-on work planned.

## Framing

**Grill decisions; investigate facts.** This skill adapts Matt Pocock's [`grilling`](https://github.com/mattpocock/skills/blob/main/skills/productivity/grilling/SKILL.md) design-tree pattern for scientific research intake. Model the analysis as a **design tree** — each decision branches into the decisions that depend on it. Work the tree in **rounds**. At any moment the **frontier** is the set of decisions whose prerequisites are settled enough to discuss now. Recompute the frontier after every round.

Two adaptations for research:

1. **Facts vs. decisions is a stronger split than in generic grilling.** When a frontier question needs a fact — a published efficacy estimate, whether HIVsim has a mechanism, when a policy changed, what enrolled populations a trial had, what data exist for a country — do not ask the researcher. Investigate: dispatch sub-agents, read the code, search the literature, look at the data. Only *scientific judgment* is on the researcher's plate. Their time is expensive; agent time is not.

2. **Some branches must block on evidence, not on the researcher.** In generic grilling, the goal is to exhaust the tree before acting. Research works differently: some branches cannot be resolved until we investigate the world. Each branch has three states — **settled**, **open**, **blocked on evidence** — and intake progresses when the frontier of *currently-answerable-with-researcher-input* questions is exhausted, even if evidence-blocked branches remain. The typical arc: initial framing → literature / data discovery → evidence extraction → revised framing → model specification → analysis.

The intake should feel like a conversation with a strong scientific collaborator, not a form. Rigorous without being exhausting — only ask a question when its answer will change something.

## Entry paths — detect from the opening

**Question-first** ("How should I distribute long-acting PrEP among sub-populations in my country?"): start with the question, work outward — decision context, populations, intervention definition, comparator, outcomes, time horizon, equity dimensions, constraints, evidence. Do not immediately translate the first sentence into a simulation specification.

**Data-first** ("We have new efficacy / uptake / adherence data — can we do something with it?"): start with what the evidence measures, then explore scientific implications — what mechanism it affects, whether HIVsim represents that mechanism, what existing assumptions it might replace, whether the evidence generalizes to the target setting, what decision becomes answerable because these data now exist. Do not assume interesting data warrants a simulation analysis.

**Tool-first** ("I'd like to use HIVsim"): help discover the research question. Ask what setting, intervention, policy problem, or scientific uncertainty motivated the interest. Explain briefly what HIVsim is good at if useful. Do not immediately start collecting configuration parameters.

## Instructions

1. **Before the first round, assess HIVsim experience and desired guidance level.** Ask enough to establish: first HIVsim analysis? Familiar with agent-based / individual-based models? How much explanation do they want? Do they want to make every modeling decision explicitly, or would they prefer the agent to recommend conventional choices and surface only consequential decisions? If it is a first HIVsim analysis, offer a short crash course — do not force it. Experienced users should be able to say "skip the background". Record the desired guidance level so later skills do not re-ask or re-explain.

2. **Run rounds of 2–5 substantive questions.** Do not dump the whole frontier if it contains fifteen questions — group related questions and prioritize the highest information value. Follow the round format:

    ```
    ❓ **Q1** — **<question title>**: <question body, may include options>

    ➡️ <recommended answer with reasoning>

    ---

    ❓ **Q2** — **<question title>**: …
    ```

    After each round: update the tree, mark decisions as **settled / open / blocked on evidence**, record decisions in the living analysis brief, recompute the frontier, ask the next round.

3. **Investigate facts yourself.** When the frontier includes a factual question — what ART data exist for this country, what a trial's enrolled populations were, what a paper reported, when a policy changed, whether HIVsim supports a mechanism, how a mechanism is currently represented — mark the branch **blocked on evidence** in the tree, dispatch a sub-agent or read the code / literature / data, and continue asking the *other* frontier questions in parallel. Bring evidence back to the researcher when a decision is now unblocked.

4. **Once the provisional research question, populations, and decision purpose are settled enough to classify, invoke `analysis-selector` before any HIVsim-specific configuration.** This is the point at which the intake either confirms HIVsim is the right tool or redirects the researcher to statistics / ML / causal inference / a compartmental model / another method. Record the methodological assessment in the analysis brief. Redirection away from HIVsim is a successful intake outcome, not a failure — do not skip this step because HIVsim is the tool the researcher walked in with.

5. **Recommendations should expose reasoning, not just answers.** Better:

    > **Recommendation:** Model allocation as a constrained number of PrEP initiations rather than unconstrained coverage percentages, because your stated question is about allocating a fixed program supply. We should verify whether the available program data support that representation.

    …than "counts or coverage?". The reasoning lets the researcher push back on the premise, not just the choice.

6. **Only ask questions whose answers will change something.** Skip a question if its answer will not change the research question, analysis design, evidence needs, HIVsim configuration, or project execution. Grilling means rigorous, not exhaustive. Use the researcher's answers to prune irrelevant branches.

7. **Establish practical constraints once the scientific shape is clear enough to make them meaningful.** Timeline (final and intermediate deadlines, meetings, presentations, manuscripts, policy decision points); deliverables (figures, manuscript analyses, presentations, policy-facing outputs, methods documentation, reusable code, calibrated models); stakeholders (who consumes, reviews, or decides using the work); compute environment (laptop, VM, cluster, cloud, multi-machine); collaborators (individual or collaborative — if collaborative, invoke `project-memory` and git-hygiene skills immediately).

8. **Hand off to `project-memory` for durable-memory strategy** once the project's substantive shape is clear. A multi-month HIVsim analysis should not depend on one enormous agent conversation. If the project is collaborative or likely to span machines, treat memory-strategy as an early design decision, not a cleanup task.

9. **Maintain a living analysis brief.** The durable output of intake is a `references/analysis-brief-template.md`-shaped document (see references), not a completed questionnaire. Every field carries a status: **Decided / Provisional / Needs evidence / Needs researcher decision**. Do not invent values for unresolved fields — leaving a field open is more useful than filling it with a placeholder that later gets treated as if it were settled.

10. **Reframe the research question after evidence lands.** Return to the initial framing and ask: did the evidence change what is answerable? Did we discover important populations we had omitted? Are originally proposed distinctions unsupported by data? Did we discover an important comparator? Are the outcomes still appropriate? Does HIVsim represent the required mechanisms? Should the scope expand or contract? Preserve both the original and revised question in the analysis brief when the change is scientifically meaningful — the evolution of the question is part of the research record, not something to overwrite.

## Stopping condition — initial intake

Do not require the design tree to be empty before useful research begins. Initial intake is complete when:

- the provisional research question is clear enough to guide evidence discovery,
- the major purpose and scope are understood,
- immediate evidence / data needs are identifiable,
- the researcher's desired working style is understood,
- and the project has enough operational structure to proceed safely.

At that point, summarize the current shared understanding, explicitly identify unresolved branches (which are still **open** vs. **blocked on evidence**), and move into literature / data discovery. Later, after evidence lands, return to the tree and continue grilling the branches that are now unblocked.

## Checks before completion (initial intake)

1. A living analysis brief exists in the project's chosen durable-memory location.
2. Every field in the brief has an explicit status: **Decided / Provisional / Needs evidence / Needs researcher decision**.
3. HIVsim experience and desired guidance level are recorded so downstream skills do not re-elicit.
4. `project-memory` has been consulted; the memory strategy is recorded.
5. If collaborative, git-hygiene and shared-memory decisions are recorded.
6. Evidence-blocked branches have concrete investigation tasks queued.
7. The next-round frontier is identified and either scheduled or handed off to the evidence-gathering skill.

## The overall principle

**Grill decisions; investigate facts.** Use the researcher's time for the decisions that require their scientific judgment. Use agents and tools to discover everything that can reasonably be discovered independently. The goal is not to finish intake — it is to progressively build a shared, evidence-informed specification of a good research analysis.
