# Analysis brief — living template

An analysis brief captures the current provisional specification of an HIV Sim analysis. It is **living** — fields are added, changed, and moved between statuses throughout the project. It is not a form to complete once. It is the durable output of `analysis-intake` and the input to every subsequent skill.

## Status markers

Every field carries a status, marked explicitly:

- **Decided** — the researcher has confirmed, in writing, that this is the answer.
- **Provisional** — best current guess or recommendation; may change with evidence or further discussion. Recommendations from the agent live here until confirmed.
- **Needs evidence** — cannot be answered until literature / data / code investigation completes. Should always have a concrete investigation task attached.
- **Needs researcher decision** — the evidence is in and the choice now requires scientific judgment.

Do not invent values for unresolved fields. Leaving a field open is more useful than filling it with a placeholder that later gets treated as if it were settled.

## Fields

### Scientific frame

- **Provisional research question** — one paragraph, plain language
- **Decision / scientific purpose** — what the analysis is intended to inform
- **Setting** — country, region, or explicitly-theoretical
- **Population** — who the analysis is about, including any sub-populations
- **Intervention(s)** — what is being modelled as the change
- **Comparator(s)** — what the intervention is compared against
- **Outcomes** — primary and secondary; include distributional / equity outcomes if relevant
- **Time horizon** — from when to when
- **Scenarios** — the specific comparisons the analysis will run

### Evidence

- **Known evidence / data** — what is already in hand, with provenance
- **Evidence / data still needed** — what will be gathered before or during analysis
- **Questions blocked on research** — each with a concrete investigation task
- **Important uncertainties** — the ones expected to matter for the conclusion

### Constraints

- **Deliverables** — figures, manuscript analyses, presentations, policy-facing outputs, methods docs, reusable code, calibrated models
- **Timeline** — final deadline; intermediate deadlines (meetings, presentations, manuscript milestones, policy decision points); whether preliminary results are needed
- **Stakeholders** — who will consume, review, or make decisions using the work
- **Compute environment** — laptop / VM / cluster / cloud / multi-machine
- **Collaborators** — individual or collaborative; if collaborative, who and in what role

### Working style

- **Researcher HIV Sim experience** — first analysis / some experience / expert
- **Desired guidance level** — explicit-every-decision / surface-only-consequential / minimal
- **Memory strategy** — what holds curated project memory (per the `project-memory` skill)
- **Repository / git strategy** — main repo, working branch, PR flow

### Open items

- **Open design decisions** — with recommended answers and reasoning where available
- **Superseded framings** — original questions or scoping that were revised after evidence; preserve so the evolution of the research is part of the record

## Status transitions to expect

- **Provisional** → **Decided** when the researcher confirms explicitly.
- **Needs evidence** → **Needs researcher decision** when investigation returns.
- **Needs researcher decision** → **Decided** when the researcher chooses.
- **Decided** → **Superseded** when a later evidence-informed reframing changes the answer. Do not delete — mark superseded, preserve the prior value, and record the reason for the change.

## Where the brief lives

The brief lives in the project's chosen durable-memory location (per `project-memory`). For repo-based memory, `analysis-brief.md` at the workspace root is the default. For projects using a session-indexing or agent-memory provider, the brief still needs a durable canonical location — the memory provider is a retrieval layer, not the source of truth.
