---
name: project-memory
description: Use during initial project intake to establish a durable memory strategy for a research project, and throughout the project to enforce that important decisions live outside any single agent session. Directs research projects to a deliberate memory mechanism (repo-based, session-indexing, agent-memory provider, or hybrid) rather than letting the current chat session become the de facto project memory.
---

# Project memory

## When to use

- Initial intake for a new research project — before code, before analysis, before the first meaningful session.
- The researcher is starting to feel that "everything important lives in this session" — interrupt and set up durable memory.
- Reviewing a project's memory setup mid-flight because the researcher is switching machines, resuming after a long break, adding a collaborator, or handing off.
- Trigger phrases: "start a new project", "set up this project", "we're just going to keep this session going", "let me paste all the context", "I need to remember X for later", "how do you remember what we did last time", "hand this off to Y".

## When NOT to use

- The user is running a truly one-off analysis that will not continue beyond this session — a quick calculation, a one-shot plot, an exploratory question. Setting up durable memory infrastructure is overhead.
- The user is inside an established session and asks to record one specific fact — that goes into the existing project memory mechanism directly. Invoke this skill only if no mechanism has been set up yet.

## Framing

**A chat session is working memory, not project memory.** A research project may last months or years, span many sessions, move between laptop and VM, and involve multiple collaborators. If the project's real memory is one enormous agent session, none of those transitions work: switching machines is painful, resuming after a break requires re-reading a transcript, sharing with a collaborator is impossible, and context compression can silently discard important detail.

The specific memory technology matters less than satisfying the requirements. A well-chosen mechanism is:

- **Durable** — survives sessions and context resets
- **Project-scoped** — one project's memory does not mix with unrelated work
- **Portable** — moves with the researcher across machines and collaborators
- **Searchable** — findable without pasting the whole project's history into every prompt
- **Provenance-preserving** — you can tell "we decided X" from "the code happens to do X"
- **Inspectable** — humans can read, correct, delete; not opaque embeddings alone
- **Shareable when needed** — with an explicit collaborator-vs-personal distinction
- **Secure** — appropriate for unpublished results, confidential data, credentials, drafts

Project memory sits inside the same sharing theme as `comment-hygiene` and `extending-stisim`: a project whose memory lives only in one person's chat sessions is siloed by construction — invisible to collaborators, to future you on a different machine, and to the next agent session that starts cold.

## Instructions

1. **In initial project intake, raise memory strategy explicitly.** Ask something like:

    > "This project may span many agent sessions and potentially multiple machines or collaborators. What would you like to use for durable project memory? If you already have a system, we can use it. Otherwise I can help you pick between a project-local approach, a session-indexing system such as Funes, an agent-memory provider (e.g. Hermes Agent supports several backends), or a hybrid."

    Do not require one specific provider. If the researcher declines a specialised system, establish at minimum a version-controlled project-memory structure in the repository (see step 2).

2. **The default fallback: repo-based markdown memory.** For projects that decline a specialised system, propose a small set of committed markdown files. Do not create all of these mechanically — a single concise file may suffice. Only create what the project needs:

    - `PROJECT.md` — stable project purpose and scope
    - `DECISIONS.md` — methodological or software decisions with rationale
    - `STATUS.md` — current state and immediate next steps
    - `ASSUMPTIONS.md` — scientific assumptions that should stay visible
    - `DATA_SOURCES.md` — provenance and notes about important datasets

    Transparent, portable, version-controlled, available anywhere the repo is cloned. Weakness: requires active curation and does not automatically make old sessions searchable — pair with a session-indexing layer if that matters (step 3).

3. **Consider the hybrid pattern for substantial projects: curated project memory + searchable session history.** The curated layer holds the project's durable understanding; the session layer allows deeper recall of why a decision was made or what was tried weeks earlier. Treat session memory as a retrieval layer over history, not as the authoritative record — important decisions should still be promoted into curated memory when they stabilise.

4. **Teach the memory hierarchy explicitly.** Not all context deserves the same permanence:

    | Layer | What lives here | Examples |
    |---|---|---|
    | **1. Authoritative artifacts** | Source of truth | Code, config, input data, manuscript, outputs, protocols, committed docs |
    | **2. Curated project memory** | Durable context not readable from the artifacts | Why a methodological choice was made, rejected alternatives, conventions, assumptions, unresolved questions, collaborator decisions |
    | **3. Searchable session history** | Detailed working context | Debugging history, exploratory ideas, abandoned approaches, undistilled rationale |
    | **4. Session working memory** | Disposable | Temporary plans, intermediate thoughts, scratch calculations, hypotheses under test |

    Memory should point at layer 1 rather than duplicate it. Layer 4 should not be promoted wholesale into layer 2.

5. **Promote important information upward.** When a session concludes with a project decision — e.g. "we will calibrate ART coverage directly to UNAIDS estimates rather than treating ART initiation rates as independently identified" — record the conclusion in curated memory (concise, dated if useful, with pointers to the relevant artifacts). The 40-turn conversation that led to it can stay in session history.

6. **Address cross-machine, collaborator, and security questions explicitly.**
    - **Machines:** if memory must follow the researcher across machines, use a deliberate sync mechanism (git for markdown, a shared project store, or a memory provider designed for cross-machine access). Do not discover halfway through a project that all useful memory lives in a hidden directory on a terminated VM.
    - **Collaborators:** ask *"should this be shared among collaborators or is it your personal working memory?"*. Some information describes the project (methodology, assumptions, provenance, status) and should be shared; some is personal (individual scratch reasoning, task organisation, incomplete thoughts) and need not be. Configure the memory system accordingly.
    - **Security:** research memory may contain unpublished results, confidential data, credentials, or drafts. Assess a cloud-hosted memory service against the project's data-handling requirements before adopting it. Never rely on automated secret filtering as the sole protection for sensitive material.

7. **What belongs in durable memory.** Good candidates: major modelling decisions; important scientific assumptions; data inclusion / exclusion decisions; provenance of unusual inputs; rationale for non-obvious implementation choices; calibration strategy; important rejected approaches and why they were rejected; known limitations; recurring debugging discoveries; project-specific conventions; collaborator decisions; current unresolved questions.

    **What does not:** every command run, every failed test, every temporary debugging theory, generic facts easy to rediscover, large copies of source files, verbose per-session summaries, conversational filler. Memory should reduce future work, not create another corpus to wade through.

8. **Anti-pattern: the eight-million-token session.** If a session is becoming the project's de facto memory, interrupt: recommend moving durable information into the chosen memory mechanism and starting fresh. Symptoms include: reasoning depends on inaccessible earlier turns, switching machines is painful, sharing is impossible, context compression is imminent.

9. **Never use shared-package source comments as project memory.** If a downstream research project discovered an important project-specific fact, put it in project memory. If it revealed a general property of the shared software, update the software documentation or implementation. Full rule in `comment-hygiene`.

10. **Periodically maintain durable memory.** Long-running projects need memory maintenance the way they need code maintenance: review for stale conclusions, reversed decisions, duplicated entries, resolved questions, assumptions since anchored to data. Where the memory system supports it, mark decisions as superseded rather than silently rewriting history.

11. **Every significant session should end recoverable.** Before closing a substantive session, verify: important decisions recorded in durable memory; current status recoverable; unresolved questions captured; session-indexing synchronised if applicable; a competent new agent can resume from durable memory plus the repo without needing the preceding conversation pasted in. Per-session mechanics in `session-close`; this skill is about the memory infrastructure `session-close` writes into.

## Checks before completion (initial intake)

By the end of project setup, the researcher and agent should have agreed on:

1. What system will hold curated project memory.
2. Whether session histories will also be indexed / searchable.
3. Whether memory needs to work across machines.
4. Whether memory is personal or collaborative.
5. Where sensitive information is allowed to be stored.
6. How important decisions will be promoted into durable memory.
7. How a new agent session should recover context.

Record the agreed choice somewhere durable in the project (typically in `PROJECT.md` or its equivalent), so subsequent agents can see how project memory is supposed to work without having to re-elicit the strategy from the researcher.

## The standard a good memory setup meets

**A competent new agent should be able to resume the project from durable memory plus the repository, without needing the preceding conversation pasted into its prompt.** If that standard fails on any given day, the memory system needs attention before the next substantive session.
