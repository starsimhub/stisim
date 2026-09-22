---
name: comment-hygiene
description: Use when writing or editing docstrings and inline comments in shared library code (stisim, starsim, or any package used by people outside the current project). Enforces the rule that library comments explain the software itself, not the history of the project that motivated the change or the agent's working context.
---

# Comment hygiene for shared library code

## When to use

- Writing or editing docstrings, inline comments, or module-level narration in shared library code (stisim, starsim, or any other `pip install -e`-shared package).
- Reviewing an inherited file for comments that read like agent scratch space or downstream-project memory.
- Preparing a PR against a shared package — check that no comments carry references only meaningful in the originating conversation.
- Trigger phrases: "add a docstring", "document this function", "add a comment explaining", "clean up the comments", "why is there a comment about X", "explain what this does".

## When NOT to use

- Writing comments in downstream project or experiment code, where project-specific references ("exp 015", "the ANC screening handoff") are stable and meaningful within that repo's scope. This skill is scoped to shared library code specifically.
- Writing test names or test docstrings that reference a specific historical scenario — a regression test can and should preserve the concrete failing case; production code should stick to the underlying invariant. See instruction 6.

## Framing

Comments in shared library code have a different audience than comments in project code. A stisim docstring is read by every stisim user, forever. A comment in an experiment folder is read by you and maybe one collaborator, all of whom share the same project context.

The rule for shared code: **explain the software, not the history of the project that motivated the change.** A comment like *"Special case added to fix Experiment 5"* makes sense to the agent that wrote it during one downstream analysis, but a future stisim maintainer has no idea what "Experiment 5" refers to, whether it is still relevant, or why it mattered.

AI makes this failure mode especially common. Agents run inside task-specific conversations, so their default is to describe *why they are making the change right now* — the current task, the current experiment, the user who requested it. That context is invisible to everyone who reads the code later. Library comments must be written for a maintainer who has no access to the originating conversation.

This is the same siloing dynamic that motivates `extending-stisim` and `editable-dep-hygiene`, applied to the documentation layer: a comment only the originating agent can decode is unshared knowledge, even though it lives in shared code.

Three kinds of information a comment might carry — classify before writing:

| Kind | Belongs where |
|---|---|
| Intrinsic to the software (units, invariants, non-obvious contracts) | Docstring or inline comment, if useful |
| Historical rationale for a durable design decision | Concise comment, or a linked issue / PR / ADR / citation |
| Context for the current project, analysis, debugging session, or conversation | **Project memory or issue tracker — never in shared source** |

## Instructions

1. **Docstrings describe the public contract, not the implementation.** Include: what the function does, non-obvious arguments, non-obvious return value, side effects, invariants, units, exceptions callers must handle. Exclude: line-by-line narration; restating information already obvious from the name, type annotations, or code; long prose written because an agent can generate it. Prefer short and precise over comprehensive.

2. **Inline comments explain the WHY, not the WHAT.** Good reasons for a comment: a non-obvious mathematical or epidemiological assumption; a subtle invariant; a workaround for a clearly identified external limitation; a compatibility constraint; an implementation choice where an obvious alternative would be wrong; the provenance of an important constant, with a durable citation. If a comment restates the code, delete it.

3. **Never encode downstream-project memory in shared library code.** Explicit anti-patterns to strip: "Experiment 5", "the canonical case", "Cassie's analysis", "the Zimbabwe run we're debugging", "the notebook downstream", "this fixes the current task", "we changed this because the user requested it". If a change was motivated by a downstream application, **translate the motivation into the general software reason** for the change:

    - Bad: `# Special case added to fix Experiment 5.`
    - Better: `# Allow coverage targets to change representation over time, since historical datasets may switch between counts and proportions.`

    The second explains a general capability and remains meaningful to any future user.

4. **Project-specific memory goes outside shared source code.** Inspect the downstream repo first for its convention — agent memory files, analysis README, experiment SUMMARY, issue tracker, PR description, commit message, or a project-level decisions file — and use it. If the project has an agent-memory system, update that memory when project-specific context needs to persist across sessions. Do not invent a new mechanism, and do not encode the memory into library comments.

5. **When editing a file, notice obviously problematic nearby comments and clean them up if low-risk.** Stale, misleading, excessively verbose, tied to an inaccessible downstream project, or functioning as agent scratch memory — remove or rewrite. Do not aggressively rewrite unrelated documentation, and do not add comments merely to make the diff look more documented.

6. **Prefer tests over comments for regression cases.** If a behavior exists because of a specific past failure, a regression test preserves the concrete case more reliably than a comment. The production code documents the underlying invariant; the test documents the case that once failed. Test names or test docstrings may reference a specific historical scenario — that scoping is legitimate; the same language does not belong in the core module.

## Checks before completion

1. No comment in shared library code references a downstream experiment, notebook, analysis, or conversation by name or by an unnamed pointer ("the current task", "this fix", "the user requested").
2. No docstring narrates the implementation line by line or restates what the function name and types already convey.
3. Every comment that survived exists because it explains WHY something non-obvious is happening; anything that only restated WHAT the code does has been deleted.
4. If a change was motivated by a downstream analysis, that motivation was recorded in the project's designated memory location — not embedded in the library comment.

## Writing style

Comments and docstrings in shared library code should be concise, factual, durable, and understandable without access to the agent's conversation history. Write for a future maintainer who has never met you, never read the originating chat, and has no idea what project motivated the change.
