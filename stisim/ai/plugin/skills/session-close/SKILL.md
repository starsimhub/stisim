---
name: session-close
description: Use at the end of a STIsim work session to produce a handoff summary, verify the analysis brief is current, and surface any unresolved items or uncommitted state.
---

# STIsim session close

## When to use

- The user is about to end a work session and wants to leave a clean handoff.
- Trigger phrases: "wrap up", "closing for today", "call it", "let's stop here", "session done", "handoff".
- Also whenever a substantive block of work has landed and the user wants state captured before context is lost.

## When NOT to use

- Mid-session, when work is still in flight. This skill assumes the session is at a natural pause.
- To close a whole research project (archival, cherry-pick to main). That is a different pattern and belongs in a later, distinct skill.

## Instructions

*Scaffold — procedural content to be written by the lead developer.*

1. Establish what the session covered. Ask the user for a summary if the conversation does not make it obvious.
2. If an `analysis-brief.md` exists in the workspace, verify that its sections are current with what actually happened this session. If it does not exist yet (research-workflow has not shipped, or the user is doing informal work), skip this step.
3. Check for uncommitted changes in the workspace repo. Ask the user how to handle them: commit, WIP branch, or leave with a note. Do not commit unilaterally.
4. Produce a handoff summary — see `references/handoff-template.md` for the format. Offer to append it to `SESSION_LOG.md` at the workspace root, newest first.
5. Surface any items flagged as unresolved — from `analysis-brief.md` if present, from the session dialogue otherwise. State them plainly rather than silently closing.

## How this skill is evaluated

A session closed with this skill should be resumable cold by a different person (or a fresh Claude session) reading only the handoff summary and the workspace state. Judged manually across the two example analyses as they mature.

## Checks before completion

1. Handoff summary written, or explicitly declined by the user.
2. Uncommitted work committed, on a WIP branch, or explicitly noted as intentional.
3. Unresolved items surfaced to the user, not silently accepted.
