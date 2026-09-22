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

1. **Establish what the session covered.** Summarise the session back to the user in 3–5 bullets so they can confirm or correct. If the conversation was long or the arc is unclear, ask directly: *"What are the two or three most important things that happened this session?"* Distinguish work done from work only discussed — both matter but they are recorded differently.

2. **Update `analysis-brief.md` if present.** Check the workspace root for `analysis-brief.md`. If absent, skip to step 3. If present, verify that its *Intended model changes* and *Unresolved decisions* sections reflect this session. Propose the specific edits and ask the user to confirm before writing. Do not restructure the brief — only add, update, or mark items resolved.

3. **Check workspace state.** Run `git status` in the workspace root and report the state. If there are uncommitted changes, ask the user which of the following applies:
   - Commit now with a message the user provides or approves.
   - Move the changes to a WIP branch and commit there.
   - Leave uncommitted and note the intent explicitly in the handoff.

   Do not commit unilaterally; this is a material state change.

4. **Produce a handoff summary.** Fill in the template at `references/handoff-template.md`. Ask the user whether to append the block to `SESSION_LOG.md` at the workspace root (newest at the top; existing entries not modified) or leave the summary in the conversation only.

5. **Surface unresolved items.** Enumerate anything from `analysis-brief.md`'s *Unresolved* section (if the brief exists) or from the session dialogue. State each item plainly. Do not resolve them silently or defer them without recording the deferral in the handoff.

6. **Confirm session is closed.** Restate: the workspace's current state, where the handoff lives, and what to pick up next time. Ask whether anything else needs capturing before the user disconnects.

## How this skill is evaluated

A session closed with this skill should be resumable cold by a different person (or a fresh Claude session) reading only the handoff summary and the workspace state. Judged manually across the two example analyses as they mature.

## Checks before completion

1. Handoff summary written, or explicitly declined by the user.
2. Uncommitted work committed, on a WIP branch, or explicitly noted as intentional.
3. Unresolved items surfaced to the user, not silently accepted.
