---
name: model-primer
description: Use before answering STIsim-specific questions, writing STIsim code, or scoping an analysis with the model. Provides an index and mental model of STIsim's architecture and points at canonical sources.
---

# STIsim model primer

## When to use

- Before answering STIsim-specific questions about module semantics, API shape, or how components compose.
- Before writing or modifying STIsim code — load this first and consult `references/architecture.md`.
- Before scoping an analysis, to check what the model does and does not represent.

## When NOT to use

- General Python or scientific-Python questions.
- Calibration workflow — hand off to the `calib:` plugin.
- HIV or STI epidemiology questions that do not involve the model — use standard sources, not this primer.

## What this skill provides

An index and mental model of STIsim's architecture. It points at canonical sources (source code, docs, tests) and describes how the pieces fit together. It does not duplicate the user guide, the source, or the tests.

## Instructions

1. Read `references/architecture.md` for the mental model of how STIsim's components compose.
2. For a specific area, use `references/canonical-sources.md` to locate the authoritative implementation, docs, and tests.
3. If the task requires deep semantics (parameter meanings, algorithm details, edge behavior), read the canonical source directly rather than paraphrasing from this primer.

## How this skill is evaluated

Design-time evaluation: whether invoking this skill measurably improves Claude's factual answers to STIsim-specific questions relative to not invoking it. Judged against a small held-out set of questions drawn from Example #1 (Zambia HIV recency + PN) as that example matures.

## Checks before completion

Reference skills produce no artifacts. Not applicable.
