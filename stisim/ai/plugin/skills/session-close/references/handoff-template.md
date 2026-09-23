# Handoff template

A handoff entry lets someone opening the workspace cold — a colleague, a fresh Claude session, or the same user two weeks from now — pick up without re-deriving context from `git log` alone.

Append entries to `SESSION_LOG.md` at the workspace root, newest at the top, each under a `## YYYY-MM-DD` heading. Do not modify existing entries; each session gets its own block.

## Template

```markdown
## YYYY-MM-DD

**Done this session:**
- <bullet>
- <bullet>

**State of the analysis:**
- Working: <what is in a good state>
- Not yet: <what is planned but not done, or partly done>
- Ruled out: <what was tried and abandoned — only include when relevant>

**Next steps:**
1. <first thing to pick up>
2. <second thing>

**Unresolved:**
- <open question, deferred decision, blocking dependency — or the single
  word "None." if the session left no open items>

**Workspace state:**
- Branch: <name>
- Uncommitted: <one-line summary, or "Clean.">
```

## Style

Fits on one screen. If a handoff runs longer, either the session did too many things or the summary is over-detailed. Plain sentences over bullet paragraphs. Names, not "the thing we discussed"; concrete file paths and commit shas where they help a cold reader.
