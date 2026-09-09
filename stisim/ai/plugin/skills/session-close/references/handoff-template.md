# Handoff template

*Scaffold — content to be written by the lead developer.*

A session handoff summary lets someone opening the workspace cold (a colleague, a fresh Claude session, or you two weeks from now) pick up without re-deriving context from git log alone.

Suggested sections (rearrange, cut, or extend as the pattern earns its keep):

- **Session date** — YYYY-MM-DD.
- **What was done this session** — a few bullets, plain language.
- **Where the analysis stands** — what is working, what is not, what has been ruled out.
- **Next steps** — what to pick up next time, in order.
- **Unresolved** — decisions deferred, questions still open, blocking dependencies.
- **Workspace state** — branch, uncommitted files (or "clean"), any temp scaffolding worth knowing about.

Written where: append to `SESSION_LOG.md` at the workspace root, newest at the top, each entry under a `## <date>` heading. Existing entries are not modified; each session gets a fresh block.

Length target: short. If the block runs longer than the screen, the session either did too many things or the summary is too detailed.
