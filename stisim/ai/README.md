# stisim.ai — AI-assisted STIsim research (WIP)

**Status:** early scaffold, actively iterating. Structure, skill inventory,
packaging and roadmap are all subject to change during this branch.

## What this is

A Claude Code plugin that ships alongside the `stisim` Python package. Its
goal is to install guardrails, checkpoints, and workflow discipline around
STIsim research so that model-assisted analyses are more defensible, more
reproducible, and less prone to known failure modes. The value proposition
is **guardrails, not speed** — the plugin may sometimes slow an experienced
user relative to vanilla Claude Code. Better work, not faster work.

## Proposed structure

Three layers under a `stisim:` plugin namespace:

- **Meta layer** — orchestration and process. Intake (`getting-started`),
  the research workflow, session hygiene (`session-close`), and topic
  teach-and-quiz (`model-exam`).
- **Specialist layer** — invoked by the research workflow at each step.
  `model-writer`, `extend-model`, `scenario-runner`, `results-analyst`,
  and more.
- **Reference layer** — loaded on demand. Architecture index
  (`model-primer`), plus per-topic references as they earn their keep.

Three skills are scaffolded so far:

| Skill | Layer | Status |
|---|---|---|
| `stisim:model-primer` | Reference | Scaffolded; `references/canonical-sources.md` populated, `references/architecture.md` filled in |
| `stisim:session-close` | Meta | Scaffolded; procedural steps and handoff template filled in |
| `stisim:model-writer` | Specialist | Scaffolded; procedural steps and authoring checklist filled in |

## Directory layout

```
stisim/ai/
├── __init__.py
├── __main__.py                     # python -m stisim.ai install / uninstall / status
├── README.md                       # this file
└── plugin/                         # the Claude Code plugin tree
    ├── .claude-plugin/
    │   ├── plugin.json             # plugin manifest
    │   └── marketplace.json        # local marketplace descriptor
    └── skills/
        ├── model-primer/
        │   ├── SKILL.md
        │   └── references/
        ├── model-writer/
        │   ├── SKILL.md
        │   └── references/
        └── session-close/
            ├── SKILL.md
            └── references/
```

## Activation

```
pip install stisim
python -m stisim.ai install
```

Then reload the Claude Code session:
- VS Code / Positron: Command Palette → *Developer: Reload Window*
- CLI: exit and relaunch `claude`

`/stisim:*` skills are then available in every Claude Code session.

The bootstrap edits `~/.claude/settings.json` to add `stisim@stisim-local`
under `enabledPlugins` and register a directory-sourced marketplace at the
pip-installed plugin path. The write is atomic, idempotent, and preserves
all other settings. `python -m stisim.ai uninstall` reverses it;
`python -m stisim.ai status` reports current registration.


## AI collaboration

This subpackage was scaffolded in a collaborative design session with
Claude Code (Anthropic Claude Opus 4.7). The plugin skeleton, bootstrap
CLI, initial skill scaffolding, and this README were drafted by Claude
under human review. All architectural and scoping decisions —
package layout, invariants, install mechanism, roadmap ordering, and what
belongs in MVP versus later — were made by humans on the STIsim team.
Skill body content is authored by AI and humans; areas that still contain
AI-drafted placeholder text are marked as such in the files.
