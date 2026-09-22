# stisim.ai — AI-assisted STIsim research

A Claude Code plugin that ships alongside the `stisim` Python package. Its
goal is to install guardrails, checkpoints, and workflow discipline around
STIsim / HIVsim research so that model-assisted analyses are more
defensible, more reproducible, and less prone to known failure modes. The
value proposition is **guardrails, not speed** — the plugin may sometimes
slow an experienced user relative to vanilla Claude Code. Better work,
not faster work.

## Skill inventory

Fourteen skills, in three functional clusters.

### Analysis workflow — scaffold and structure a research project

| Skill | Purpose |
|---|---|
| `analysis-intake` | Grilling-style design-tree intake that turns a vague research idea into a provisional analysis specification. Adapts Matt Pocock's [`grilling`](https://github.com/mattpocock/skills/blob/main/skills/productivity/grilling/SKILL.md) with facts-vs-decisions split and `blocked-on-evidence` as a third branch state. |
| `analysis-selector` | Runs very early in intake, before any tool is assumed. Classifies the research question into one of six analytical objectives (description / association / structure / causal / forecasting / dynamic-mechanistic) and only routes to HIVsim if dynamic mechanistic transmission is actually required. Redirection away from HIVsim is a successful outcome. |
| `project-memory` | Establishes a durable memory strategy at project intake (repo-based, session-indexing, agent-memory provider, or hybrid) so a multi-month project can survive session ends, machine switches, and collaborator handoffs. A chat session is working memory, not project memory. |
| `session-close` | Handoff summariser at session end, feeding into whichever memory mechanism `project-memory` established. |

### Model authoring and calibration — invoked once the analysis is scoped

| Skill | Purpose |
|---|---|
| `model-primer` | Reference-layer architectural index of stisim — assembly, disease hierarchy, networks, interventions, connectors, analyzers, demographics, care-seeking, timestep discipline. Includes the `calibration-knobs.md` reference for what stisim actually exposes. |
| `model-writer` | Composes a new Sim from a scoped research question. |
| `hiv-interventions` | Design-interview for the ART / testing / VMMC / PrEP set, with per-intervention data-source references. |
| `network-data` | DHS-focused sexual-network calibration input. |
| `calibration-strategy` | HIVsim / STIsim-specific what/why/whether-to-calibrate: parameter classification (data-informed / literature-fixed / country-specific-uncertain / behavioural / tuning / intervention-assumption), target-vs-knob distinction, ART/testing structural dependency, failure diagnosis. Delegates algorithm / sampler / likelihood / diagnostics to the `calib:*` plugin. |

### Software quality — cross-cutting, unified by "share your work"

All in this cluster reinforce the theme *share your work rather than
accumulate private forks of shared code or private caches of project
knowledge* — a pattern that becomes especially easy to fall into when
code generation is cheap and code review is not.

| Skill | Purpose |
|---|---|
| `extending-stisim` | Before subclassing / monkey-patching stisim, classify the change as (a) fix, (b) opt-in research knob, or (c) project-specific data preprocessing, and enforce PR-upstream for fixes, no-op default for knobs, cleanup of downstream artifacts once their upstream fix lands, and file-an-issue as the fallback. |
| `editable-dep-hygiene` | Any edit to an editable `pip install -e` dependency gets committed and PR'd immediately; before a version bump or branch switch, the dep checkout is verified clean of unmerged local work. |
| `comment-hygiene` | Shared-library docstrings and comments explain the software itself, not project-specific memory. Strips anti-patterns like "Experiment 5", "the current task", "we changed this because the user requested it". |
| `result-extraction` | Directs to `ss.Result` / `ss.Results` methods (`annualize`, `resample`, `to_df`) over hand-rolled `df.groupby('year').mean() / .sum()`, which silently mishandles the flow-vs-stock distinction. |
| `writing-tests` | Prefer a small number of scientifically meaningful tests over comprehensive enumeration. Every proposed test should have a one-sentence answer to "what meaningful bug would this catch?"; if the answer is "confirms a value we assigned still has that value", don't add it. |

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
        ├── analysis-intake/
        ├── analysis-selector/
        ├── calibration-strategy/
        ├── comment-hygiene/
        ├── editable-dep-hygiene/
        ├── extending-stisim/
        ├── hiv-interventions/
        ├── model-primer/
        ├── model-writer/
        ├── network-data/
        ├── project-memory/
        ├── result-extraction/
        ├── session-close/
        └── writing-tests/
```

Each skill directory contains `SKILL.md`; skills with deeper reference
material also carry a `references/` subdirectory.

## Activation

### From a released stisim (pip)

```
pip install stisim
python -m stisim.ai install
```

### From an editable clone (contributors, and anyone tracking `main` or a feature branch)

```
git clone git@github.com:starsimhub/stisim.git
cd stisim
pip install -e .
python -m stisim.ai install
```

`stisim.ai install` resolves the plugin path from the installed `stisim`
package location, so an editable install activates the plugin *from your
git checkout* — any local edits to `stisim/ai/plugin/skills/` are picked
up immediately on the next Claude Code reload, without a reinstall.
This is the workflow to use if you are developing or editing the skills
themselves.

If you `git pull` on the editable stisim checkout, the plugin content
updates on the next Claude Code reload; no reinstall needed. See the
`editable-dep-hygiene` skill for the git protocol around editable
dependencies.

### After install (either path)

Reload the Claude Code session:
- VS Code / Positron: Command Palette → *Developer: Reload Window*
- CLI: exit and relaunch `claude`

`/stisim:*` skills are then available in every Claude Code session.

The bootstrap edits `~/.claude/settings.json` to add `stisim@stisim-local`
under `enabledPlugins` and register a directory-sourced marketplace at the
stisim plugin path. The write is atomic, idempotent, and preserves all
other settings. `python -m stisim.ai uninstall` reverses it;
`python -m stisim.ai status` reports current registration.

## Companion plugins

`stisim.ai` is intentionally narrow — it covers HIVsim / STIsim
specifically. The plugins below cover adjacent concerns (generic
calibration, generic starsim / disease modelling, project memory,
IDM-wide engineering standards) and are designed to compose. We
strongly encourage installing them alongside `stisim.ai`.

- **[`calib`](https://github.com/InstituteforDiseaseModeling/calib-plugin)**
  — generic calibration machinery: algorithm choice, prior predictive,
  re-identification, workflow sequencing, method selection,
  likelihood design, plotting. `calibration-strategy` decides *what*
  to calibrate; `calib:*` skills handle *how*. Install:
  `/plugin marketplace add https://github.com/InstituteforDiseaseModeling/calib-plugin`.

- **[`starsim_ai`](https://github.com/starsimhub/starsim_ai)** —
  three Claude Code plugins covering the layer under stisim:
  `starsim-ai` (Starsim + Sciris MCP tools and modelling skills),
  `disease-modeling` (general disease-modelling skills that apply
  beyond HIVsim/STIsim), and `project-improver` (engineering-quality
  review). Install:
  `/plugin marketplace add https://github.com/starsimhub/starsim_ai`.

- **[`canonize`](https://github.com/emiliasimmons/canonize)** — agent
  skills for durable decision capture in computational modelling
  projects. *"Sources feed the wiki. Decisions are internal sources.
  Collaborators browse the wiki, not the sources."* Complements
  `project-memory` with a concrete, opinionated implementation of the
  curated-project-memory layer.

- **[`idm_standards`](https://github.com/InstituteforDiseaseModeling/idm_standards)**
  — IDM's central hub for software-quality standards, engineering
  practice, style, and documentation. Covers what `comment-hygiene`,
  `writing-tests`, and the other software-quality skills touch, but
  from a broader institutional standards perspective. Install:
  `/plugin marketplace add https://github.com/InstituteforDiseaseModeling/idm_standards`.

## AI collaboration

This subpackage was scaffolded in a collaborative design session with
Claude Code (Anthropic Claude Opus 4.7). The plugin skeleton, bootstrap
CLI, initial skill scaffolding, and this README were drafted by Claude
under human review. All architectural and scoping decisions —
package layout, invariants, install mechanism, roadmap ordering, and what
belongs in MVP versus later — were made by humans on the STIsim team.
Skill body content is authored by AI and humans.
