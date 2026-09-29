# AI plugin

STIsim ships with `stisim.ai`, a [Claude Code](https://claude.com/claude-code) plugin that adds STIsim- and HIVsim-specific skills to Claude Code. Its goal is guardrails rather than speed: the skills add checkpoints and workflow discipline around research with STIsim so that AI-assisted analyses are more defensible and reproducible, and less prone to known failure modes.

## Installation

The plugin is included in the `stisim` package, so there is nothing extra to download. After installing STIsim, register the plugin with Claude Code:

```sh
pip install stisim
python -m stisim.ai install
```

Then reload Claude Code so it picks up the plugin:

- **CLI:** exit and relaunch `claude`.
- **VS Code / Positron:** Command Palette → *Developer: Reload Window*.

The installer adds the plugin to `~/.claude/settings.json` (as `stisim@stisim-local`, via a local directory-based marketplace pointing at your installed `stisim` package). It preserves all other settings and is safe to run more than once. To check or undo the registration:

```sh
python -m stisim.ai status     # Check whether the plugin is registered
python -m stisim.ai uninstall  # Remove the registration
```

### Editable installs

If you installed STIsim from a clone with `pip install -e .`, the plugin is loaded directly from your checkout. Edits to the skills under `stisim/ai/plugin/skills/`, or updates from `git pull`, take effect the next time Claude Code reloads, with no need to reinstall.

## Usage

Once installed, the skills are available in every Claude Code session. Claude will invoke them automatically when a task matches a skill's description (for example, asking it to set up a new HIV model for a country will trigger `model-writer`). You can also invoke a skill explicitly by typing `/stisim:<skill-name>`, for example:

```
/stisim:analysis-intake
```

A typical workflow for a new analysis is:

1. **Scope the question** with `/stisim:analysis-intake`, which interviews you to turn a research idea into a provisional analysis specification. It calls `analysis-selector` early on to check whether a dynamic transmission model is actually the right tool; if not, it will say so.
2. **Set up project memory** with `/stisim:project-memory`, so decisions survive across sessions, and use `/stisim:session-close` at the end of each session to write a handoff summary.
3. **Build the model** with `/stisim:model-writer`, drawing on `model-primer` (STIsim architecture), `hiv-interventions` (testing, ART, VMMC, PrEP), and `network-data` (sexual network inputs from DHS data).
4. **Plan the calibration** with `/stisim:calibration-strategy`, which helps decide what to calibrate and to which targets.

## Skills

### Analysis workflow

| Skill | Purpose |
|---|---|
| `analysis-intake` | Interview-style intake that turns a vague research idea into a provisional analysis specification. |
| `analysis-selector` | Classifies the research question by analytical objective, and only routes to HIVsim/STIsim if a dynamic transmission model is actually needed. |
| `project-memory` | Sets up a durable memory strategy so a multi-month project survives session ends, machine switches, and collaborator handoffs. |
| `session-close` | Writes a handoff summary at the end of a session. |

### Model authoring and calibration

| Skill | Purpose |
|---|---|
| `model-primer` | Architectural reference for STIsim: sim assembly, diseases, networks, interventions, connectors, analyzers, demographics, and calibration knobs. |
| `model-writer` | Composes a new sim from a scoped research question. |
| `hiv-interventions` | Design interview for HIV testing, ART, VMMC, and PrEP, with data-source references. |
| `network-data` | Sexual network calibration inputs from DHS data. |
| `calibration-strategy` | Decides what, why, and whether to calibrate: parameter classification, targets vs. knobs, and failure diagnosis. |

### Software quality

| Skill | Purpose |
|---|---|
| `extending-stisim` | Before subclassing or monkey-patching STIsim, classifies the change and steers fixes upstream rather than into private forks. |
| `editable-dep-hygiene` | Keeps edits to editable (`pip install -e`) dependencies committed and shared. |
| `comment-hygiene` | Keeps shared-library comments and docstrings about the code, not about a particular project or session. |
| `result-extraction` | Uses `ss.Result` methods (`annualize`, `resample`, `to_df`) rather than hand-rolled aggregation, which can mishandle flows vs. stocks. |
| `writing-tests` | Favors a small number of scientifically meaningful tests over exhaustive enumeration. |

## Companion plugins

`stisim.ai` is intentionally narrow. These plugins cover adjacent topics and are designed to be used alongside it:

- **[starsim_ai](https://github.com/starsimhub/starsim_ai)**: general Starsim and Sciris skills, disease-modeling skills, and engineering-quality review. Install from within Claude Code with `/plugin marketplace add https://github.com/starsimhub/starsim_ai`.
- **[calib](https://github.com/InstituteforDiseaseModeling/calib-plugin)**: generic calibration methods (algorithms, likelihoods, diagnostics); `calibration-strategy` decides *what* to calibrate, and `calib` handles *how*. Install with `/plugin marketplace add https://github.com/InstituteforDiseaseModeling/calib-plugin`.
- **[canonize](https://github.com/emiliasimmons/canonize)**: skills for durable decision capture in modeling projects, complementing `project-memory`.
- **[idm_standards](https://github.com/InstituteforDiseaseModeling/idm_standards)**: IDM software-quality, style, and documentation standards. Install with `/plugin marketplace add https://github.com/InstituteforDiseaseModeling/idm_standards`.

For the full plugin layout and contributor notes, see [`stisim/ai/README.md`](https://github.com/starsimhub/stisim/blob/main/stisim/ai/README.md).
