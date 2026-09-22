---
name: extending-stisim
description: Use when about to subclass a stisim / starsim class, override a library method, or write a downstream workaround for behavior that seems wrong upstream. Enforces the upstream-vs-downstream decision — real bugs go upstream as PRs, opt-in project knobs are fine downstream with a no-op default, and downstream workarounds get deleted the moment their upstream fix lands.
---

# Extending stisim

## When to use

- About to write a subclass of `sti.*` or `ss.*` to change library behavior.
- About to override a library method, monkey-patch, or replace a library-provided intervention with a local variant.
- Reviewing an inherited repo for downstream subclasses that may be dead code (upstream fix has landed but the subclass was never deleted).
- Trigger phrases: "subclass sti.HIV", "override this stisim behavior", "the upstream VMMC/ART/… is wrong", "let me patch stisim locally", "workaround for a stisim bug", "in-repo intervention subclass".

## When NOT to use

- The user is writing an `ss.Analyzer` — analyzers are always downstream, that is what they are for.
- The user is writing country-specific data preprocessing (mortality reconstruction, coverage cleaning) — that is a data-construction script, belongs downstream by definition. See `mortality_construction.py`-shape scripts for the canonical pattern.
- The user is adding a parameter that already exists on the parent — consult `model-primer/references/calibration-knobs.md` first for what stisim already exposes.

## Framing

**The goal is to share fixes back to the community that benefits from them.** When you find a bug or gap in stisim, the default action is an upstream PR — every stisim user benefits, and future you gets code review from the maintainers rather than maintaining an untested private fork alone. If you spot a bug you cannot fix in this session, file an issue upstream so at least the knowledge is shared; silent knowledge of a bug is also siloing.

**The cost of a downstream workaround is not that a version bump might wipe it — it is that it stays siloed.** Every other stisim user hits the same bug until the fix lands upstream. Every downstream subclass is code that only your team can maintain, review, and improve. AI makes downstream subclassing especially tempting (fast to draft, feels productive) and especially dangerous (unshared code accumulates faster than the org can review it, and code review — not code production — is the binding constraint on org throughput now). The version-bump-wipes-the-fix problem is a symptom of the deeper anti-pattern: a private fork of shared code.

Two questions decide where a change lives. Ask them in this order:

1. **Am I fixing broken behavior, or adding a project-specific opt-in knob?** If the change corrects a bug or gap that every stisim user would want fixed, it belongs upstream. If it is a research-specific handle that other users would not want as the default, it can live downstream — but only as an opt-in with a no-op default.

2. **If it is a fix, when does the upstream PR happen — before, after, or with the local hotfix?** With. A downstream workaround with no matching PR is a fix that will silently vanish at the next version bump (see `editable-dep-hygiene` for why). The workaround is a hotfix, not a fix, until the PR is merged.

The `model-primer/references/calibration-knobs.md` "Framing" section already states the three modes (**feed data / calibrate / upstream fix**) that a candidate change falls into. Consult it before writing anything.

Common failure modes this skill prevents:

- **Reinvented existing knob.** `HIVMortalityMultiplier` was written as a `sti.HIV` subclass exposing a CD4-death-rate multiplier; `pars.rel_death` already existed for exactly this. Sixty-line subclass, one-line kwarg.
- **Method override when a callable parameter suffices.** `AgeDependentSurvival` overrode `set_prognoses` to make `dur_latent` depend on age at infection; starsim distributions accept callable parameters ([`starsim/distributions.py:939`](../../../../../../starsim/starsim/distributions.py)) — `dur_latent = ss.lognorm_ex(mean=age_func, …)` achieves the same result in ~10 lines.
- **Real fix kept downstream forever.** `VMMCPrevalenceTarget` was a legitimate workaround for a real stisim VMMC bug; the exp README rationalised keeping it in-repo "so a stisim git pull cannot silently wipe it again". The correct response was the PR that eventually landed as #535 — months later, with every other stisim user hitting the same bug in the interim.
- **Never cleaned up after upstream landed.** `VLSStockTarget` was upstreamed in commit `45a8b21`; the downstream subclass is still the default in the exp repo. Dead code that will confuse the next reader and drift out of sync with the upstream API.

## Instructions

1. **Classify the change.** State plainly whether this is (a) a fix to broken or wrong behavior in stisim / starsim, (b) an opt-in research knob, or (c) a project-specific data or workflow concern that does not belong in the library at all. Do not proceed until this is clear.

2. **Before writing any subclass, grep the parent for an existing par or an equivalent mechanism.** Read `model-primer/references/calibration-knobs.md` for the canonical index. Check `define_pars(…)`, `self.<name> = …` attributes in `__init__`, and — for anything that draws from a distribution based on per-agent state — starsim's callable-parameter mechanism. If the mechanism already exists, use it and do not subclass.

3. **If it is a fix (a)** — write the change as a PR against the dep, not as a downstream subclass. A downstream hotfix is acceptable *only* if the PR is opened in the same session and the exp / SUMMARY records the PR URL. If the fix is beyond what you can PR in this session (needs a maintainer conversation, a design decision, or expertise you don't have), **file an issue in the dep repo describing the bug and the workaround you used**, and link the issue from the exp — the knowledge is then shared even if the fix is not yet.

4. **If it is a research knob (b)** — ensure the subclass has a documented no-op default (a same-seed sim with the knob at its no-op value must be bit-identical to the base class). Opt-in via `hiv_class=` or an explicit `interventions=[…]` entry — never make it the default in a shared config.

5. **When an upstream PR that supersedes a downstream workaround merges**, delete the downstream artifact in the same commit that bumps the dep version. Grep the repo for uses of the old subclass and switch them to the upstream API. This is the point most often skipped; it is what left `VLSStockTarget` alive in the exp repo months after the upstream fix.

6. **On any inherited repo**, sweep for downstream `sti.*` / `ss.*` subclasses and ask, for each: does this reinvent an existing par? Is it still needed given current upstream? Is the "kept in-repo so git pull cannot wipe it" rationale in its docstring — because that framing itself is the anti-pattern.

## Checks before completion

1. If the change is a fix: an upstream PR exists (URL recorded in the exp / SUMMARY / commit message), or the exp explicitly names the downstream-only status as technical debt with a deadline.
2. If the change is a research knob: it defaults to an exact no-op, is opt-in via constructor arg, and does not shadow an existing par.
3. If an upstream PR that supersedes a local workaround has merged, the local workaround has been deleted and its call sites switched to the upstream API — in the same commit that bumps the dep version.
