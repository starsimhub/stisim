---
name: editable-dep-hygiene
description: Use when about to edit a file inside an editable pip install (stisim, starsim, or any other `pip install -e` dependency), or about to bump the version or switch the branch of such a dependency. Enforces the rule that any edit to an editable dep must be committed, pushed, and PR'd immediately — and that no version bump or branch switch happens with unmerged local work in the dep checkout.
---

# Editable-dep hygiene

## When to use

- About to edit any file inside an editable pip install (stisim, starsim, fpsim, hpvsim, etc. — anything installed with `pip install -e`).
- About to bump the version of an editable dep, switch its branch, or pull upstream.
- Reviewing a work-in-progress state where an editable dep's working tree or branch state is unclear.
- Trigger phrases: "let me patch stisim locally", "edit the starsim source", "upgrade stisim", "switch stisim to rc…", "the editable checkout of X", "pull the latest stisim", "bump the dep version".

## When NOT to use

- The user is editing a file inside the project workspace itself (not inside an editable dep) — that is normal project work, no dep hygiene concern.
- The user is editing a released, non-editable dep (which they shouldn't be — that is a signal to escalate before continuing).

## Framing

An editable pip install (`pip install -e ../stisim`) points python at a git checkout. Any edit lives in that checkout's working tree — outside the project repo's history and outside the dep's shared history until it is committed and pushed. Two consequences follow:

1. **Local edits are invisible to anyone else** — including future you on a different machine, and any experiment that runs after a version bump replaces the file.
2. **Version bumps and branch switches routinely discard uncommitted work.** `git reset --hard <tag>`, `git checkout <branch>` after `git stash`, `pip install --force-reinstall`, or an environment recreate all silently remove uncommitted changes. Even committed changes on a local-only branch are lost from view when the checkout switches to a different branch.

The canonical failure: the exp-005 VMMC fix. A local edit to the editable stisim checkout worked for months. A later stisim upgrade replaced the file. The bug returned. It took months to trace back to the version bump because the fix had never been in stisim's git history to begin with — nothing was "wiped", it just was never there to survive. Full account in `experiments/015_vmmc_prevalence_target/README.md` of the hivsim_eswatini repo.

The rule this skill enforces: **any edit to an editable dep gets committed and PR'd immediately, and no version bump happens with unmerged work in the dep checkout.**

## Instructions

1. **Before editing a file inside an editable dep, name that you are doing so.** State the dep, the file, and — from consulting `extending-stisim` — whether this belongs upstream at all. If it belongs upstream, the edit is a temporary hotfix on the road to a PR, not a persistent local change.

2. **Immediately after the edit compiles**, commit to a branch in the dep's checkout with a message referencing the downstream context, and open a PR against the dep. Not "when I'm done for the day", not "after I test it in the sim" — immediately, so the change survives any subsequent version bump or branch switch. The PR can be a draft; it exists to make the change visible and versioned.

3. **Record the dep-side PR URL in the downstream exp / SUMMARY / commit** that depends on the fix. This makes the coupling explicit: "this experiment's result depends on unmerged dep PR #NNN". Without this, a fresh checkout that has the merged upstream version will silently behave differently.

4. **Before bumping the dep version, switching its branch, or pulling upstream**, verify in the dep checkout:
   ```
   cd $(python -c "import <pkg>; import os; print(os.path.dirname(os.path.dirname(<pkg>.__file__)))")
   git status                        # no uncommitted changes
   git log ..origin/<current-branch> # no unpushed local commits
   git branch --show-current         # on the branch you think you are on
   ```
   If any of these show local work, resolve it first — either finish the PR flow (commit + push + PR), or explicitly acknowledge that the work is being discarded. Never bump through unmerged work.

5. **After a version bump**, sweep for downstream artifacts that were workarounds for now-fixed upstream bugs (see `extending-stisim` step 5). A version bump is the natural moment to delete a downstream subclass whose upstream fix has landed.

## Checks before completion

1. No file inside an editable dep has been edited without a corresponding commit + PR in the dep's repo.
2. Any downstream exp / commit that depends on an unmerged dep PR records the PR URL.
3. Any version bump / branch switch of an editable dep was preceded by a clean-state check on the dep checkout, and any local work resolved (merged, PR'd, or explicitly discarded).
