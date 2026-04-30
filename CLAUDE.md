The refactored devtest file shows clarifications to the design philosophy and API sketches — no new gotchas or implementation patterns discovered. The file is now more readable with better inline comments explaining:

- Why `on_prep` lives on the HIV object (coordination without a manager)
- The gap-proportional allocation logic for supply constraints
- Open questions on Clark's PR #432 approach (manager dependency, CRN breakage, disabled reuptake)

This is pure documentation/clarity work. No new patterns or blockers emerged.

SKIP
