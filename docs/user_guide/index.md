# User guide

STIsim is built on top of [Starsim](https://docs.starsim.org). For understanding the core model architecture -- agents, modules, parameters, time handling, results, and the simulation loop -- the [Starsim user guide](https://docs.starsim.org/user_guide/intro_starsim.html) is your primary reference.

This user guide documents what STIsim adds on top of Starsim:

- **[Diseases](diseases/index.md)** -- State diagrams and parameter tables for each STI module (HIV, syphilis, chlamydia, gonorrhea, trichomoniasis, BV, GUD).
- **[Networks](networks/structured_sexual.md)** -- The structured sexual network: risk groups, partnership types, age mixing, sex work, and condom use. See also [pair matching](networks/matching.md) and [MSM networks](networks/msm.md).
- **[Interventions](interventions/interventions.md)** -- Testing, treatment, HIV-specific interventions (ART, VMMC, PrEP), syndromic management, ANC/infant (PMTCT) screening, partner notification, and pregnancy-driven risk reduction.
- **[Analyzers](analyzers.md)** -- Observe and record derived results: coinfection, sex-work transmission, and network/partnership structure.
- **[Care-seeking](care_seeking.md)** -- Cross-disease care-seeking propensity shared across interventions.
- **[Connectors](connectors.md)** -- Coinfection interactions between co-circulating diseases.
- **[Calibration](calibration.md)** -- Fitting model parameters to data with Optuna.
- **[HIVsim](hivsim.md)** -- The HIV-focused convenience wrapper around `sti.Sim`.
