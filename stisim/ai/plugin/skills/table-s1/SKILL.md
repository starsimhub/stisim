---
name: table-s1
description: Use when the user needs a parameter table for an STIsim or HIVsim model — a manuscript supplement ("Table S1"), a reference table for a collaborator, or a comparison of a model's parameters against the STIsim defaults or another model. Assembles the table from the user-guide parameter tables, adds the values used by the user's sim(s), and flags which values differ from the defaults and still need a source.
metadata:
  version: "0.1"
  versiondate: "2026-10-01"
---

# Table S1: model parameter table

## When to use

- The user asks for a parameter table, "Table S1", supplementary parameter table, or a list of the parameters and sources used in a model.
- A collaborator asks for STIsim's parameters for comparison with another model.
- Before submitting a manuscript that uses STIsim, to document every parameter the analysis relies on.

## When NOT to use

- The user wants the meaning of a single parameter. Point them to the relevant user-guide page or the [parameter reference](https://docs.stisim.org/parameters.html).
- The user wants calibration diagnostics or posterior plots. Use `calibration-strategy` or the `calib:` plugin.

## Framing

The STIsim user-guide pages are the single source for parameter defaults and their sources: each module's page has tables with the columns `Parameter | Default | Description | Source`. `scripts/table_s1.py` parses those tables and adds one value column per sim. Values that differ from the class default are marked with `*`. The docs page [Parameters](https://docs.stisim.org/parameters.html) is the default-only version of the same table.

The format follows the "Table 1" convention of descriptive tables in epidemiology, as implemented in the R package [table1](https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html) by Benjamin Rich. Rows are grouped under labeled module and section headers, values carry their units, and the table has a caption and footnotes. table1 shows one column per stratum. Here there is one column per simulation instead (for example the defaults alongside a calibrated model, or one column per country).

A default's Source is only as good as what the developers recorded. Many defaults are "Not recorded", and a calibrated model overrides many defaults. **Every `*` value in the user's table needs a source that only the user can supply.** Never invent or guess citations.

## Instructions

1. **Identify the sim(s).** Ask which sim(s) the table should describe: the final calibrated sim, one per scenario or setting, or defaults only. Get a function or script that builds each sim; do not reconstruct it from memory. Record `sti.__version__`.

2. **Build the base table.** Use local docs when working in an STIsim checkout; otherwise fetch the pages for the installed version:

   ```python
   import sys; sys.path.insert(0, '<path to this skill>/scripts')
   import stisim as sti, table_s1
   rows = table_s1.parse(ref=f'v{sti.__version__}')   # or parse(docs='path/to/stisim/docs')
   for label, sim in sims.items():
       rows = table_s1.add_sim(rows, sim, label=label)
   rows = [r for r in rows if any(v != 'not in sim' for v in r['values'].values())]
   ```

   The command-line form (`python table_s1.py --format csv -o table_s1.csv`) gives the default-only table.

3. **Add what the docs do not cover.** Add these as extra rows, using the same columns:
   - **Interventions:** for each module in `sim.interventions`, the pars the analysis sets (coverage, efficacy, eligibility, timing). The docs do not list these because they are set per analysis.
   - **Calibrated parameters:** the prior range and the posterior summary (or best-fit value) for each calibrated parameter, with "Calibrated" as the source. Get these from the calibration outputs, not from the sim.
   - **Demographics:** the data source (for example, UN WPP for the country and year range) and any scaling (`total_pop`).

4. **Resolve sources for changed values.** List every `*` value and every added row, then ask the user for each source: data set, publication, calibration, or assumption. For each citation the user gives, check the DOI or PubMed ID before including it. If the user cannot give a source, write "Assumption" rather than leaving the cell blank.

5. **Flag weak defaults.** List the defaults that are still "Not recorded" and that the model is likely to be sensitive to: transmission probabilities, natural-history durations, and network structure. Note that reviewers may ask about them, and suggest a source search, a sensitivity analysis, or calibration. Do not fill them in yourself.

6. **Choose the scope.** Ask whether the user wants the full table or a compact one. A compact table drops unchanged "Not applicable" rows and modules the analysis does not depend on. In either case, keep every row with a `*`.

7. **Write the output.** Use `table_s1.to_markdown(rows, caption=...)` or `table_s1.to_csv(rows, path)`. For Word, convert the markdown with pandoc or the `docx` skill. The caption names the model, version, and setting. Keep the footnotes, which define the notation, the source vocabulary, and the `*` marker. If the journal uses numbered references, convert the linked sources to a numbered list at the end.

## Reporting

Give the user the table file, the list of `*` values still awaiting a source, and the "Not recorded" defaults flagged in step 5.
