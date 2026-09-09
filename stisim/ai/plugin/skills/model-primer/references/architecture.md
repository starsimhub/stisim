# STIsim architecture — mental model

The purpose of this file is to give the reader a mental model of how STIsim's pieces compose, in the reader's own working memory, so that they can navigate the codebase and reason about the model without holding the whole thing in their head at once.

*Scaffold; content to be written by the lead developer.*

## Orientation

*One-paragraph "what is STIsim" — how it relates to Starsim, what it's for, what it isn't.*

## The Sim as the assembly point

*How `sti.Sim` composes diseases, networks, connectors, interventions, analyzers, and demographics. Why the assembly is a list of modules, not a dict. What happens on `init` vs on `step`.*

## Diseases

*The base `STI` class and the disease hierarchy. What all STIs share (natural history, symptomatic behaviour, care-seeking hooks). What HIV and syphilis do that the bacterial STIs don't. How BV is different from the sexually-transmitted diseases in the same tree.*

## Networks

*The sexual-network model: MF, MSM, FSW, and how they layer. What edges represent. Duration, concurrency, and how partners are drawn.*

## Interventions

*The intervention abstraction. When interventions run. How they compose with diseases (screening → treatment → care-seeking).*

## Connectors

*What connectors are for. Concrete examples: HIV-syph coupling, pregnancy-STI.*

## Analyzers

*How analyzers observe the sim without changing it. When they run.*

## Demographics

*Background population, births, deaths, migration. How they interact with disease dynamics.*

## Care-seeking

*The care-seeking model as a first-class module rather than a per-disease detail. Why.*

## Timestep and time discipline

*Monthly timestep. What that means for pregnancy windows, waning immunity, event scheduling. Traps for the unwary.*

## Common idioms and gotchas

*A short list — not exhaustive — of the things that trip people up. Populated from institutional memory as they arise.*
