"""
PrEP API vignettes — design sandboxing for issue #431 / v1.6 milestone
=======================================================================

Purpose
-------
These vignettes define what the *user-facing* PrEP API should look like before
we commit to any implementation. The goal is to write the code we wish existed,
then work backwards to figure out what needs to change in stisim. Read this file
top-to-bottom: it starts from what runs today and escalates in complexity toward
the Kenya CE research question that is motivating the whole design effort.

The vignettes also serve as a reference when evaluating Clark's draft PR #432,
which proposes one implementation path. That approach is shown in V3b for direct
comparison.

Run status
----------
  V1, V2a   — run today against rc1.5.4 with no changes to stisim
  V2b–V5    — aspirational API; commented out or pass-bodied (marked NOTIMPLEMENTED)


Guiding principles
------------------
1. "Common things should be simple."  The one-liner case must be a one-liner.
   Complexity should be opt-in, not the default.

2. Eligibility is always explicit.  The current Prep class targets FSW by default,
   which is not self-explanatory and will surprise anyone not already familiar with
   the codebase. The new API requires the user to say who gets PrEP. If truly
   universal coverage is desired, the eligibility lambda makes that visible:
       eligibility=lambda sim: ~sim.diseases.hiv.infected

Clark response: My preference is no defaults and users must provide eligibility and coverage data to use. However,
I used the current default-based setup because I was asked to.

3. Self-documenting code.  A vignette should read like a description of the
   intervention, not a configuration puzzle. Labels, group names, and product
   names should appear in results so figures are interpretable without the source.

4. Results first.  Every new feature must produce something accessible in
   sim.results so it can be plotted, exported, and fed into a CE analysis without
   post-hoc arithmetic on internal state.

Clark response: Sure, as noted below, I would be happy to discuss regarding what default outputs would be available
and where..


Research context — why this matters now
----------------------------------------
Our collaborator at CEMA (Kenya) has asked:
  "Estimate the price threshold at which Kenya would be willing to pay for
   lenacapavir. We'd also like to use HIVsim to evaluate the current rollout."

Lenacapavir (brand name Sunlenca) is a 6-monthly injectable with ~99.9% efficacy
in the PURPOSE 1 trial. At current Gilead pricing (~$40,000/year) it is far above
any plausible WTP threshold for Kenya. The question is where the threshold lies
and how the answer changes with targeting strategy.

Target groups under active consideration: FSW, AGYW (women aged 15–24), PBFW
(pregnant and breastfeeding women), MSM. These map onto different sub-networks
and demographic filters in STIsim; see eligibility lambdas in V5b.


Current state of sti.Prep (rc1.5.4)
-------------------------------------
File: stisim/interventions/hiv_interventions.py:435

What it does:
  - Per-timestep Bernoulli coverage model (coverage = probability of being on PrEP
    at each step, not a stock).
  - Applies rel_sus *= (1 - eff_prep) to newly covered agents each step. There is
    no clean reversal on dropout; agents who fall off PrEP retain the reduced
    susceptibility until the next step.
  - on_prep state lives on the intervention object as a BoolArr.
  - Eligibility is overridable via a lambda, but defaults to FSW (see above —
    this default should be removed).

What it lacks:
  - No concept of a product (oral vs. injectable vs. lenacapavir).
  - No cost tracking, so no CE analysis is possible.
  - No duration model (relevant for injectables where re-dosing has a schedule).
  - No multi-group API; targeting multiple groups requires separate instances,
    which cannot share a supply budget.
  - on_prep on the intervention means two concurrent Prep instances have no way
    to coordinate — agents could be double-enrolled and double-discounted.


Clark's PR #432 — overview and open questions
----------------------------------------------
Branch: 330 (draft, not yet merged)
PR: https://github.com/starsimhub/stisim/pull/432

Clark's PR introduces a supply-chain abstraction:
  Product   — name, type, delivery_mode, cost, efficacy curve (eff_by_ti[])
  Supply    — quantity of a given Product; tracks quantity remaining + accrued cost
  Supplies  — container for one or more Supply objects
  SuppliedPrep — new Prep class using the above
  PrepManager  — a separate intervention that manages dropout timing

The Product/Supply/Supplies utility classes are well-designed and well-tested.
The concerns are in the integration layer:

  CONCERN 1 — PrEP state on the disease object.
    PR #432 adds 8 new states to hiv.py: on_prep, prep_source, ti_prep_start,
    ti_prep_end, ti_prep_drop, prep_n_reuptake, prep_product, prep_eff. By Starsim
    convention, state that belongs to an intervention should live on the
    intervention, not on the disease. These states encode the PrEP program
    logic, not the HIV natural history. The exception — see below — is a single
    shared coordination flag (on_prep).

Clark response: The reason for the new states, which are only added if a PrEP intervention has been added to the
simulation, is that they are required for inter-PrEP communication and coordination. These are "global" values, not
specific to any one PrEP intervention, because there can be > 1 PrEP-related intervention in a simulation and they
can be active at the same time.

  CONCERN 2 — PrepManager as a required singleton.
    PrepManager must be added to sim.interventions alongside every SuppliedPrep;
    if omitted, agents never leave PrEP (the dropout step is never called). This
    is a hidden dependency that will catch users out. The ordering dependency
    (PrepManager.start_step must fire before SuppliedPrep.step) is also fragile.

Clark response: False, PrepManager is not added "alongside every SuppliedPrep". It is auto-added at most once for a
simulation, if at least one PrEP intervention has been detected. A user does not need to know about its existence
unless they wish to subclass it for different behavior. This "not need to know" functionality will need the
prep uptake() staticmethod is moved elsewhere. Perhaps an input to the SuppliedIntervention __init__, but this is
open for discussion. Also, "The ordering dependency ... is also fragile" is False. The class is ONLY auto-added, and
always added at the *front* of the intervention list, which *guarantees* the PrepManager.start_step will run before all
functions of all prep interventions. The primary, and potentially only, purpose of PrepManager is to force agents off
of prep at the end of their current enrollments, to allow them to subsequently, with no temporal delay, choose to
enroll in any prep intervention, regardless of intervention order, if they are eligible and care seeking. Unless a user
specifically sabotages their own auto-generated sim object, it is impossible for ordering to be wrong. At that point,
it is their active fault, not the fault of the model/intervention.

  CONCERN 3 — random.sample breaks CRN.
    SuppliedPrep.step() uses Python's random.sample for agent selection. This
    bypasses Starsim's Common Random Numbers system. All agent selection should
    go through ss.bernoulli or ss.choice.

Clark response: Easy enough to fix if there is a direct slot-in using ss functionality. For a draft PR, I simply knew
about random.sample, but if ss.choice does the same thing, happy to fix.


  CONCERN 4 — Reuptake is fully disabled.
    The re-enrollment logic (~70 lines) is commented out with a "TODO: pending
    clarification" note. This means agents who complete a lenacapavir injection
    cycle never re-enrol, which is epidemiologically wrong for any realistic
    PrEP program.

Clark response: Reuptake can easily be re-enabled. It was disabled due to what I interpreted as uncertainty whether
it was needed or provided value given that we may or may not have data to set it properly. It is currently disabled
to reduce complexity until if/when we need it. The assertion that agents "completing a lenavapavir injection cycle"
cannot re-enroll is False. They simply have to: be eligible and care seeking to reuptake an intervention at the end
of their previous enrollment. Any agent can switch PrEP interventions or reenroll in the same one with no temporal
gaps in the current design.

  CONCERN 5 — Hardcoded care-seeking rate.
    base_care_seeking_rate = 0.8 is labelled "TODO, arbitrary" in the source.
    The care-seeking mechanism itself (agents have a scalar care-seeking draw
    that modulates uptake probability) may be worth keeping, but the value needs
    a defensible basis or a parameter.

Clark response: Happy to adjust, 0.8 is, as noted, "completely arbitrary" (but necessary for creating functioning code)
and awaiting more information from Adam. The full care seeking calculation is awaiting details of how it should be
calculated, hence a draft PR and notation in my presentations as needing research input to finalize.

These concerns are fixable; the question is whether it's easier to fix them
within Clark's architecture or to redesign with his utility classes as a foundation.

Clark response, Here is what needs to be done to address any of the above:
- Care seeking does need research input to define it exactly, which is already known and has previously been communicated.
- Slight refactoring to move (non-critical) functionality currently in PrepManager is probably necessary to enable
users to 1) not worry about instantiating or ordering PrepManager properly, and 2) enable them to subclass it to
provide alternate drop() functionality, if needed. But this is a simple update.
- random sampling should be simple to update, provided existing ss functionality to probabilistically select N items from
a list of M items (N <= M, each with a unique weight for selection). If not, then new ss functionality to do the
equivalent is necessary, but happy to have it use ss code here.
- new prep states on hiv module are required, due to the potential, and possible desirability, of > 1 prep intervention
that can be active at any given time. Minimum requirement (I think): on_prep, prep_source. The others can be discussed
case-by-case, but we should always keep in mind that there can be > 1 prep intervervention and if interventionA might
need to know "global" prep state (information from other interventions).


Proposed design for the new Prep class
----------------------------------------
The sketches below (V2b–V5) assume the following architecture. These are proposals,
not decisions — the point of this file is to generate discussion.

  A. on_prep on HIV (one state only).
     A single hiv.on_prep BoolArr shared across all Prep interventions solves
     the coordination problem without a PrepManager. On uptake, each Prep sets
     hiv.on_prep[uids] = True and applies rel_sus. On dropout, it clears the flag
     and restores rel_sus. Each Prep instance additionally owns a per-program
     BoolArr (e.g. self.enrolled) to track who is in *its* program specifically.
     Only hiv.on_prep needs to be on the disease; everything else stays on the
     intervention.

  B. groups= list of dicts on Prep.
     Rather than instantiating one Prep per group, the user passes a list of
     group dicts. Each has: label (str), eligibility (callable), coverage (float
     or time-series). The Prep class tracks per-group actual coverage and, when
     supply= is set, allocates doses proportionally to the coverage gap. This is
     essentially Clark's gap-proportional logic, but contained within one class.

  C. product= dict (or dataclass) on Prep.
     Bundles: name, efficacy, cost_per_year (or cost_per_dose), and optionally
     duration. If duration is set, the model tracks when each agent's dose
     expires (ti_prep_end on the intervention, not on HIV) and handles re-enrolment.
     If duration is omitted, behaviour is as today: continuous coverage, on/off
     each step. Cost is accumulated in init_results/step and exposed via
     sim.results.<name>.total_cost.

  D. supply= int/float on Prep.
     Total doses available per year (or per timestep, scaled). When set, total
     uptake across all groups is capped. Allocation across groups is proportional
     to gap_i = (target_coverage_i - actual_coverage_i) * group_size_i.
     If supply= is not set (the common case), there is no cap.

  E. No PrepManager.
     Dropout and re-enrolment logic lives in Prep.step(), ordered correctly
     within that single method. The shared hiv.on_prep flag means multiple Prep
     instances can coexist safely.


Open questions
--------------
These are the decisions that need to be made before implementation begins.
Comments / votes welcome — add your initials next to any position you support.

  Q1. Where does on_prep live?
      Option A: hiv.on_prep — one shared BoolArr on the disease object (proposed above).
      Option B: intervention-side only — each Prep owns its own BoolArr; multi-
                intervention coordination via lazy union across instances. Cleaner
                architecture but requires iterating sim.interventions each step.
      Option C: Clark's approach — 8 states on hiv.py, fully centralised.

Clark response: See above for the need for centralizing at least some of this data. Some of the data is potentially
up for discussion, again, see above.

  Q2. Should the product= interface accept a plain dict, or a typed dataclass/class?
      A plain dict is the Starsim way (cf. how pars are passed). A class (like
      Clark's Product) provides tab-completion and validation, but adds a new
      public type for users to learn. Possible middle ground: Product is an
      internal implementation detail that Prep constructs from a dict internally.

Clark response: I would suggest that it is not necessarily the "starsim way" to only pass dicts of things. For example,
HIVPars objects can be passed as parameter sets instead of kwargs (effective dict) during instantiation. This creates
a "new interface" for a user to use if they wish. The current pattern of Supplies follows this and allows for an
extensible object (a Supplies object) that can contain additional future behavior (or be subclassed for specific
additional features). A simple dict of Supply objects does not allow for this naturally.

  Q3. Is gap-proportional allocation the right supply-constrained default?
      Clark's implementation allocates remaining doses proportional to the coverage
      gap in each group. This is reasonable but privileges groups with higher
      targets. Alternatives: proportional to group size, priority ordering, equal
      split. Should the allocation rule be user-configurable?

Clark response:
It is reasonable to ask if this is a reasonable implementation. One important item to note here is that in the limit
where available supply goes from 0 to "no restriction this step", all group targets are *simultaneously* hit when there
are exactly enough doses to cover enough agents to hit the targets.

"privileges groups with higher targets": False
The current proposed implementation does not privilege "groups with higher targets". It privileges
groups with a larger gap from their current target. Example, assuming supply limitation is hit this step:
groupA: target: 80%, current: 60%, gap: 80 - 60 = 20%
groupB: target: 50%, current: 10%, gap: 50 - 10 = 40%

Proportionality of limited supply distribution:
groupA: 20 / (20 + 40) = 33.33%
groupB: 40 / (20 + 40) = 66.67%

In this case, as supplies are (conceptually) increased from "not enough" to "enough" for this step, groupA and groupB
will hit their targets (80% and 50%, respectively, at exactly the same moment (when supplies are "just enough"). This
means that the current algorithm is prioritizing the groups exactly as needed to keep supply shortages from impacting any
one group in a way to disadvantage it "unfairly" from reaching its target coverage.

Note, the "lower target group" gets prioritized (proportionally by %) in this case, because it farther from its coverage
target in coverage/% units.


  Q4. How should dose duration / re-enrolment be modelled?
      For lenacapavir (6-month injection), an agent who starts a dose at time t
      is protected until t+6 months. At expiry they can either: (a) automatically
      re-enrol (if still eligible and supply allows), (b) drop out permanently,
      or (c) drop out with some probability and re-enrol otherwise. Do we need a
      per-agent re-enrolment probability, or is a population-level dropout rate
      sufficient?

Clark response: Again, as noted above, I disabled auto-reenroll due to it being unclear that such additional more
complicated logic was necessary or desired. It would be easy enough to restore auto-reenroll if desired.

  Q5. How do we handle PrEP + ART interaction?
      An agent who seroconverts while on PrEP should come off PrEP and (ideally)
      be prioritised for ART. The current Prep.step() doesn't explicitly handle
      this; agents who seroconvert remain in the BoolArr until the next step. Is
      that acceptable, or do we need an explicit seroconversion handler?

Clark response: It is a feature of the current implementation that agents stay "on_prep" until PrEP durability end after
seroconverting. This explicitly allows for detecting agents who are on a type of PrEP AND HIV+ (after receiving
PrEP),  which can be a source for drug resistance. Agents should never change their behavior until they have been
*diagnosed*, and even then, a lenacapavir-injected seroconverter cannot have the drug removed from their system until
it "wears off" naturally. Additionally, there will be cases where an agent is *diagnosed* while still on_prep (especially
injectable), and a healthcare provider will need to determine if they can be given ART (now, or later). The current
design allows for these cases to be captured.

  Q6. What results should Prep expose?
      Minimum for CE analysis: total_cost (cumulative), n_on_prep (timeseries),
      n_new_initiations (timeseries). Per-group breakdowns are needed for V5b.
      Should these be in sim.results.<prep_name> or sim.results.hiv?

Clark response: Happy to record all kinds of things like this on whatever the appropriate result object is

  Q7. AGYW and PBFW are not pre-defined convenience properties.
      Currently, sim.networks.structuredsexual exposes fsw and client as boolean
      properties. AGYW (female, age 15–24) and PBFW (pregnant or breastfeeding)
      require lambda expressions (see vignettes). Should we add convenience
      properties, and if so, where? (network? demographics? a new targeting module?)

Clark response: Um, this is prep intervention feedback? I guess I'm missing the point.

  Q8. Should care-seeking be part of the PrEP uptake model?
      Clark's SuppliedPrep modulates uptake probability by hiv.care_seeking, a
      per-agent continuous trait drawn from N(1, 0.5). This captures heterogeneity
      in health-seeking behaviour and is epidemiologically defensible. It adds
      complexity. Is it in scope for v1.6, or a v2.0 feature?

Clark response: If no care seeking calculation is performed, then *all* agents eligible will *immediately* jump on
the intervention in one giant step-jump of on_prep. Also, I believe we *do* want to capture heterogenerity in whether
agents do something, rather than just bulk-putting-agents-on-prep. That approach may work for some questions and use
cases, but very likely not all. Moreover, the default care seeking rate and calculation does need scientific input
from Adam. We could very easily set the default to be "all are seeking" as well.

"""

import sciris as sc
import starsim as ss
import stisim as sti


# =============================================================================
# V1 — Simple single-group PrEP
# Runs today against rc1.5.4 with explicit eligibility.
# =============================================================================

def v1_fsw():
    """
    Simplest useful case: one group, constant coverage, default oral efficacy.

    Compared with the current rc1.5.4 Prep, the only change is that eligibility
    is explicit. The FSW default has been removed — see design note above.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep = sti.Prep(
        coverage=0.5,
        eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    )
    sim = ss.Sim(n_agents=500, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    sim.run()
    return sim


def v1b_agyw():
    """
    Same structure; different group and time-varying coverage ramp.

    AGYW = adolescent girls and young women (here: female, age 15–24).
    No convenience property exists yet — see open question Q7.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep = sti.Prep(
        coverage={'year': [2015, 2020], 'value': [0, 0.4]},
        eligibility=lambda sim: (
            sim.people.female & ~sim.diseases.hiv.infected
            & (sim.people.age >= 15) & (sim.people.age < 25)
        ),
    )
    sim = ss.Sim(n_agents=500, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    sim.run()
    return sim


# =============================================================================
# V2 — Multi-group targeting
#
# Two sub-vignettes showing the same scenario two ways:
#
#   V2a: one Prep instance per group (works today).
#        Limitation: each instance manages coverage independently. If supply is
#        shared across groups, this approach cannot allocate it correctly.
#
#   V2b: single Prep with groups= list (NOTIMPLEMENTED).
#        One instance sees all groups and can enforce shared supply awareness.
#        This is the preferred API; V2a is the interim workaround.
# =============================================================================

def v2a_multigroup_instances():
    """
    One Prep per group. Runs today; name= is required to avoid key collision.

    Adequate when supply is unlimited. Breaks down when a dose budget needs to
    be allocated proportionally across groups (see V4).
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep_fsw = sti.Prep(
        name='prep_fsw',
        coverage=0.6,
        eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    )
    prep_agyw = sti.Prep(
        name='prep_agyw',
        coverage=0.3,
        eligibility=lambda sim: (
            sim.people.female & ~sim.diseases.hiv.infected
            & (sim.people.age >= 15) & (sim.people.age < 25)
        ),
    )
    sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net,
                 interventions=[prep_fsw, prep_agyw], dt=1/12)
    sim.run()
    return sim


def v2b_multigroup_groups_api():
    """
    NOTIMPLEMENTED — preferred multi-group API.

    A single Prep instance manages all groups. Per-group coverage is tracked
    separately; when supply= is set, doses are allocated proportional to the
    per-group coverage gap (see V4). Results exposed as sim.results.prep.fsw.*
    and sim.results.prep.agyw.* (exact structure TBD — see Q6).

    PBFW eligibility uses sim.demographics.pregnancy.pregnant — this requires
    the Pregnancy module to be in the sim. If it's absent, the lambda would
    raise AttributeError; some defensive handling or documentation is needed.
    """
    # NOTIMPLEMENTED: groups= kwarg on Prep
    #
    # prep = sti.Prep(
    #     groups=[
    #         dict(label='fsw',
    #              eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    #              coverage=0.6),
    #         dict(label='agyw',
    #              eligibility=lambda sim: sim.people.female & ~sim.diseases.hiv.infected
    #                                     & (sim.people.age >= 15) & (sim.people.age < 25),
    #              coverage=0.3),
    #         dict(label='pbfw',
    #              eligibility=lambda sim: sim.demographics.pregnancy.pregnant & ~sim.diseases.hiv.infected,
    #              coverage=0.5),
    #     ],
    # )
    # sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    # sim.run()
    pass


# =============================================================================
# V3 — Product-aware PrEP (NOTIMPLEMENTED)
#
# "Product" bundles together the clinical and economic properties of a specific
# PrEP formulation: efficacy, cost per year (or per dose), and optionally a
# dose duration for injectables.
#
# Products of interest:
#   oral TDF/FTC   efficacy ~80%,  cost ~$50/year,      continuous (daily pill)
#   CAB-LA         efficacy ~90%,  cost ~$200/year,     duration ~2 months
#   lenacapavir    efficacy ~99.9%, cost ~$40,000/year, duration ~6 months
#
# Duration matters for:
#   (a) re-enrolment logistics (see Q4): once injected, an agent cannot "drop off"
#       mid-cycle the way they can with a daily pill.
#   (b) cost attribution: the cost of one lenacapavir injection covers 6 months,
#       so cost-per-timestep accounting needs to spread it correctly.
#
# Two sub-vignettes are shown:
#   V3a — proposed minimal API (product= as a dict on Prep)
#   V3b — Clark's Product/Supply/Supplies class API (PR #432), shown for comparison
# =============================================================================

def v3a_product_dict_api():
    """
    NOTIMPLEMENTED — product= dict on Prep.

    The dict is the Starsim-conventional way to pass structured parameters.
    Prep constructs a Product object internally if needed; the user never imports
    a Product class. Efficacy from product= overrides the eff_prep scalar.

    Cost is accumulated per timestep and written to sim.results.<name>.total_cost.
    Duration, if provided, activates the re-enrolment model (see Q4); if omitted,
    coverage is continuous (oral behaviour).

    Coordination note: with two concurrent Prep instances, hiv.on_prep (a single
    shared BoolArr on the disease — see design proposal C above) prevents double-
    enrolment without requiring a PrepManager. Each instance checks ~hiv.on_prep
    before offering uptake; on uptake it sets the flag and adjusts rel_sus; on
    dropout it clears the flag and restores rel_sus.
    """
    # NOTIMPLEMENTED: product= on Prep
    #
    # oral = sti.Prep(
    #     name='oral',
    #     coverage=0.5,
    #     eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    #     product=dict(name='oral_tfv_ftc', efficacy=0.80, cost_per_year=50),
    # )
    # lnc = sti.Prep(
    #     name='lnc',
    #     coverage=0.5,
    #     eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    #     product=dict(name='lenacapavir', efficacy=0.999, cost_per_year=40_000,
    #                  duration=ss.dur(6, 'months')),
    # )
    # sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net, interventions=[oral, lnc], dt=1/12)
    # sim.run()
    # print('oral total cost: ', sim.results.oral.total_cost[-1])
    # print('lnc  total cost: ', sim.results.lnc.total_cost[-1])
    pass


def v3b_product_clark_api():
    """
    NOTIMPLEMENTED — Clark's Product/Supply/Supplies API from PR #432.

    Shown here for direct comparison with V3a. These modules exist on branch
    '330' only; they are not in rc1.5.4.

    What Clark's approach does well:
      - Product, Supply, Supplies are clean, well-tested utility classes.
      - Supply.use() enforces quantity bounds and accrues cost atomically.
      - Multiple Supply objects (e.g. oral + injectable) can coexist in one
        Supplies container and be queried by name or type.

    Open concerns (see full list in module docstring):
      1. PrepManager is a required but implicit dependency. If the user forgets
         to add it, agents are never dropped from PrEP — silent wrong behaviour.
      2. 8 PrEP states added to hiv.py couple the disease to intervention logic.
      3. Agent selection via random.sample breaks Starsim's CRN system.
      4. Re-enrolment is fully disabled (commented out, pending design decision).
      5. base_care_seeking_rate = 0.8 is hardcoded with a TODO comment.
    """
    # from stisim.interventions.product import Product
    # from stisim.interventions.supply import Supply
    # from stisim.interventions.supplies import Supplies
    # from stisim.interventions.prep import SuppliedPrep
    # from stisim.interventions.prep_manager import PrepManager
    #
    # product  = Product(name='lenacapavir', type='prep', delivery_mode='injectable',
    #                    cost=40_000, eff_by_ti=[0.999])
    # supply   = Supply(quantity=1000, product=product)
    # supplies = Supplies(supplies=[supply])
    #
    # prep_manager = PrepManager()   # required; ordering matters; easy to forget
    # prep = SuppliedPrep(
    #     name='lnc_fsw',
    #     eligibilities=[lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected],
    #     coverages=[0.6],
    #     supplies=supplies,
    # )
    # sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net,
    #              interventions=[prep_manager, prep], dt=1/12)
    # sim.run()
    pass


# =============================================================================
# V4 — Supply-constrained allocation (NOTIMPLEMENTED)
#
# Scenario: a fixed annual supply of lenacapavir doses must be allocated across
# two or more target groups. The allocation algorithm distributes doses
# proportionally to the per-group coverage gap:
#
#   gap_i         = (target_coverage_i - actual_coverage_i) * group_size_i
#   doses_group_i = total_supply × gap_i / Σ gap_j
#
# This is the same logic as Clark's SuppliedPrep.step(). The question is
# whether it belongs inside Prep (proposed here) or in a separate Supplies
# helper (Clark's approach). Both work; keeping it in Prep reduces the number
# of classes the user needs to know about.
#
# Edge cases to handle:
#   - All groups at or above target coverage: no doses distributed.
#   - One group has zero gap: all doses go to remaining groups.
#   - supply= is exhausted mid-year: stop distributing, log a warning.
# =============================================================================

def v4_supply_constrained():
    """
    NOTIMPLEMENTED — fixed annual supply of lenacapavir doses, gap-proportional allocation.

    500 doses/year across FSW (target 60%) and AGYW (target 30%).
    If FSW are at 50% and AGYW at 10%, the gaps are 10pp and 20pp respectively,
    so AGYW receive 2/3 of the available supply that timestep.

    The result sim.results.prep.doses_remaining should track supply drawdown.
    """
    # NOTIMPLEMENTED: groups=, product=, supply= on Prep
    #
    # prep = sti.Prep(
    #     groups=[
    #         dict(label='fsw',
    #              eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    #              coverage=0.6),
    #         dict(label='agyw',
    #              eligibility=lambda sim: sim.people.female & ~sim.diseases.hiv.infected
    #                                     & (sim.people.age >= 15) & (sim.people.age < 25),
    #              coverage=0.3),
    #     ],
    #     product=dict(name='lenacapavir', efficacy=0.999, cost_per_year=40_000,
    #                  duration=ss.dur(6, 'months')),
    #     supply=500,   # doses per year; gap-proportional allocation across groups
    # )
    # sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    # sim.run()
    # print('doses used:  ', prep.doses_used)
    # print('total cost:  ', prep.total_cost)
    # print('n FSW on PrEP:', sim.results.prep.fsw.n_on_prep[-1])
    pass


# =============================================================================
# V5 — Kenya lenacapavir CE: price threshold sweep (NOTIMPLEMENTED)
#
# This is the primary research question motivating the whole design effort.
#
# The CEMA Kenya collaborator has asked us to:
#   (a) Estimate the price per person-year at which Kenya would be willing to
#       pay for lenacapavir PrEP (vs. oral PrEP or no PrEP).
#   (b) Evaluate the current lenacapavir rollout in Kenya (FSW-first prioritisation).
#
# Kenya's WTP threshold is approximately $865/DALY (using GNI-per-capita rule),
# though values up to ~$3,000/DALY are sometimes cited for HIV interventions.
# The ICER we want is: incremental cost / infections averted (or DALYs averted,
# which requires additional conversion).
#
# V5a — price sweep: hold coverage constant, vary cost per person-year.
# V5b — targeting comparison: hold cost and total doses constant, vary who gets PrEP.
#
# Both require:
#   - product= cost tracking (V3 feature)
#   - infections-averted vs. a no-PrEP or oral-PrEP baseline
#   - Kenya-calibrated HIV parameters — these should come from hiv_kenya or
#     be passed as explicit pars; see make_kenya_sim() placeholder below.
# =============================================================================

def make_kenya_sim(interventions):
    """
    Placeholder: build a Kenya-calibrated HIV sim. In practice this would use
    hiv_kenya parameters (or a subset) rather than defaults.

    NOTIMPLEMENTED: hiv_kenya localization not yet wired in here.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    return ss.Sim(n_agents=5000, diseases=hiv, networks=net,
                  interventions=interventions, dt=1/12)


def v5a_kenya_price_sweep():
    """
    NOTIMPLEMENTED — sweep lenacapavir cost per person-year; compute ICER vs. oral PrEP.

    We hold coverage and targeting fixed (FSW at 30% — approximate current rollout)
    and vary the cost of lenacapavir from $50 (hypothetical generic price) to
    $40,000 (current Gilead price). For each cost point we compute:

        ICER = (cost_lnc - cost_oral) / (infections_averted_lnc - infections_averted_oral)

    The output is an ICER curve. The intersection with Kenya's WTP threshold
    (~$865/DALY, or a range of plausible thresholds) gives the maximum acceptable
    price per person-year.

    Note on stochasticity: at n_agents=5000 the number of infections averted
    per run is noisy. For publication-quality results we should run n=50–100
    replicates per cost point and use the median ICER. For a first look,
    a single replicate per point is fine.
    """
    # costs_per_year = [50, 200, 500, 1_000, 2_000, 5_000, 10_000, 40_000]
    # fsw_elig = lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected
    #
    # # Baseline: oral PrEP at 30% FSW coverage
    # baseline = make_kenya_sim(interventions=[
    #     sti.Prep(name='oral', coverage=0.3, eligibility=fsw_elig,
    #              product=dict(name='oral', efficacy=0.80, cost_per_year=50)),
    # ])
    # baseline.run()
    # baseline_infections = baseline.results.hiv.new_infections[:].sum()
    # baseline_cost       = baseline.results.oral.total_cost[-1]
    #
    # results = []
    # for cost in costs_per_year:
    #     sim = make_kenya_sim(interventions=[
    #         sti.Prep(name='lnc', coverage=0.3, eligibility=fsw_elig,
    #                  product=dict(name='lenacapavir', efficacy=0.999, cost_per_year=cost,
    #                               duration=ss.dur(6, 'months'))),
    #     ])
    #     sim.run()
    #     infections = sim.results.hiv.new_infections[:].sum()
    #     averted    = baseline_infections - infections
    #     incr_cost  = sim.results.lnc.total_cost[-1] - baseline_cost
    #     results.append(sc.objdict(
    #         cost_per_year = cost,
    #         averted       = averted,
    #         incr_cost     = incr_cost,
    #         icer          = incr_cost / averted if averted > 0 else float('inf'),
    #     ))
    #
    # # Plot: ICER vs. cost/year; horizontal line at Kenya WTP threshold
    # import matplotlib.pyplot as plt
    # costs  = [r.cost_per_year for r in results]
    # icers  = [r.icer          for r in results]
    # plt.semilogx(costs, icers, 'o-')
    # plt.axhline(865, color='red', linestyle='--', label='Kenya WTP (~$865/DALY)')
    # plt.xlabel('Lenacapavir cost per person-year (USD)')
    # plt.ylabel('ICER (USD per infection averted)')
    # plt.title('Kenya: lenacapavir price threshold')
    # plt.legend(); plt.show()
    pass


def v5b_kenya_targeting_comparison():
    """
    NOTIMPLEMENTED — fixed supply of lenacapavir; compare targeting strategies.

    Holds total doses constant at 10,000/year. Asks: which targeting strategy
    averts the most infections for that fixed budget?

    Strategies:
      fsw_only   — 80% coverage of FSW (current rollout priority)
      agyw_only  — 40% coverage of AGYW (alternative targeting)
      msm_only   — 60% coverage of MSM (requires AgeMatchedMSM network)
      pbfw_only  — 60% coverage of PBFW (requires Pregnancy module)
      fsw_agyw   — 50% FSW + 30% AGYW (split; gap-proportional allocation)

    MSM eligibility requires sim.networks.agematchedmsm to be present.
    PBFW eligibility requires sim.demographics.pregnancy to be present.
    If those modules are absent, the lambdas raise AttributeError — Prep
    should validate eligibility groups at init_pre time and raise a clear error.
    """
    # fsw_elig  = lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected
    # agyw_elig = lambda sim: (sim.people.female & ~sim.diseases.hiv.infected
    #                          & (sim.people.age >= 15) & (sim.people.age < 25))
    # msm_elig  = lambda sim: sim.networks.agematchedmsm.msm & ~sim.diseases.hiv.infected
    # pbfw_elig = lambda sim: sim.demographics.pregnancy.pregnant & ~sim.diseases.hiv.infected
    #
    # strategies = {
    #     'fsw_only':  [dict(label='fsw',  eligibility=fsw_elig,  coverage=0.8)],
    #     'agyw_only': [dict(label='agyw', eligibility=agyw_elig, coverage=0.4)],
    #     'msm_only':  [dict(label='msm',  eligibility=msm_elig,  coverage=0.6)],
    #     'pbfw_only': [dict(label='pbfw', eligibility=pbfw_elig, coverage=0.6)],
    #     'fsw_agyw':  [dict(label='fsw',  eligibility=fsw_elig,  coverage=0.5),
    #                   dict(label='agyw', eligibility=agyw_elig, coverage=0.3)],
    # }
    # results = {}
    # for label, groups in strategies.items():
    #     prep = sti.Prep(
    #         groups=groups,
    #         product=dict(name='lenacapavir', efficacy=0.999, cost_per_year=40_000,
    #                      duration=ss.dur(6, 'months')),
    #         supply=10_000,
    #     )
    #     sim = make_kenya_sim(interventions=[prep])
    #     sim.run()
    #     results[label] = sim.results.hiv.new_infections[:].sum()
    #
    # sc.pp(results)
    pass


# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    sc.heading('V1: simple FSW PrEP (runs today)')
    v1_fsw()
    v1b_agyw()

    sc.heading('V2a: multi-group via multiple instances (runs today)')
    v2a_multigroup_instances()

    sc.heading('V2b–V5: API sketches (NOTIMPLEMENTED — see comments and module docstring)')
