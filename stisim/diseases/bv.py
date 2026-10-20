"""
BV modules
Includes:
    - a simple BV model for generating background prevalence of vaginal discharge
    - a detailed model of the vaginal microbiome including community state types (CSTs)
"""

import numpy as np
import starsim as ss
import sciris as sc
import stisim.utils as ut
from stisim.diseases.sti import BaseSTI

__all__ = ['SimpleBV', 'BV']


class SimpleBV(ss.Disease):
    """
    Simple bacterial vaginosis model for generating background BV prevalence.

    Models BV as a spontaneous condition (not sexually transmitted) using a
    logistic regression model with risk factors (douching, menstrual hygiene,
    number of partners). Includes symptomatic/asymptomatic states, care-seeking,
    and natural clearance or persistence.

    Args:
        pars (dict): Override default parameters.
        name (str): Module name. Default: ``'bv'``.
        **kwargs: Additional parameters passed to ``update_pars``.
    """

    def __init__(self, pars=None, name='bv', **kwargs):
        super().__init__(name=name)

        self.define_pars(
            dt='month',
            include_care=True,
            p_symp=ss.bernoulli(p=0.5),
                    dur_presymp=ss.uniform(ss.weeks(1), ss.weeks(8)),  # Duration of presymptomatic period
        dur_asymp2clear=ss.uniform(ss.weeks(1), ss.weeks(18)),  # Duration of asymptomatic infection
        dur_symp2clear=ss.uniform(ss.weeks(1), ss.weeks(18)),  # Duration of symptoms
            init_prev=ss.bernoulli(p=0.025),

            # Spontaneous occurrence parameters. These will be used within a logistic regression
            # model to calculate the probability of spontaneous occurrence. The model is flexible
            # but should always include an intercept term.
            p_bv=ss.bernoulli(p=0.01),  # Probability of BV in the general population. Overwritten by the model below
            p_douching=ss.bernoulli(p=0),  # Share of population douching
            p_poor_menstrual_hygiene=ss.bernoulli(p=0),  # Share of population with poor menstrual hygiene
            p_base=0.45,                 # Used to calculate the baseline (intercept) probability of spontaneous occurrence
            p_spontaneous=sc.objdict(
                douching=3,             # OR of BV for douching
                poor_menstrual_hygiene=2,    # OR of BV for poor menstrual hygiene
            ),

            # Care-seeking and clearance
            dur_symp2care=ss.lognorm_ex(ss.weeks(4), ss.weeks(4)),
            p_symp_care=ss.bernoulli(p=0.3),
            p_clear=ss.bernoulli(p=0.5),
            dur_persist=ss.constant(ss.years(100)),
        )
        self.update_pars(pars, **kwargs)

        # States that elevate risk of BV
        self.define_states(
            ss.BoolState('susceptible', default=True, label='Susceptible'),
            ss.BoolState('bv_prone', label='Prone to BV', default=ss.bernoulli(p=0.5)),  # Percentage of women "prone" to BV
            ss.BoolState('infected', label='Infected'),
            ss.BoolState('asymptomatic', label='Asymptomatic'),
            ss.BoolState('symptomatic', label='Symptomatic'),
            ss.FloatArr('rel_sus', default=1.0, label='Relative susceptibility'),
            ss.FloatArr('rel_trans', default=1.0, label='Relative transmissibility'),  # NOT USED
            ss.FloatArr('ti_seeks_care', label='Time of care seeking'),
            ss.BoolArr('seeking_care', label='Seeking care'),
            ss.FloatArr('ti_infected', label='Time of infection'),
            ss.FloatArr('ti_clearance', label='Time of clearance'),
            ss.FloatArr('ti_symptomatic', label='Time of symptoms'),
            ss.FloatArr('dur_inf', label='Duration of infection'),
            ss.BoolArr('douching', label='Douching'),
            ss.BoolArr('poor_menstrual_hygiene', label='Poor menstrual hygiene'),
        )
        self.sex_keys = {'': 'alive', 'f': 'female', 'm': 'male'}

        return

    @property
    def treatable(self):
        """ Responds to treatment """
        return self.infected

    def init_results(self):
        """ Initialize results """
        super().init_results()
        results = [
            ss.Result('prevalence', scale=False, label="Prevalence"),
            ss.Result('symp_prevalence', scale=False, label="Symptomatic prevalence", auto_plot=False),
            ss.Result('incidence', scale=False, label="Incidence"),
            ss.Result('new_infections', dtype=int, label="New infections"),
            ss.Result('new_symptomatic', dtype=int, label="New symptomatic", auto_plot=False),
        ]

        if self.pars.include_care:
            sexkeys = ['', 'f', 'm']
            for sk in sexkeys:
                skk = '' if sk == '' else f'_{sk}'
                skl = '' if sk == '' else f' - {sk.upper()}'
                results += [
                    ss.Result('new_care_seekers'+skk, dtype=int, label="New care seekers"+skl, auto_plot=False),
                    ss.Result('new_false_pos'+skk, dtype=int, label="New false positives"+skl, auto_plot=False),
                    ss.Result('new_true_pos'+skk, dtype=int, label="New true positives"+skl, auto_plot=False),
                    ss.Result('new_false_neg'+skk, dtype=int, label="New false negatives"+skl, auto_plot=False),
                    ss.Result('new_true_neg'+skk, dtype=int, label="New true negatives"+skl, auto_plot=False),
                    ss.Result('new_treated_success'+skk, dtype=int, label="Successful treatments"+skl, auto_plot=False),
                    ss.Result('new_treated_failure'+skk, dtype=int, label="Unsuccessful treatments"+skl, auto_plot=False),
                    ss.Result('new_treated_unnecessary'+skk, dtype=int, label="Unnecessary treatments"+skl, auto_plot=False),
                    ss.Result('new_treated'+skk, dtype=int, label="Treatments"+skl, auto_plot=False),
                ]

        self.define_results(*results)
        return

    def _get_uids(self, upper_age=None):
        """ Get uids of females younger than upper_age """
        people = self.sim.people
        if upper_age is None: upper_age = 1000
        within_age = people.age < upper_age
        return (within_age & people.female).uids

    def bv_sus(self):
        return self.sim.people.female & (self.sim.people.age > 15) & self.susceptible & self.bv_prone

    def set_hygiene_states(self, upper_age=None):
        """ Set vaginal hygiene states """
        f_uids = self._get_uids(upper_age=upper_age)
        self.douching[f_uids] = self.pars.p_douching.rvs(f_uids)
        self.poor_menstrual_hygiene[f_uids] = self.pars.p_poor_menstrual_hygiene.rvs(f_uids)
        return

    def init_post(self):
        """ Initialize with sim properties """
        ss.Module.init_post(self)  # Skip the disease init_post() since we create infections in a different way
        self.set_hygiene_states()
        self.infect()
        return

    def spontaneous(self, uids):
        """ Create new cases via spontaneous occurrence """
        # Set intercept
        p = sc.dcp(self.pars.p_spontaneous)
        intercept = -np.log(1/self.pars.p_base-1)  # Use a transformation consistent with the logistic regression
        rhs = np.full_like(uids, fill_value=intercept, dtype=float)

        # Add all covariates
        for term, val in p.items():
            rhs += val * getattr(self, term)[uids]

        # Calculate the probability of spontaneous occurrence of BV
        p_bv = 1 / (1+np.exp(-rhs))

        return p_bv

    def infect(self):
        # Create new cases via spontaneous occurrence
        uids = self.bv_sus().uids
        p_bv = self.spontaneous(uids)
        self.pars.p_bv.set(p_bv)
        bv_cases = self.pars.p_bv.filter(uids)

        # Set prognoses
        if len(bv_cases):
            self.set_prognoses(bv_cases)
        return

    def step(self):
        self.set_hygiene_states()
        self.infect()

    def clear_infection(self, uids):
        self.infected[uids] = False
        self.symptomatic[uids] = False
        self.asymptomatic[uids] = False
        self.susceptible[uids] = True
        self.ti_clearance[uids] = self.ti

    def step_state(self):
        """ Updates for this timestep """
        ti = self.ti

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Presymptomatic -> symptomatic
        new_symptomatic = (self.asymptomatic & (self.ti_symptomatic <= ti)).uids
        if len(new_symptomatic):
            self.asymptomatic[new_symptomatic] = False
            self.symptomatic[new_symptomatic] = True
            self.ti_symptomatic[new_symptomatic] = ti

        # Clear infections
        new_cleared = (self.infected & (self.ti_clearance <= ti)).uids
        self.clear_infection(new_cleared)

        # Symptomatic/PID care seeking
        old_seekers = self.seeking_care
        self.seeking_care[old_seekers] = False
        self.ti_seeks_care[old_seekers] = np.nan  # Remove the old
        new_seekers = (self.infected & (self.ti_seeks_care <= ti)).uids
        self.seeking_care[new_seekers] = True
        self.ti_seeks_care[new_seekers] = ti
        if self.pars.include_care:
            self.results['new_care_seekers'][ti] = np.count_nonzero(self.ti_seeks_care == ti)

        return

    def wipe_dates(self, uids):
        """ Clear all previous dates """
        self.ti_infected[uids] = np.nan
        self.ti_symptomatic[uids] = np.nan
        self.ti_clearance[uids] = np.nan
        self.dur_inf[uids] = np.nan
        return

    def set_infection(self, uids):
        self.susceptible[uids] = False
        self.infected[uids] = True
        self.asymptomatic[uids] = True
        self.ti_infected[uids] = self.ti
        return

    def set_symptoms(self, uids):
        p = self.pars
        symp, asymp = p.p_symp.split(uids)
        dur_presymp = self.pars.dur_presymp.rvs(symp)
        self.ti_symptomatic[symp] = self.ti_infected[symp] + dur_presymp
        return symp, asymp

    def set_care_seeking(self, symp):
        symp_care = self.pars.p_symp_care.filter(symp)
        dur_symp2care = self.pars.dur_symp2care.rvs(symp_care)
        self.ti_seeks_care[symp_care] = self.ti_symptomatic[symp_care] + dur_symp2care
        return

    def set_duration(self, symp, asymp):
        """ Set duration of infection """
        # Extract people who have persistent infection and set their clearance time
        clear, persist = self.pars.p_clear.split(symp | asymp)
        self.ti_clearance[persist] = self.ti_infected[persist] + self.pars.dur_persist.rvs(persist)

        # Set clearance times for the rest
        symp_clear = symp.intersect(clear)
        asymp_clear = asymp.intersect(clear)
        if len(symp_clear):
            dur_inf_symp = self.pars.dur_symp2clear.rvs(symp_clear)
            self.ti_clearance[symp_clear] = dur_inf_symp + self.ti_symptomatic[symp_clear]
        if len(asymp_clear):
            dur_inf_asymp = self.pars.dur_asymp2clear.rvs(asymp_clear)
            self.ti_clearance[asymp_clear] = dur_inf_asymp + self.ti_infected[asymp_clear]

        return

    def set_prognoses(self, uids, source_uids=None):
        """
        Set initial prognoses for adults newly infected
        """
        self.wipe_dates(uids)  # Clear prior dates
        self.set_infection(uids)  # Set infection
        symp, asymp = self.set_symptoms(uids)  # Set symptoms & presymptomatic duration
        self.set_care_seeking(symp)  # Determine who seeks care and when
        self.set_duration(symp, asymp)

        # Determine overall duration of infection
        self.dur_inf[uids] = self.ti_clearance[uids] - self.ti_infected[uids]

        if (self.dur_inf[uids] < 0).any():
            errormsg = 'Invalid durations of infection'
            raise ValueError(errormsg)

        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        women = (self.sim.people.age >= 15) & self.sim.people.female
        self.results['prevalence'][ti] = ut.cond_prob(self.infected, women)
        self.results['symp_prevalence'][ti] = ut.cond_prob(self.symptomatic, women)
        self.results['new_infections'][ti] = ut.count(self.ti_infected == (ti - 1))
        self.results['new_symptomatic'][ti] = ut.count(self.ti_symptomatic == (ti - 1))

        return


adult_f = ss.ppl.female & (ss.ppl.age >= 15)
stable = lambda *vals: ss.by(stable=list(vals)) # A value for each stable CST: 1, 3, 4


class BV(BaseSTI):
    """
    Bacterial vaginosis, as the community state type (CST) of the vaginal microbiome

    Women enter the model at 15 (`debut`), at their stable CST: CST 1 (*L. crispatus*),
    CST 3 (*L. iners*), or CST 4 (BV, anaerobes). Each woman's stable CST (the `stable`
    dimension) sets her risk of worsening, her chance of symptoms, and the CST she returns to
    after BV clears.

    States: U (no CST: men and girls), CST1, CST3, CST4a (asymptomatic BV), CST4p
    (presymptomatic), CST4s (symptomatic), and P (after treatment).

        U -> CST1 | CST3 | CST4p | CST4a           on turning 15, by stable CST
        CST1 -> CST3                               p_cst_change_cst1, per month (shift13)
        CST3 -> CST1                               after dur2clear_cst3, in stable-CST1 women (revert31)
        CST3 -> CST4p (p_symp) | CST4a             p_cst_change_cst3, per month (onset)
        CST4p -> CST4s                             after dur_presymp
        CST4s, CST4a -> CST1 (stable 1) | CST3     after dur2clear_cst4
        CST4* -> P                                 treatment (tx): moved only by sti.treat_BV
        P -> CST1 (stable 1) | CST3                after dur_post_tx (settle)
        P -> CST4p (p_symp) | CST4a                recur_rate*exp(-P.dwell/dur_recur_decay) (recur)

    The rates of worsening (shift13, onset, and recur) are multiplied by the risk factors,
    each a trait drawn once for each woman and a dimension of the disease: the stable CST
    (`rr_stable_cst1/3/4`), douching (`p_douching`, `rr_douching`), poor menstrual hygiene
    (`p_poor_menstrual_hygiene`, `rr_poor_menstrual_hygiene`), and concurrent partners
    (`n_partners`, `rr_concurrency`; with `count_partners=True`, from the sexual network
    instead). Onset and recurrence are also multiplied by `rr_uncircumcised` for women
    with an uncircumcised partner in a sexual network (`p_circumcised` of men).

    With the default `dur_post_tx` of 0, treated women return to their CST at once, as in
    STIsim 3; a longer window (e.g. `ss.months(6)`) models recurrence after treatment.

    Results (among women 15+): `prevalence` (CST 4), `symp_prevalence`, `cst1_prevalence`,
    `cst3_prevalence`, `cst4_prevalence`, `n_infected`, `new_infections` (onset and
    recurrence), `new_symptomatic`, and the counts and flows of every state and transition,
    with `.select()` and `.by()` by stable CST, risk factor, and the sim's age groups.

    Args:
        name (str): the module name
        kwargs (dict): parameter values; those of the risk factors and the stable CST are read when the module is made

    Examples:
        ```python
        bv = sti.BV(p_poor_menstrual_hygiene=0.3, dur_post_tx=ss.months(6))
        sim = sti.Sim(diseases=bv, interventions=sti.treat_BV())
        ```
    """
    U, CST1, CST3, CST4a, CST4p, CST4s, P = ss.states('U CST1 CST3 CST4a CST4p CST4s P')
    transmits = False

    # Natural history
    p_cst_change_cst1 = ss.par(ss.probpermonth(0.10), doc='Rate of change from CST 1 to CST 3')
    p_cst_change_cst3 = ss.par(ss.probpermonth(0.05), doc='Rate of change from CST 3 to CST 4 (BV)')
    p_symp = ss.par(stable(0.8, 0.7, 0.1), doc='Probability of symptomatic BV, by stable CST')
    dur_presymp = ss.par(ss.uniform(ss.weeks(1), ss.weeks(2)), doc='Time from BV to symptoms, for those who develop them')
    dur2clear_cst3 = ss.par(ss.uniform(ss.weeks(1), ss.weeks(18)), doc='Duration of a transient CST 3, in stable-CST1 women')
    dur2clear_cst4 = ss.par(ss.uniform(ss.weeks(5), ss.weeks(50)), doc='Duration of untreated BV (https://doi.org/10.1016/S0002-9378(03)01042-1)')
    dur_post_tx = ss.par(ss.constant(0), doc='Time after treatment before return to the stable CST (0: at once)')
    recur_rate = ss.par(ss.permonth(0.06), doc='Rate of recurrence just after treatment')
    dur_recur_decay = ss.par(ss.months(3), doc='Time over which the rate of recurrence falls by a factor of e')

    # Risk factors: read when the module is made (they are the adjustments of its dimensions)
    stable_cst_distribution = ss.par([0.10, 0.30, 0.60], doc='Shares with stable CST 1, 3, 4 (https://doi.org/10.1186/s40168-021-01096-9)')
    rr_stable_cst1 = ss.par(0.1, doc='Relative rate of worsening, stable CST 1')
    rr_stable_cst3 = ss.par(1.0, doc='Relative rate of worsening, stable CST 3 (reference)')
    rr_stable_cst4 = ss.par(2.0, doc='Relative rate of worsening, stable CST 4')
    p_douching = ss.par(0.64, doc='Share of women douching (Nigeria, Ghana)')
    rr_douching = ss.par(1.21, doc='Relative rate, douching (https://pmc.ncbi.nlm.nih.gov/articles/PMC2574994/)')
    p_poor_menstrual_hygiene = ss.par(0.55, doc='Share of women with poor menstrual hygiene (https://pmc.ncbi.nlm.nih.gov/articles/PMC9817285/)')
    rr_poor_menstrual_hygiene = ss.par(4.1, doc='Relative rate, poor menstrual hygiene (Kenya calibration)')
    n_partners = ss.par(0.05, doc='Mean number of concurrent partners beyond the first (Poisson), for the share with concurrent partners')
    count_partners = ss.par(False, doc='Count concurrent partners from the sexual network instead of drawing them')
    rr_concurrency = ss.par(1.28, doc='Relative rate, concurrent partners (https://pmc.ncbi.nlm.nih.gov/articles/PMC5429208/)')
    p_circumcised = ss.par(0.4, doc='Share of men who are circumcised')
    rr_uncircumcised = ss.par(4.0, doc='Relative rate of BV onset with an uncircumcised partner (Kenya calibration)')

    # Transitions
    debut    = ss.progression('U -> CST1 | CST3 | CST4p | CST4a', on=adult_f, p=[stable(1, 0, 0), stable(0, 1, 0), stable(0, 0, 1)*p_symp, None])
    shift13  = ss.Transition(CST1, CST3, rate=p_cst_change_cst1)
    revert31 = ss.Transition(CST3, CST1, dur=dur2clear_cst3)
    onset    = ss.Transition(CST3, [CST4p, CST4a], rate=p_cst_change_cst3, p=[p_symp, None])
    symptoms = ss.Transition(CST4p, CST4s, dur=dur_presymp)
    clear_s  = ss.Transition(CST4s, [CST1, CST3], dur=dur2clear_cst4, p=[stable(1, 0, 0), None])
    clear_a  = ss.Transition(CST4a, [CST1, CST3], dur=dur2clear_cst4, p=[stable(1, 0, 0), None])
    tx       = ss.Transition(CST4a + CST4p + CST4s, P, undriven=True)
    settle   = ss.Transition(P, [CST1, CST3], dur=dur_post_tx, p=[stable(1, 0, 0), None])
    recur    = ss.Transition(P, [CST4p, CST4a], rate=recur_rate*np.exp(-P.dwell/dur_recur_decay), p=[p_symp, None])

    # Results, among women 15+ (everyone with a CST)
    n_infected       = CST4a + CST4p + CST4s
    n_women          = CST1 + CST3 + CST4a + CST4p + CST4s + P
    prevalence       = (CST4a + CST4p + CST4s)/(CST1 + CST3 + CST4a + CST4p + CST4s + P)
    cst4_prevalence  = (CST4a + CST4p + CST4s)/(CST1 + CST3 + CST4a + CST4p + CST4s + P)
    cst1_prevalence  = CST1/(CST1 + CST3 + CST4a + CST4p + CST4s + P)
    cst3_prevalence  = CST3/(CST1 + CST3 + CST4a + CST4p + CST4s + P)
    symp_prevalence  = CST4s/(CST1 + CST3 + CST4a + CST4p + CST4s + P)
    new_infections   = onset + recur
    new_symptomatic  = ss.outcome(symptoms)

    def __init__(self, name='bv', **kwargs):
        for key,lam in [('n_partners', 'lam'), ('stable_cst_distribution', 'p')]: # Starsim 3 distributions
            val = kwargs.get(key)
            if isinstance(val, ss.Dist):
                kwargs[key] = val.pars[lam]
        super().__init__(name=name, **kwargs)
        v = self.value
        rr = lambda key: [1, v(key)]
        worsen = ['shift13', 'onset', 'recur']
        p_conc = 1 - np.exp(-v('n_partners'))
        self.stratify(stable=['1', '3', '4'], split=list(v('stable_cst_distribution')),
                      adjust={k:[v('rr_stable_cst1'), v('rr_stable_cst3'), v('rr_stable_cst4')] for k in worsen} | dict(revert31=[1, 0, 0])) # A transient CST 3 only in stable-CST1 women
        for dim,p_key,rr_key in [('douching', 'p_douching', 'rr_douching'), ('poor_menstrual_hygiene', 'p_poor_menstrual_hygiene', 'rr_poor_menstrual_hygiene'),
                                 ('concurrency', None, 'rr_concurrency')]:
            p = p_conc if p_key is None else v(p_key)
            if dim == 'concurrency' and v('count_partners'):
                continue # From the network instead: see modifiers()
            if v(rr_key) != 1 and 0 < p < 1: # A risk factor that changes something: a dimension
                self.stratify(**{dim:[False, True]}, split=[1 - p, p], adjust={k:rr(rr_key) for k in worsen})
        self.define_states(ss.BoolArr('circumcised', default=ss.bernoulli(v('p_circumcised')), label='Circumcised (men)'))
        return

    def modifiers(self, sim):
        """ Risk factors that come from the sexual network: an uncircumcised partner, and (with count_partners) concurrent partners """
        nets = [name for name,net in sim.networks.items() if isinstance(net, ss.SexualNetwork)]
        if not nets:
            return []
        name = self.name
        v = self.value
        def uncirc_partner(sim):
            """ Women with an uncircumcised partner """
            bv = sim.diseases[name]
            out = np.zeros(len(sim.people.uid.raw), dtype=bool)
            uncirc = (sim.people.male & ~bv.circumcised).uids
            for net in nets:
                out[ss.uids(sim.networks[net].find_contacts(uncirc))] = True
            return out
        def concurrent(sim):
            """ Women with two or more partners in the sexual networks """
            out = np.zeros(len(sim.people.uid.raw), dtype=bool)
            auids = sim.people.auids
            n = np.zeros(len(auids))
            for net in nets:
                n += sim.networks[net].degrees(auids) if hasattr(sim.networks[net], 'degrees') else sim.networks[net].partners[auids]
            out[auids] = n > 1
            return out
        mods = [ss.modify(getattr(self, k), multiply=v('rr_uncircumcised'), where=uncirc_partner, name=f'{name}_uncirc_{k}') for k in ['onset', 'recur']]
        if v('count_partners'):
            mods += [ss.modify(getattr(self, k), multiply=v('rr_concurrency'), where=concurrent, name=f'{name}_concurrency_{k}') for k in ['shift13', 'onset', 'recur']]
        return mods

    def flat_results(self):
        return [] # Women only: no results by sex

    # States as STIsim 3 named them, as agents, for connectors and custom code
    @property
    def cst1(self): return self.CST1
    @property
    def cst3(self): return self.CST3
    @property
    def cst4(self): return self.CST4a | self.CST4p | self.CST4s
    @property
    def infected(self): return self.cst4
    @property
    def symptomatic(self): return self.CST4s
    @property
    def asymptomatic(self): return self.CST4a | self.CST4p
    @property
    def susceptible(self): return self.CST1 | self.CST3 | self.P
