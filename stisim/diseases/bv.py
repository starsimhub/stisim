"""
Bacterial vaginosis (BV) models:
    - SimpleBV: a spontaneous (not sexually transmitted) condition, for a background prevalence of vaginal discharge
    - BV: a model of the vaginal microbiome by community state type (CST)
"""

import functools
import numpy as np
import starsim as ss
import sciris as sc
import stisim.utils as ut
from stisim.networks.base import BaseNetwork
from stisim.diseases.sti import STIResults, sexual_beta, default_age_bins, default_sex_keys

__all__ = ['SimpleBV', 'BV']


class SimpleBV(ss.Disease, STIResults):
    """
    Bacterial vaginosis as a spontaneous condition, for a background prevalence of vaginal discharge, declared

    BV is not sexually transmitted here. Each month, women over 15 who are prone to BV (`bv_prone`, a share
    `p_bv_prone` drawn at birth) develop it with a probability from a logistic model: an intercept of
    `p_base`, raised by douching and poor menstrual hygiene (`p_spontaneous`, odds ratios). A share `p_symp`
    become symptomatic after `dur_presymp`; the rest stay asymptomatic. A share `p_symp_care` of the
    symptomatic seek care, `dur_symp2care` after their symptoms start. Infections clear `dur_symp2clear`
    after symptoms start, or `dur_asymp2clear` after onset without symptoms, except a share `1 - p_clear`
    that persists for `dur_persist`.

    States: `susceptible`, `presymptomatic`, `asymptomatic`, `symptomatic`; `infected` is a condition over
    the last three.

    Transitions (each is an event that interventions can be triggered by, e.g. `who=bv.seek_care`):

    - `onset`: spontaneous BV, to `presymptomatic` (with probability `p_symp`) or `asymptomatic`
    - `symptoms`: `dur_presymp` after onset
    - `seek_care`: a share `p_symp_care` of the symptomatic seek care, `dur_symp2care` after their symptoms start; they stay symptomatic
    - `clear_symp`, `clear_asymp`: clearance; `persist_symp`, `persist_asymp`: the share `1 - p_clear` that does not clear

    The results are those of the other STIs (`sti.STIResults`): `prevalence` is over adults of both sexes,
    so among women it is `prevalence_f`.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'bv'``)
        age_bins (list): bin edges for age-stratified results; ``None`` skips them
        sex_keys (dict): sex strata for results; ``None`` skips sex-stratified results
        include_care (bool): whether to define the results of care
        kwargs (dict): parameter values, or `label`
    """
    suggested_solver = ss.solvers.abm(step='wholestep') # As agents, rates over whole steps, as the other STIs
    def __init__(self, pars=None, name='bv', age_bins=default_age_bins, sex_keys=default_sex_keys, include_care=True, **kwargs):

        # Natural history
        p_symp = ss.par(0.5, doc='Probability that BV becomes symptomatic')
        dur_presymp = ss.par(ss.uniform(ss.weeks(1), ss.weeks(8)), doc='Time from onset to symptoms')
        p_symp_care = ss.par(0.3, doc='Probability that a woman with symptoms seeks care')
        dur_symp2care = ss.par(ss.lognormal(ss.weeks(4), ss.weeks(4)), name='dur_symp2care', doc='Time from symptoms to seeking care')
        dur_symp2clear = ss.par(ss.uniform(ss.weeks(1), ss.weeks(18)), doc='Time from symptoms to clearance')
        dur_asymp2clear = ss.par(ss.uniform(ss.weeks(1), ss.weeks(18)), doc='Time from onset to clearance, without symptoms')
        p_clear = ss.par(0.5, doc='Probability that BV clears (the rest persist)')
        dur_persist = ss.par(ss.constant(ss.years(100)), doc='Duration of persistent BV')

        # States
        susceptible, presymptomatic, asymptomatic, symptomatic = ss.states('susceptible presymptomatic asymptomatic symptomatic')

        # Spontaneous onset, with a probability per month from the risk factors (see p_onset())
        onset = ss.Transition(susceptible, {presymptomatic: p_symp, asymptomatic: None}, rate=functools.partial(self.p_onset))
        symptoms = presymptomatic >> symptomatic @ dur_presymp
        seek_care = ss.Transition(symptomatic, {symptomatic: p_symp_care, ss.stay: None}, when=ss.on(symptoms, delay=dur_symp2care)) # A share seek care, and stay symptomatic
        clear_symp = ss.Transition(symptomatic, susceptible, dur=dur_symp2clear, p=None) # Takes those who do not persist
        persist_symp = ss.Transition(symptomatic, susceptible, p=1 - p_clear, dur=dur_persist)
        clear_asymp = ss.Transition(asymptomatic, susceptible, dur=dur_asymp2clear, p=None)
        persist_asymp = ss.Transition(asymptomatic, susceptible, p=1 - p_clear, dur=dur_persist)

        # Parameters read by its own code, and (as for sti.SEIS) the time after an event, set once the graph is built
        values = sc.mergedicts(dict(init_prev=0), pars, kwargs) # Nobody is seeded: onset starts on the first step
        code_pars = dict(
            p_bv_prone = 0.5, # Share of women prone to BV
            p_douching = ss.bernoulli(p=0), # Share of women douching
            p_poor_menstrual_hygiene = ss.bernoulli(p=0), # Share of women with poor menstrual hygiene
            p_base = 0.45, # Probability per month of BV for a prone woman without risk factors (the intercept of the logistic model)
            p_spontaneous = sc.objdict(douching=3, poor_menstrual_hygiene=2), # Odds ratios of BV for each risk factor
        )
        delays = dict(dur_symp2care=dur_symp2care)
        given = {key:values.pop(key) for key in list(code_pars) + list(delays) if key in values}

        super().__init__(
            onset=onset, symptoms=symptoms, seek_care=seek_care,
            clear_symp=clear_symp, persist_symp=persist_symp, clear_asymp=clear_asymp, persist_asymp=persist_asymp,
            new_symptomatic=1*symptoms, # Result of the flow
            name=name, **sc.mergedicts(dict(dt='month'), values), # Monthly unless dt is given
        )
        self.define_pars(**code_pars, **delays)
        self.pars.update(given)

        # Each woman's own risk factors
        self.define_states(
            ss.BoolArr('bv_prone', default=ss.bernoulli(p=self.pars.p_bv_prone), label='Prone to BV'),
            ss.BoolArr('douching', label='Douching'),
            ss.BoolArr('poor_menstrual_hygiene', label='Poor menstrual hygiene'),
        )
        self.set_strata([15, 100], age_bins, sex_keys)
        self.include_care = include_care
        return

    @property
    def infected(self):
        """ Women with BV """
        return self.presymptomatic | self.asymptomatic | self.symptomatic

    @property
    def treatable(self):
        """ Responds to treatment """
        return self.infected

    def p_onset(self, sim, uids):
        """ The rate of spontaneous BV: a probability per month for prone women over 15, from the logistic model """
        p = self.pars
        rhs = np.full(len(uids), -np.log(1/p.p_base - 1))
        for term, val in p.p_spontaneous.items():
            rhs += val*getattr(self, term)[uids]
        prob = 1/(1 + np.exp(-rhs))
        eligible = sim.people.female[uids] & (sim.people.age[uids] > 15) & self.bv_prone[uids]
        return ss.probpermonth(np.where(eligible, np.minimum(prob, 1 - 1e-9), 0.0))

    def set_hygiene_states(self):
        """ Draw each woman's vaginal hygiene risk factors (each step, as they may change) """
        f_uids = self.sim.people.female.uids
        self.douching[f_uids] = self.pars.p_douching.rvs(f_uids)
        self.poor_menstrual_hygiene[f_uids] = self.pars.p_poor_menstrual_hygiene.rvs(f_uids)
        return

    def clear_infection(self, uids):
        """ Cure people (e.g. treated): they are susceptible again """
        self.move(self.susceptible, uids=uids)
        return

    def step_state(self):
        """ At the start of the step: the risk factors (the declared transitions are taken in step()) """
        self.set_hygiene_states()
        return

    def new_infected(self):
        """ Who developed BV in this step """
        return self.sim.people.uid.isin(self.onset.uids)

    def init_results(self):
        """ Initialize results """
        super().init_results()
        self.define_sti_results(include_care=self.include_care)
        self.define_results(ss.Result('symp_prevalence', scale=False, label='Symptomatic prevalence', auto_plot=False))
        return

    def record_states(self):
        super().record_states()
        t = self.tally()
        self.record_sti_stocks(t)
        women = t.vals(self.sim.people.female) & (np.asarray(self.sim.people.age.values) >= 15)
        self.results.symp_prevalence[self.ti] = t.div(t.count(t.vals(self.symptomatic) & women), t.count(women))
        return

    def record_events(self):
        super().record_events()
        t = self.tally()
        self.record_sti_flows(t)
        if self.include_care:
            seekers = t.events(self.seek_care)
            for skk, skl, attr in self.sex_strata():
                self.results[f'new_care_seekers{skk}'][self.ti] = t.count(seekers & t.vals(self.sim.people[attr]))
        return

    def finalize_results(self):
        super().finalize_results()
        self.finalize_sti_results()
        return


class BV(ss.Disease, STIResults):
    """
    Bacterial vaginosis, by the community state type (CST) of the vaginal microbiome, declared

    Each woman has a stable CST (1, 3, or 4), drawn at birth, and is in CST 1, CST 3, or CST 4 (BV) at any time.
    Women over 15 move from CST 1 to CST 3, and from CST 3 to CST 4, with a probability per month
    (`p_cst_change`) raised by behavioral risk factors (douching, poor menstrual hygiene, concurrent partners,
    and, in CST 3, an uncircumcised partner) and by her stable CST. CST 4 is symptomatic for a share
    (`p_symp`, by stable CST) after `dur_presymp`, and clears after `dur_cst4` (from the symptoms, for those
    with them) to the CST of `spontaneous_clearance`; a woman whose stable CST is 1 returns from CST 3 to it
    after `dur_cst3`. BV-associated microbes are sexually transmitted if `beta_m2f` is given: men carry them
    (`colonized`) for `dur_male`.

    The CSTs and their transitions are declared; the risk factors are each person's own (arrays), and the
    probabilities of a worse CST that they give are worked out by code (`p_worse()`), as agents. The values that
    depend only on the stable CST, which each woman has from birth, are marked `ss.lifelong()`, so they are worked
    out once for each woman.

    States: `susceptible` (men; and girls, at birth, before they take their stable CST), `cst1`, `cst3`, the
    states of CST 4 (`presymptomatic`, `asymptomatic`, `symptomatic`), and `colonized` (men). `cst4` and
    `infected` are conditions over them.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'bv'``)
        age_bins (list): bin edges for age-stratified results; ``None`` skips them
        sex_keys (dict): sex strata for results; ``None`` skips sex-stratified results
        include_care (bool): whether to define the results of care
        kwargs (dict): parameter values, or `label`
    """
    suggested_solver = ss.solvers.abm(step='wholestep') # As agents, rates over whole steps, as STIsim's diseases took them when written as code (faster; a sim's solver=ss.solvers.abm(step='midstep') overrides it)
    def __init__(self, pars=None, name='bv', age_bins=default_age_bins, sex_keys=default_sex_keys, include_care=True, **kwargs):

        # Natural history
        dur_cst3 = ss.par(ss.uniform(ss.weeks(1), ss.weeks(18)), doc='Time in CST 3 before returning to a stable CST 1')
        dur_cst4 = ss.par(ss.uniform(ss.weeks(5), ss.weeks(50)), doc='Time in CST 4 before it clears, from the symptoms for those with them') # https://www.sciencedirect.com/science/article/abs/pii/S0002937803010421
        dur_presymp = ss.par(ss.uniform(ss.weeks(1), ss.weeks(2)), doc='Time from CST 4 to symptoms')

        # Sexual transmission of BV-associated microbes, if beta_m2f is given
        beta_m2f = ss.par(0.0, doc='Probability of transmission per act from a man to a woman')
        rel_beta_f2m = ss.par(0.5, doc='Transmission from a woman to a man, relative to from a man to a woman')
        dur_male = ss.par(ss.uniform(ss.weeks(1), ss.weeks(18)), doc='Time that a man carries BV-associated microbes')

        # States
        partial = functools.partial
        lifelong = lambda func, **kw: ss.lifelong(partial(func, **kw)) # By the stable CST alone, which each woman has from birth, so worked out once for each
        susceptible, cst1, cst3, presymptomatic, asymptomatic, symptomatic, colonized = ss.states('susceptible cst1 cst3 presymptomatic asymptomatic symptomatic colonized')
        cst4 = presymptomatic + asymptomatic + symptomatic
        female = ss.person.female

        # Girls take their stable CST at birth, and everyone at the start (see p_assign())
        assign = ss.Transition(susceptible, {cst1: partial(self.p_assign, cst=1), cst3: partial(self.p_assign, cst=3), presymptomatic: partial(self.p_assign, cst=4), asymptomatic: None}, dur=0, who=female)

        # Women move to a worse CST, with a probability per month raised by their risk factors (see p_worse())
        worse1 = ss.Transition(cst1, cst3, rate=partial(self.p_worse, cst=1))
        worse3 = ss.Transition(cst3, {presymptomatic: lifelong(self.p_symp_now), asymptomatic: None}, rate=partial(self.p_worse, cst=3))
        recover3 = ss.Transition(cst3, cst1, dur=dur_cst3, who=lifelong(self.stable_is, cst=1)) # Only those whose stable CST is 1 leave CST 3 on their own

        # CST 4: symptoms, and clearance to the CST that follows it (see p_clear_to())
        symptoms = presymptomatic >> symptomatic @ dur_presymp
        clear_symp = ss.Transition(symptomatic, {cst1: lifelong(self.p_clear_to, cst=1), cst3: None}, dur=dur_cst4)
        clear_asymp = ss.Transition(asymptomatic, {cst1: lifelong(self.p_clear_to, cst=1), cst3: None}, dur=dur_cst4)

        # Sexual transmission: men from women in CST 4, and women in CST 1 or 3 from men
        sexual = sexual_beta(beta_m2f, beta_m2f*rel_beta_f2m)
        male_infection = ss.Transition(susceptible, colonized, rate=sexual*cst4, who=~female)
        female_infection = ss.Transition(cst1 + cst3, {presymptomatic: lifelong(self.p_symp_now), asymptomatic: None}, rate=sexual*colonized)
        male_clearance = colonized >> susceptible @ dur_male

        # Parameters read by its own code
        values = sc.mergedicts(dict(init_prev=0), pars, kwargs) # Nobody is seeded: women start in their stable CST
        code_pars = dict(
            eff_condom = 0.45, # Read by the sexual networks (see sti.BaseNetwork.net_beta()); https://journals.lww.com/epidem/Fulltext/2007/11000/Condom_Use_and_its_Association_With_Bacterial.9.aspx
            stable_cst_distribution = ss.choice([1, 3, 4], p=[0.10, 0.30, 0.60]), # https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01096-9
            p_symp = sc.objdict(stable_cst1=0.8, stable_cst3=0.7, stable_cst4=0.1), # Probability that CST 4 is symptomatic, by stable CST
            spontaneous_clearance = sc.objdict(stable_cst1=1, stable_cst3=3, stable_cst4=3), # The CST reached when CST 4 clears, by stable CST: https://pmc.ncbi.nlm.nih.gov/articles/PMC9387550/
            p_cst_change = sc.objdict(cst1=0.1, cst3=0.05), # Probability per month of moving from CST 1 to 3, and from 3 to 4
            rr_stable_cst1 = 0.1, # Relative risk of a CST change for those with stable CST 1
            rr_stable_cst4 = 2, # ... stable CST 4 (stable CST 3 is the reference)
            p_douching = ss.bernoulli(p=0.64), # Nigeria (https://www.sciencedirect.com/science/article/abs/pii/S1083318820302400), Ghana (https://pmc.ncbi.nlm.nih.gov/articles/PMC6368746/)
            p_poor_menstrual_hygiene = ss.bernoulli(p=0.55), # https://pmc.ncbi.nlm.nih.gov/articles/PMC9817285/
            n_partners = ss.poisson(mean=0.05), # Number of current partners
            count_partners = False, # Count current partners from the sexual network to determine concurrency
            rr_douching = 1.21, # https://pmc.ncbi.nlm.nih.gov/articles/PMC2574994/
            rr_poor_menstrual_hygiene = 4.1, # From the best fit of the Kenya calibration
            rr_concurrency = 1.28, # https://pmc.ncbi.nlm.nih.gov/articles/PMC5429208/
            p_circumcised = ss.bernoulli(0.4), # Share of men who are circumcised
            rr_uncircumcised = 4.0, # Insertive sex with an uncircumcised penis, from the best fit of the Kenya calibration
            or_ptb = {1: 3, 2: 2, 3: 1}, # Odds ratio of preterm birth by trimester: https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2020.567885/full, https://pmc.ncbi.nlm.nih.gov/articles/PMC8117784/
        )
        given = {key:values.pop(key) for key in code_pars if key in values}

        super().__init__(
            assign=assign, worse1=worse1, worse3=worse3, recover3=recover3, symptoms=symptoms,
            clear_symp=clear_symp, clear_asymp=clear_asymp,
            male_infection=male_infection, female_infection=female_infection, male_clearance=male_clearance,
            new_sexually_acquired_infections=1*female_infection, new_symptomatic=1*symptoms, new_infections_sex=female_infection + male_infection, # Results of the flows
            name=name, **sc.mergedicts(dict(dt='month'), values), # Monthly unless dt is given
        )
        self.define_pars(**code_pars)
        self.pars.update(given)

        # Each person's own: the stable CST and the risk factors, and treatment (see sti.treat_BV)
        self.define_states(
            ss.FloatArr('stable_cst', default=self.pars.stable_cst_distribution, label='Stable CST'), # From birth (read only for women)
            ss.BoolArr('douching', label='Douching'),
            ss.BoolArr('poor_menstrual_hygiene', label='Poor menstrual hygiene'),
            ss.BoolArr('concurrency', default=False, label='Concurrent sexual partners'),
            ss.BoolArr('circumcised', label='Circumcised'),
            ss.BoolArr('has_traits', default=False, label='Risk factors drawn'),
            ss.FloatArr('rel_cst_change', default=1.0, label='Relative risk of a worse CST'),
            ss.FloatArr('ti_treated', label='Time of treatment'),
            ss.BoolArr('on_tx', default=False, label='On treatment'),
        )
        self.set_strata([15, 50], age_bins, sex_keys)
        self.include_care = include_care
        return

    # The states and conditions that are not states
    @property
    def cst4(self):
        """ Women over 15 in CST 4 (BV) """
        return (self.presymptomatic | self.asymptomatic | self.symptomatic) & self.adult_women

    @property
    def infected(self):
        """ Women with BV, and men carrying BV-associated microbes """
        return self.presymptomatic | self.asymptomatic | self.symptomatic | self.colonized

    @property
    def treatable(self):
        """ Responds to treatment """
        return self.infected

    @property
    def adult_women(self):
        return self.sim.people.female & (self.sim.people.age > 15)

    @property
    def stable_cst1(self):
        return self.adult_women & (self.stable_cst == 1)

    @property
    def stable_cst3(self):
        return self.adult_women & (self.stable_cst == 3)

    @property
    def stable_cst4(self):
        return self.adult_women & (self.stable_cst == 4)

    # The values of the declared transitions that depend on each woman's stable CST and risk factors: functions of the
    # agents, as methods (with functools.partial for an argument), so that they read the sim's own copy of the disease
    def stable_is(self, sim, cst):
        """ Whether each woman's stable CST is this one (as who= reads it, for everyone) """
        return self.stable_cst == cst

    def p_assign(self, sim, uids, cst):
        """ The CST each girl takes at birth (and each woman at the start): her stable CST, except that a girl whose stable CST is 4 starts in CST 3; for CST 4, the share symptomatic """
        stable = self.stable_cst[uids]
        adult = sim.people.age[uids] > 15
        if cst == 3:
            return (stable == 3) | ((stable == 4) & ~adult)
        elif cst == 4:
            return ((stable == 4) & adult)*self.p_symp_now(sim, uids)
        return stable == cst

    def p_symp_now(self, sim, uids):
        """ The probability that CST 4 is symptomatic, by stable CST """
        p = self.pars.p_symp
        stable = self.stable_cst[uids]
        return np.select([stable == 1, stable == 3, stable == 4], [p.stable_cst1, p.stable_cst3, p.stable_cst4], 0.0)

    def p_clear_to(self, sim, uids, cst):
        """ Whether CST 4 clears to this CST, by stable CST (see spontaneous_clearance) """
        to = self.pars.spontaneous_clearance
        stable = self.stable_cst[uids]
        dest = np.select([stable == 1, stable == 3, stable == 4], [to.stable_cst1, to.stable_cst3, to.stable_cst4], 3)
        return dest == cst

    def p_worse(self, sim, uids, cst):
        """ The rate of moving from this CST to the next worse one, of women over 15: a probability per month, raised by their risk """
        p = self.pars.p_cst_change[f'cst{cst}']*self.rel_cst_change[uids]
        p = np.where(sim.people.female[uids] & (sim.people.age[uids] > 15), np.minimum(p, 1 - 1e-9), 0.0)
        return ss.probpermonth(p)

    # Each person's own risk factors
    def set_traits(self, uids):
        """ At birth (and for everyone at the start): the risk factors of women, and the circumcision of men (the stable CST is drawn as each person is made) """
        female = uids[self.sim.people.female[uids]]
        male = uids[self.sim.people.male[uids]]
        self.douching[female] = self.pars.p_douching.rvs(female)
        self.poor_menstrual_hygiene[female] = self.pars.p_poor_menstrual_hygiene.rvs(female)
        self.concurrency[female] = (self.pars.n_partners.rvs(female) + 1) > 1
        self.circumcised[male] = self.pars.p_circumcised.rvs(male)
        self.has_traits[uids] = True
        return

    def set_risk(self):
        """ Each woman's risk of a worse CST, and each person's susceptibility to sexual transmission """
        p = self.pars
        risk = self.rel_cst_change
        risk[:] = 1
        risk[self.douching.uids] *= p.rr_douching
        risk[self.poor_menstrual_hygiene.uids] *= p.rr_poor_menstrual_hygiene
        if p.count_partners: # Concurrency from the sexual network
            for net in self.sim.networks.values():
                if isinstance(net, (ss.SexualNet, BaseNetwork)):
                    female = self.sim.people.female.uids
                    self.concurrency[female] = net.partners[female] > 1
        risk[self.concurrency.uids] *= p.rr_concurrency
        risk[self.stable_cst1.uids] *= p.rr_stable_cst1
        risk[self.stable_cst4.uids] *= p.rr_stable_cst4

        # Women in CST 3 with an uncircumcised partner
        uncirc = (self.sim.people.male & ~self.circumcised).uids
        for net in self.sim.networks.values():
            if isinstance(net, (ss.SexualNet, BaseNetwork)):
                partners = ss.uids(net.find_contacts(uncirc)) & self.cst3.uids
                risk[partners] *= p.rr_uncircumcised

        # Susceptibility to sexual transmission: of women who douche, and uncircumcised men
        self.rel_sus[:] = 1
        self.rel_sus[self.douching.uids] *= p.rr_douching
        self.rel_sus[uncirc] *= p.rr_uncircumcised
        return

    def update_pregnancy(self, uids):
        """
        Raise the risk of preterm birth of women who acquire BV while pregnant, by trimester (``or_ptb``)

        ``ss.Pregnancy`` draws each pregnancy's duration at conception from the woman's relative risk of
        preterm birth (``preterm_risk()``), so this redraws it with that risk times the odds ratio.
        """
        if 'pregnancy' not in self.sim.demographics:
            return
        preg = self.sim.demographics.pregnancy
        uids = uids.intersect(preg.P.uids)
        if len(uids):
            trimester = np.digitize(preg.P.time_in[uids].weeks, [13, 27]) + 1
            preg.preterm_risk(uids, rel=[self.pars.or_ptb[t] for t in trimester]) # For this redraw, on top of her own risk
        return

    def clear_infection(self, uids):
        """ Cure people (e.g. treated): women return to the CST that follows CST 4, and men are no longer colonized """
        to_cst1 = uids[self.p_clear_to(self.sim, uids, cst=1)]
        self.move(self.cst1, uids=to_cst1)
        self.move(self.cst3, self.susceptible, uids=uids.remove(to_cst1))
        return

    def init_post(self):
        super().init_post()
        self.set_traits(self.sim.people.auids)
        return

    def step_state(self):
        """ At the start of the step: traits of the newborn, and everyone's risk (the declared transitions are taken in step()) """
        self.set_traits((~self.has_traits).uids)
        self.set_risk()
        return

    def step(self):
        """ The declared transitions; then the preterm risk of the pregnant women who acquired BV """
        super().step()
        self.update_pregnancy(self.worse3.uids | self.female_infection.uids)
        return

    def new_infected(self):
        """ Who entered CST 4 (or, for men, was colonized) in this step """
        new = self.worse3.uids | self.female_infection.uids | self.male_infection.uids
        if self.ti == 0:
            new = new | self.infected.uids
        return self.sim.people.uid.isin(new)

    def init_results(self):
        """ Initialize results """
        super().init_results()
        self.define_sti_results(include_care=self.include_care)
        self.define_results(
            ss.Result('cst1_prevalence', scale=False, label='CST 1 Prevalence', auto_plot=False),
            ss.Result('cst3_prevalence', scale=False, label='CST 3 Prevalence', auto_plot=False),
            ss.Result('cst4_prevalence', scale=False, label='CST 4 Prevalence'),
            ss.Result('symp_prevalence', scale=False, label='Symptomatic prevalence', auto_plot=False),
        )
        return

    def record_states(self):
        super().record_states()
        t = self.tally()
        self.record_sti_stocks(t)
        ti = self.ti
        res = self.results
        women = t.vals(self.sim.people.female) & (np.asarray(self.sim.people.age.values) >= 15)
        n_women = t.count(women)
        for key, state in [('cst1', self.cst1), ('cst3', self.cst3), ('cst4', self.cst4), ('symp', self.symptomatic)]:
            res[f'{key}_prevalence'][ti] = t.div(t.count(t.vals(state) & women), n_women)
        return

    def record_events(self):
        super().record_events()
        self.record_sti_flows()
        return

    def finalize_results(self):
        super().finalize_results()
        self.finalize_sti_results()
        return
