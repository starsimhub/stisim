"""
Bacterial vaginosis (BV): a model of the vaginal microbiome by community state type (CST)
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.diseases.sti import STIResults, sexual_beta, default_age_bins, default_sex_keys

__all__ = ['BV']


class BV(ss.Infection, STIResults):
    """
    Bacterial vaginosis, by the community state type (CST) of the vaginal microbiome

    Each woman has a stable CST (1, 3, or 4), drawn at birth, and is in CST 1, CST 3, or CST 4 (BV) at any time.
    Women over 15 move from CST 1 to CST 3, and from CST 3 to CST 4, with a probability per month
    (`p_cst_change`) raised by behavioral risk factors (douching, poor menstrual hygiene, concurrent partners,
    and, in CST 3, an uncircumcised partner) and by her stable CST. CST 4 is symptomatic for a share
    (`p_symp`, by stable CST) after `dur_presymp`, and clears after `dur_cst4` (from the symptoms, for those
    with them) to the CST of `spontaneous_clearance`; a woman whose stable CST is 1 returns from CST 3 to it
    after `dur_cst3`. BV-associated microbes are sexually transmitted if `beta_m2f` is given: men carry them
    (`colonized`) for `dur_male`.

    States: `susceptible` (men; and girls, at birth, before they take their stable CST), `cst1`, `cst3`, the
    states of CST 4 (`presymptomatic`, `asymptomatic`, `symptomatic`), and `colonized` (men); `infected` is
    CST 4 or colonized, and `cst4` the women over 15 in CST 4. The women whose symptoms started in the step
    are `symptoms` (e.g. for `sti.treat_BV`). `set_state()` puts women in a CST.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'bv'``)
        age_bins (list): bin edges for age-stratified results; ``None`` skips them
        sex_keys (dict): sex strata for results; ``None`` skips sex-stratified results
        include_care (bool): whether to define the results of care
        kwargs (dict): parameter values, or `label`
    """
    def __init__(self, pars=None, name='bv', age_bins=default_age_bins, sex_keys=default_sex_keys, include_care=True, **kwargs):
        super().__init__(name=name, label=kwargs.pop('label', None))
        self.define_pars(
            dt = 'month',
            init_prev = 0, # Nobody is seeded: women start in their stable CST

            # Natural history
            dur_cst3 = ss.uniform(ss.weeks(1), ss.weeks(18)), # Time in CST 3 before returning to a stable CST 1
            dur_cst4 = ss.uniform(ss.weeks(5), ss.weeks(50)), # Time in CST 4 before it clears, from the symptoms for those with them: https://www.sciencedirect.com/science/article/abs/pii/S0002937803010421
            dur_presymp = ss.uniform(ss.weeks(1), ss.weeks(2)), # Time from CST 4 to symptoms
            stable_cst_distribution = ss.choice([1, 3, 4], p=[0.10, 0.30, 0.60]), # https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01096-9
            p_symp = sc.objdict(stable_cst1=0.8, stable_cst3=0.7, stable_cst4=0.1), # Probability that CST 4 is symptomatic, by stable CST
            spontaneous_clearance = sc.objdict(stable_cst1=1, stable_cst3=3, stable_cst4=3), # The CST reached when CST 4 clears, by stable CST: https://pmc.ncbi.nlm.nih.gov/articles/PMC9387550/
            p_cst_change = sc.objdict(cst1=0.1, cst3=0.05), # Probability per month of moving from CST 1 to 3, and from 3 to 4

            # Risk factors
            rr_stable_cst1 = 0.1, # Relative risk of a CST change for those with stable CST 1
            rr_stable_cst3 = 1, # ... stable CST 3 (the reference)
            rr_stable_cst4 = 2, # ... stable CST 4
            p_douching = ss.bernoulli(p=0.64), # Nigeria (https://www.sciencedirect.com/science/article/abs/pii/S1083318820302400), Ghana (https://pmc.ncbi.nlm.nih.gov/articles/PMC6368746/)
            p_poor_menstrual_hygiene = ss.bernoulli(p=0.55), # https://pmc.ncbi.nlm.nih.gov/articles/PMC9817285/
            n_partners = ss.poisson(lam=0.05), # Number of current partners
            count_partners = False, # Count current partners from the sexual network to determine concurrency
            rr_douching = 1.21, # https://pmc.ncbi.nlm.nih.gov/articles/PMC2574994/
            rr_poor_menstrual_hygiene = 4.1, # From the best fit of the Kenya calibration
            rr_concurrency = 1.28, # https://pmc.ncbi.nlm.nih.gov/articles/PMC5429208/
            p_circumcised = ss.bernoulli(0.4), # Share of men who are circumcised
            rr_uncircumcised = 4.0, # Insertive sex with an uncircumcised penis, from the best fit of the Kenya calibration
            or_ptb = {1: 3, 2: 2, 3: 1}, # Odds ratio of preterm birth by trimester: https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2020.567885/full, https://pmc.ncbi.nlm.nih.gov/articles/PMC8117784/

            # Sexual transmission of BV-associated microbes, if beta_m2f is given
            beta = 0.0, # On each network, from the ones below (see validate_beta())
            beta_m2f = 0.0, # Probability of transmission per act from a man to a woman
            rel_beta_f2m = 0.5, # Transmission from a woman to a man, relative to from a man to a woman
            dur_male = ss.uniform(ss.weeks(1), ss.weeks(18)), # Time that a man carries BV-associated microbes
            eff_condom = 0.45, # Read by the sexual networks (see sti.BaseNetwork.net_beta()); https://journals.lww.com/epidem/Fulltext/2007/11000/Condom_Use_and_its_Association_With_Bacterial.9.aspx
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.BoolState('cst1', label='CST 1'),
            ss.BoolState('cst3', label='CST 3'),
            ss.BoolState('presymptomatic', label='Presymptomatic'),
            ss.BoolState('asymptomatic', label='Asymptomatic'),
            ss.BoolState('symptomatic', label='Symptomatic'),
            ss.BoolState('colonized', label='Colonized'),
            ss.BoolArr('symptoms', label='Symptoms started this step'),
            ss.TimeArr('ti_next', label='Time of the next scheduled move: symptoms, clearance, the return to CST 1, or the end of colonization'),

            # Each person's own: the stable CST and the risk factors, and treatment (see sti.treat_BV)
            ss.FloatArr('stable_cst', label='Stable CST'),
            ss.BoolArr('douching', label='Douching'),
            ss.BoolArr('poor_menstrual_hygiene', label='Poor menstrual hygiene'),
            ss.BoolArr('concurrency', default=False, label='Concurrent sexual partners'),
            ss.BoolArr('circumcised', label='Circumcised'),
            ss.FloatArr('rel_cst_change', default=1.0, label='Relative risk of a worse CST'),
            ss.FloatArr('ti_treated', label='Time of treatment'),
            ss.BoolArr('on_tx', default=False, label='On treatment'),
        )
        self.u_worse = ss.random(name='u_worse') # Whether each woman moves to a worse CST in a step
        self.u_arrived = ss.random(name='u_arrived') # ... and in the half step after she arrives in CST 1 or 3 during one
        self.u_symp = ss.random(name='u_symp') # Whether CST 4 is symptomatic
        self.set_strata([15, 50], age_bins, sex_keys)
        self.include_care = include_care
        return

    @property
    def cst4(self):
        """ Women over 15 in CST 4 (BV) """
        return (self.presymptomatic | self.asymptomatic | self.symptomatic) & self.adult_women

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

    def by_stable(self, values, uids, other):
        """ A value for each woman by her stable CST (`values` has `stable_cst1`, `stable_cst3`, and `stable_cst4`), or `other` for none """
        stable = self.stable_cst[uids]
        return np.select([stable == 1, stable == 3, stable == 4], [values.stable_cst1, values.stable_cst3, values.stable_cst4], other)

    def set_state(self, uids, state, start=None):
        """
        Put some women in a CST (`cst1`, `cst3`, or, for CST 4, `presymptomatic` or `asymptomatic`), out of the others, with its scheduled
        move: the symptoms, the clearance of CST 4 without them, or, for those whose stable CST is 1, the return from CST 3, timed from `start` (by default now)
        """
        if not len(uids):
            return
        p = self.pars
        for s in [self.cst1, self.cst3, self.presymptomatic, self.asymptomatic, self.symptomatic, self.susceptible]:
            s[uids] = False
        state[uids] = True
        cst4 = state is self.presymptomatic or state is self.asymptomatic
        self.infected[uids] = cst4
        self.ti_next[uids] = np.nan
        start = self.ti if start is None else start
        if cst4:
            self.ti_infected[uids] = self.ti
            self.ti_next[uids] = start + (p.dur_presymp if state is self.presymptomatic else p.dur_cst4).rvs(uids)
        elif state is self.cst3:
            back = uids[self.stable_cst[uids] == 1] # Only those whose stable CST is 1 leave CST 3 on their own
            self.ti_next[back] = start + p.dur_cst3.rvs(back)
        return

    def to_cst4(self, uids, start=None):
        """ Women who enter CST 4: symptoms to come for a share, by stable CST """
        symp = self.u_symp.rvs(uids) < self.by_stable(self.pars.p_symp, uids, 0.0)
        self.set_state(uids[symp], self.presymptomatic, start)
        self.set_state(uids[~symp], self.asymptomatic, start)
        self.update_pregnancy(uids)
        return

    def clear(self, uids):
        """ Women whose CST 4 clears: to the CST that follows it, by stable CST (see spontaneous_clearance) """
        to_cst1 = self.by_stable(self.pars.spontaneous_clearance, uids, 3) == 1
        self.set_state(uids[to_cst1], self.cst1)
        self.set_state(uids[~to_cst1], self.cst3)
        return

    def assign(self, uids):
        """ Girls at birth (and women at the start) take their stable CST, except that a girl whose stable CST is 4 starts in CST 3 """
        stable = self.stable_cst[uids]
        adult = self.sim.people.age[uids] > 15
        self.set_state(uids[stable == 1], self.cst1)
        self.set_state(uids[(stable == 3) | ((stable == 4) & ~adult)], self.cst3)
        self.to_cst4(uids[(stable == 4) & adult])
        return

    def uncolonize(self, uids):
        """ Men who no longer carry BV-associated microbes """
        self.colonized[uids] = False
        self.infected[uids] = False
        self.susceptible[uids] = True
        self.ti_next[uids] = np.nan
        return

    def clear_infection(self, uids):
        """ Cure people (e.g. treated): women return to the CST that follows CST 4, and men are no longer colonized """
        women = self.sim.people.female[uids]
        self.clear(uids[women])
        self.uncolonize(uids[~women])
        return

    # Each person's own risk factors
    def set_traits(self, uids):
        """ At birth (and for everyone at the start): the stable CST and risk factors of women, and the circumcision of men """
        female = uids[self.sim.people.female[uids]]
        male = uids[self.sim.people.male[uids]]
        self.stable_cst[female] = self.pars.stable_cst_distribution.rvs(female)
        self.douching[female] = self.pars.p_douching.rvs(female)
        self.poor_menstrual_hygiene[female] = self.pars.p_poor_menstrual_hygiene.rvs(female)
        self.concurrency[female] = (self.pars.n_partners.rvs(female) + 1) > 1
        self.circumcised[male] = self.pars.p_circumcised.rvs(male)
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
                if isinstance(net, ss.SexualNet):
                    female = self.sim.people.female.uids
                    self.concurrency[female] = net.partners[female] > 1
        risk[self.concurrency.uids] *= p.rr_concurrency
        risk[self.stable_cst1.uids] *= p.rr_stable_cst1
        risk[self.stable_cst4.uids] *= p.rr_stable_cst4

        # Women in CST 3 with an uncircumcised partner
        uncirc = (self.sim.people.male & ~self.circumcised).uids
        for net in self.sim.networks.values():
            if isinstance(net, ss.SexualNet):
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
        if 'pregnancy' not in self.sim.demographics or not self.sim.initialized:
            return
        preg = self.sim.demographics.pregnancy
        uids = uids.intersect(preg.P.uids)
        if len(uids):
            trimester = np.digitize(preg.P.time[uids].weeks, [13, 27]) + 1
            preg.preterm_risk(uids, rel=[self.pars.or_ptb[t] for t in trimester]) # For this redraw, on top of her own risk
        return

    def init_post(self):
        """ Everyone's traits, and each woman's CST """
        self.set_traits(self.sim.people.auids)
        super().init_post()
        self.assign(self.sim.people.female.uids)
        return

    def step_state(self):
        """ Traits of the newborn and everyone's risk; then the worse CSTs, at a probability per month, the moves due, and half a step of the worse CSTs for those who arrived in CST 1 or 3 """
        ppl = self.sim.people
        newborns = (ppl.female | ppl.male) & np.isnan(self.stable_cst) & ~self.circumcised
        newborns = newborns.uids[ppl.age[newborns.uids] < self.t.dt.years]
        self.set_traits(newborns)
        self.assign(newborns[ppl.female[newborns]])
        self.set_risk()
        self.symptoms[:] = False

        # Worse CSTs, for women over 15, drawn for those in CST 1 or 3 at the start of the step
        months = self.t.dt/ss.months(1)
        adult = self.adult_women
        def worse(key, uids, u, steps): # The women of `uids` who move to a worse CST from `key` in so many steps
            p = np.minimum(self.pars.p_cst_change[key]*self.rel_cst_change[uids], 1)
            return uids[u.rvs(uids) < 1 - (1 - p)**(months*steps)]
        to3, to4 = worse('cst1', (self.cst1 & adult).uids, self.u_worse, 1), worse('cst3', (self.cst3 & adult).uids, self.u_worse, 1)
        self.set_state(to3, self.cst3)
        self.to_cst4(to4)

        # The moves due: symptoms (then clearance, from them), clearance, the return to CST 1, and the end of colonization
        due = lambda state: (state & self.ti_next.due()).uids
        uids = due(self.presymptomatic)
        self.presymptomatic[uids] = False
        self.symptomatic[uids] = True
        self.symptoms[uids] = True
        self.ti_next[uids] += self.pars.dur_cst4.rvs(uids) # From when the symptoms started
        cleared, back = due(self.symptomatic) | due(self.asymptomatic), due(self.cst3)
        self.clear(cleared)
        self.set_state(back, self.cst1)
        self.uncolonize(due(self.colonized))

        # Those who arrived in CST 1 or 3 during the step face its worse CST for half a step (as for an arrival during a step, on average), once
        arrived = (to3 | cleared | back) & adult.uids
        to3, to4 = worse('cst1', arrived[self.cst1[arrived]], self.u_arrived, 0.5), worse('cst3', arrived[self.cst3[arrived]], self.u_arrived, 0.5)
        self.set_state(to3, self.cst3)
        self.to_cst4(to4)
        return

    def validate_beta(self):
        """ The betas of the sexual networks, from the parameters (so a fitted beta_m2f is read as it is now) """
        p = self.pars
        self.pars.beta = sexual_beta(p.beta_m2f, p.beta_m2f*p.rel_beta_f2m)
        return super().validate_beta()

    def transmit(self):
        """ Sexual transmission: men from women in CST 4, and women in CST 1 or 3 from men who carry BV-associated microbes """
        if not self.pars.beta_m2f:
            return self.finalize_infections([], [], [])
        betamap = self.validate_beta()
        source = self.presymptomatic.raw | self.asymptomatic.raw | self.symptomatic.raw | self.colonized.raw
        target = (self.susceptible.raw & ~self.sim.people.female.raw) | self.cst1.raw | self.cst3.raw
        rel_trans = self.rel_trans.asnew(source*self.rel_trans.raw, copy=False)
        rel_sus = self.rel_sus.asnew(target*self.rel_sus.raw, copy=False)
        out = [self.transmit_route(i, route, betamap[ss.utils.standardize_netkey(nkey)], rel_trans, rel_sus) for i,(nkey,route) in enumerate(self.sim.networks.items())]
        return self.finalize_infections(*zip(*out)) if out else self.finalize_infections([], [], [])

    def set_prognoses(self, uids, sources=None):
        """ Those infected sexually: women enter CST 4, and men are colonized """
        super().set_prognoses(uids, sources)
        women = self.sim.people.female[uids]
        start = self.ti + 1 if self.sim.initialized else self.ti # A course runs from when the infection happens: at the end of the step, in transmission (as a declared disease dates it), or as ss.Infection.init_post() dates those seeded at the start
        self.to_cst4(uids[women], start)
        men = uids[~women]
        self.colonized[men] = True
        self.infected[men] = True
        self.ti_infected[men] = self.ti
        self.ti_next[men] = start + self.pars.dur_male.rvs(men)
        self.results.new_sexually_acquired_infections[self.ti] += np.count_nonzero(women)
        self.results.new_infections_sex[self.ti] += len(uids)
        return

    def new_infected(self):
        """ Who entered CST 4 (or, for men, was colonized) in this step, as a boolean array """
        return self.ti_infected == self.ti

    def init_results(self):
        """ Initialize results """
        super().init_results()
        self.define_sti_results(include_care=self.include_care)
        self.define_results(
            ss.Result('cst1_prevalence', scale=False, label='CST 1 Prevalence', auto_plot=False),
            ss.Result('cst3_prevalence', scale=False, label='CST 3 Prevalence', auto_plot=False),
            ss.Result('cst4_prevalence', scale=False, label='CST 4 Prevalence'),
            ss.Result('symp_prevalence', scale=False, label='Symptomatic prevalence', auto_plot=False),
            ss.Result('new_symptomatic', dtype=int, label='New symptomatic', auto_plot=False),
            ss.Result('new_sexually_acquired_infections', dtype=int, label='New sexually acquired infections', auto_plot=False),
        )
        return

    def update_results(self):
        super().update_results()
        t = self.tally()
        self.update_sti_results(t)
        ti = self.ti
        res = self.results
        women = t.vals(self.sim.people.female) & (np.asarray(self.sim.people.age.values) >= 15)
        n_women = t.count(women)
        for key, state in [('cst1', self.cst1), ('cst3', self.cst3), ('cst4', self.cst4), ('symp', self.symptomatic)]:
            res[f'{key}_prevalence'][ti] = t.div(t.count(t.vals(state) & women), n_women)
        res.new_symptomatic[ti] = t.count(t.vals(self.symptoms))
        return
