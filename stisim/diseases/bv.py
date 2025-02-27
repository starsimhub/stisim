"""
BV module
"""

import numpy as np
import starsim as ss
import sciris as sc

__all__ = ['BV']


class BV(ss.Disease):

    def __init__(self, pars=None, name='bv', **kwargs):
        super().__init__(name=name)

        self.define_pars(
            unit='month',
            include_care=True,
            p_symp=ss.bernoulli(p=0.5),
            dur_presymp=ss.uniform(ss.dur(1, 'week'), ss.dur(8, 'week')),  # Duration of presymptomatic period
            dur_asymp2clear=ss.uniform(ss.dur(1, 'week'), ss.dur(18, 'week')),  # Duration of asymptomatic infection
            dur_symp2clear=ss.uniform(ss.dur(1, 'week'), ss.dur(18, 'week')),  # Duration of symptoms
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
                n_partners_12m=2,       # OR of BV for each additional partner in the past 12 months - not working yet
                poor_menstrual_hygiene=2,    # OR of BV for poor menstrual hygiene
            ),

            # Care-seeking and clearance
            dur_symp2care=ss.lognorm_ex(ss.dur(4, 'week'), ss.dur(4, 'week')),
            p_symp_care=ss.bernoulli(p=0.3),
            p_clear=ss.bernoulli(p=0.5),
            dur_persist=ss.constant(ss.dur(100, 'year')),
        )
        self.update_pars(pars, **kwargs)

        # States that elevate risk of BV
        self.define_states(
            ss.State('susceptible', default=True, label='Susceptible'),
            ss.State('bv_prone', label='Prone to BV', default=ss.bernoulli(p=0.5)),  # Percentage of women "prone" to BV
            ss.State('infected', label='Infected'),
            ss.State('asymptomatic', label='Asymptomatic'),
            ss.State('symptomatic', label='Symptomatic'),
            ss.FloatArr('rel_sus', default=1.0, label='Relative susceptibility'),
            ss.FloatArr('rel_trans', default=1.0, label='Relative transmissibility'),  # NOT USED
            ss.FloatArr('ti_seeks_care', label='Time of care seeking'),
            ss.BoolArr('seeking_care', label='Seeking care'),
            ss.FloatArr('ti_infected', label='Time of infection'),
            ss.FloatArr('ti_clearance', label='Time of clearance'),
            ss.FloatArr('ti_symptomatic', label='Time of symptoms'),
            ss.FloatArr('dur_inf', label='Duration of infection'),
            ss.BoolArr('douching', label='Douching'),
            ss.FloatArr('n_partners_12m', 0, label='Number of partners in the past 12 months'),
            ss.BoolArr('poor_menstrual_hygiene', label='Poor menstrual hygiene'),
        )

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
            ss.Result('symp_prevalence', scale=False, label="Symptomatic prevalence"),
            ss.Result('incidence', scale=False, label="Incidence"),
            ss.Result('new_infections', dtype=int, label="New infections"),
            ss.Result('new_symptomatic', dtype=int, label="New symptomatic"),
        ]

        if self.pars.include_care:
            sexkeys = ['', 'f', 'm']
            for sk in sexkeys:
                skk = '' if sk == '' else f'_{sk}'
                skl = '' if sk == '' else f' - {sk.upper()}'
                results += [
                    ss.Result('new_care_seekers'+skk, dtype=int, label="New care seekers"+skl),
                    ss.Result('new_false_pos'+skk, dtype=int, label="New false positives"+skl),
                    ss.Result('new_true_pos'+skk, dtype=int, label="New true positives"+skl),
                    ss.Result('new_false_neg'+skk, dtype=int, label="New false negatives"+skl),
                    ss.Result('new_true_neg'+skk, dtype=int, label="New true negatives"+skl),
                    ss.Result('new_treated_success'+skk, dtype=int, label="Successful treatments"+skl),
                    ss.Result('new_treated_failure'+skk, dtype=int, label="Unsuccessful treatments"+skl),
                    ss.Result('new_treated_unnecessary'+skk, dtype=int, label="Unnecessary treatments"+skl),
                    ss.Result('new_treated'+skk, dtype=int, label="Treatments"+skl),
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
        for state in self.states:
            if not state.initialized:
                state.init_vals()
        self.initialized = True

        # Set hygiene states and initial infections
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
        self.ti_clearance[uids] = self.sim.ti

    def step_state(self):
        """ Updates for this timestep """
        ti = self.sim.ti

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

        def cond_prob(num, denom):
            n_num = np.count_nonzero(num & denom)
            n_denom = np.count_nonzero(denom)
            return sc.safedivide(n_num, n_denom)

        self.results['prevalence'][ti] = cond_prob(self.infected, women)
        self.results['symp_prevalence'][ti] = cond_prob(self.symptomatic, women)
        self.results['new_infections'][ti] = np.count_nonzero(self.ti_infected == ti)
        self.results['new_symptomatic'][ti] = np.count_nonzero(self.ti_symptomatic == ti)

        return
