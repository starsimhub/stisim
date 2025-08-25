"""
BV modules
Includes:
    - a simple BV model for generating background prevalence of vaginal discharge
    - a detailed model of the vaginal microbiome including community state types (CSTs)
"""

import numpy as np
import starsim as ss
import sciris as sc
import stisim
from stisim.diseases.sti import BaseSTI, BaseSTIPars
import stisim.utils as ut

__all__ = ['SimpleBV', 'BV']


class SimpleBV(ss.Disease):

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
                n_partners_12m=2,       # OR of BV for each additional partner in the past 12 months - not working yet
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
            ss.FloatArr('n_partners_12m', 0, label='Number of partners in the past 12 months'),
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
        self.results['new_infections'][ti] = ut.count(self.ti_infected == ti)
        self.results['new_symptomatic'][ti] = ut.count(self.ti_symptomatic == ti)

        return


class BVPars(BaseSTIPars):
    """  Parameters for the BV model """
    def __init__(self, **kwargs):
        super().__init__()

        # Transmission
        self.sexual_transmission = False  # Flag for model to know if we are doing sexual transmission of BV microbes
        self.eff_condom = 0.45  # Condom effectiveness in reducing BV risk https://journals.lww.com/epidem/Fulltext/2007/11000/Condom_Use_and_its_Association_With_Bacterial.9.aspx
        self.p_symp = sc.objdict(  # Probability of symptomatic BV for women given stable CST state
            stable_cst1=ss.bernoulli(p=0.8),
            stable_cst3=ss.bernoulli(p=0.7),
            stable_cst4=ss.bernoulli(p=0.1),
        )
        self.init_prev = ss.bernoulli(
            p=0.23
        )  # Initial prevalence of BV (https://www.who.int/news-room/fact-sheets/detail/bacterial-vaginosis)
        self.stable_cst_distribution = ss.choice(
            [1, 3, 4], p=[0.10, 0.30, 0.60]
        )  # Distribution of stable CST states https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01096-9
        self.dur2clear = sc.objdict(  # Time until transitioning back to stable CST state
            cst3=ss.uniform(ss.weeks(1), ss.weeks(18)),  #
            cst4=ss.uniform(
                ss.weeks(5), ss.weeks(50)
            ),  # Natural clearance https://www.sciencedirect.com/science/article/abs/pii/S0002937803010421?via%3Dihub
            male=ss.uniform(ss.weeks(1), ss.weeks(18)),
        )
        self.spontaneous_clearance = sc.objdict(  # Where to transition to after spontaneous clearance
            stable_cst1=1,  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9387550/
            stable_cst3=3,  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9387550/
            stable_cst4=3,  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9387550/
        )
        self.dur_presymp = ss.uniform(ss.weeks(1), ss.weeks(2))
        self.p_spontaneous = ss.bernoulli(
            p=0.1
        )  # Placeholder for probability of spontaneous transition to worse CST state (overwritten)
        self.rr_stable_cst1=0.1  # Relative risk of CST transition for those in stable CST 1
        self.rr_stable_cst3 = 1  # Relative risk of CST transition for those in stable CST 3 (reference)
        self.rr_stable_cst4 = 2  # Relative risk of CST transition for those in stable CST 4
        self.p_cst_change = sc.objdict(
            cst1=0.1,  # Probability of transition from CST 1 to CST 3
            cst3=0.05,  # Probability of transition from CST 3 to CST 4
        )
        # Share of population douching
        #   - Nigeria (https://www.sciencedirect.com/science/article/abs/pii/S1083318820302400)
        #   - Ghana (https://pmc.ncbi.nlm.nih.gov/articles/PMC6368746/)
        # adjust OR for age (OR, 0.2; 95% CI = 0.063-0.603)
        self.p_douching=ss.bernoulli(p=0.64)

        # Share of population with poor menstrual hygiene https://pmc.ncbi.nlm.nih.gov/articles/PMC9817285/
        # adjust for ubanicity (OR, 0.33 CI = 0.25, 0.43) and SES (OR 0.46, CI = (0.30, 0.70))
        self.p_poor_menstrual_hygiene = ss.bernoulli(p=0.55)

        self.n_partners = ss.poisson(lam=0.05)  # Number of current partners
        self.count_partners = False  # Count current partners from sexual network to determine concurrency

        # Relative risks
        self.rr_douching = 1.21  # Douching CI = 1.08, 1.38 https://pmc.ncbi.nlm.nih.gov/articles/PMC2574994/
        self.rr_poor_menstrual_hygiene = 4.1  # Poor menstrual hygiene - pulled from best GOF Kenya Calibration
        self.rr_concurrency = 1.28  # Concurrent partners - https://pmc.ncbi.nlm.nih.gov/articles/PMC5429208/
        self.p_circumcised = ss.bernoulli(0.4)  # Proportion of men who are circumcised
        self.rr_uncircumcised = 4.0  # RR insertive sex with uncircumcised penis, from best GOF Kenya Calibration
        self.or_ptb = {  # PTB https://www.frontiersin.org/journals/public-health/articles/10.3389/fpubh.2020.567885/full
            1: 3,  # First trimester study on CST IV and sPTB: https://pmc.ncbi.nlm.nih.gov/articles/PMC8117784/
            2: 2,
            3: 1,
        }
        self.update(kwargs)
        return


class BV(BaseSTI):

    def __init__(self, pars=None, name="bv", **kwargs):
        super().__init__(name=name)

        # Handle parameters
        default_pars = BVPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        # Define states
        self.define_states(
            ss.BoolState("asymptomatic", default=False, label="Asymptomatic"),
            ss.BoolState("symptomatic", default=False, label="Symptomatic"),
            ss.FloatArr("ti_cst_change", label="Time of CST change"),
            ss.FloatArr("ti_symptomatic", label="Time of clearance"),
            ss.FloatArr("ti_clearance", label="Time of clearance"),
            ss.FloatArr("ti_treated", label="Time of treatment"),
            ss.BoolArr("douching", label="Douching"),
            ss.BoolArr("concurrency", default=False, label="Concurrent sexual partners"),
            ss.BoolArr("poor_menstrual_hygiene", label="Poor menstrual hygiene"),
            ss.BoolArr("condom_use", label="Consistent Condom use"),
            ss.BoolArr("circumcised", label="Circumcised"),
            ss.FloatArr("cst", default=np.nan, label="Current CST"),
            ss.FloatArr("stable_cst", default=np.nan, label="Stable CST"),
            ss.BoolArr("on_tx", default=False, label="On treatment"),
        )
        return

    def init_results(self):
        """Initialize results"""
        super().init_results()
        results = [
            ss.Result("cst1_prevalence", scale=False, label="CST 1 Prevalence"),
            ss.Result("cst3_prevalence", scale=False, label="CST 3 Prevalence"),
            ss.Result("cst4_prevalence", scale=False, label="CST 4 Prevalence"),
            ss.Result("symp_prevalence", scale=False, label="Symptomatic prevalence"),
            ss.Result(
                "new_sexually_acquired_infections",
                dtype=int,
                label="New sexually acquired infections",
            ),
            ss.Result("new_symptomatic", dtype=int, label="New symptomatic"),
            ss.Result("new_MTZ_doses", dtype=int, label="New MTZ doses administered"),
        ]
        self.define_results(*results)
        return

    def _get_fuids(self, upper_age=None):
        """Get uids of females younger than upper_age"""
        people = self.sim.people
        if upper_age is None:
            upper_age = 1000
        within_age = people.age < upper_age
        return (within_age & people.female).uids

    def _get_muids(self, upper_age=None):
        """Get uids of males younger than upper_age"""
        people = self.sim.people
        if upper_age is None:
            upper_age = 1000
        within_age = people.age < upper_age
        return (within_age & people.male).uids

    def bv_sus(self):
        return self.sim.people.female & (self.sim.people.age > 15) & ~self.cst4

    @property
    def cst1(self):
        return self.sim.people.female & (self.sim.people.age > 15) & (self.cst == 1)

    @property
    def cst3(self):
        return self.sim.people.female & (self.sim.people.age > 15) & (self.cst == 3)

    @property
    def cst4(self):
        return self.sim.people.female & (self.sim.people.age > 15) & (self.cst == 4)

    @property
    def stable_cst1(self):
        return (
            self.sim.people.female & (self.sim.people.age > 15) & (self.stable_cst == 1)
        )

    @property
    def stable_cst3(self):
        return (
            self.sim.people.female & (self.sim.people.age > 15) & (self.stable_cst == 3)
        )

    @property
    def stable_cst4(self):
        return (
            self.sim.people.female & (self.sim.people.age > 15) & (self.stable_cst == 4)
        )

    def set_cst(self, upper_age=None):
        f_uids = self._get_fuids(upper_age=upper_age)
        cst = self.pars.stable_cst_distribution.rvs(f_uids)
        self.cst[f_uids] = cst
        self.stable_cst[f_uids] = cst
        self.set_prognoses(f_uids[cst == 4], 4)
        return

    def _get_uids(self, upper_age=None):
        """ Get uids of females younger than upper_age """
        people = self.sim.people
        if upper_age is None: upper_age = 1000
        within_age = people.age < upper_age
        return (within_age & people.female).uids

    def set_hygiene_states(self, upper_age=None):
        """Set vaginal hygiene states"""
        f_uids = self._get_fuids(upper_age=upper_age)
        self.douching[f_uids] = self.pars.p_douching.rvs(f_uids)
        self.poor_menstrual_hygiene[f_uids] = self.pars.p_poor_menstrual_hygiene.rvs(
            f_uids
        )
        self.concurrency[f_uids] = (self.pars.n_partners.rvs(f_uids) + 1) > 1
        return

    def set_circumcision(self, upper_age=None):
        """Set circumcision status"""
        m_uids = self._get_muids(upper_age=upper_age)
        self.circumcised[m_uids] = self.pars.p_circumcised.rvs(m_uids)
        return

    def compute_circumcision_impact(self, spontaneous=True):
        """
        Compute the relative impact of circumcision on susceptibility for women.
        This is used for spontaneous occurence only. If sexual transmission in model,
        then we adjust man's susceptibility.
        """
        if spontaneous:
            if not self.pars.sexual_transmission:
                cst3_uids = self.cst3.uids
                for _, net in self.sim.networks.items():
                    if isinstance(net, ss.SexualNetwork):  # Skip networks that are not sexual

                        # Find partners of uncircumcised men
                        uncirc = (~self.circumcised).uids
                        f_partners = ss.uids(net.find_contacts(uncirc))

                        # Find women who are CST3 and have an uncircumcised partner
                        cst3_uncirc = f_partners & cst3_uids
                        self.rel_sus[cst3_uncirc] *= self.pars.rr_uncircumcised

                        if 'condom' in net.edges:
                            print('WARNING: condom adjustment not implemented without sexual transmission.')

        else:  # Adjust susceptibility to BV microbes for men who are uncircumcised
            m_uids = self._get_muids()
            uncircumcised_uids = m_uids[~self.circumcised[m_uids]]
            self.rel_sus[uncircumcised_uids] *= self.pars.rr_uncircumcised
        return

    def set_rel_sus(self, spontaneous=True):
        # douching impacts BV risk via spontaneous AND sexual transmission
        self.rel_sus[self.douching.uids] *= self.pars.rr_douching

        # Update circumcision based on whether sexual transmission is occurring
        self.compute_circumcision_impact(spontaneous=spontaneous)

        # update for spontaneous transitions
        if spontaneous:

            self.rel_sus[
                self.poor_menstrual_hygiene.uids
            ] *= self.pars.rr_poor_menstrual_hygiene
            if self.pars.sexual_transmission is False:
                if self.pars.count_partners:
                    # Use sexual network to determine concurrency
                    f_uids = self._get_fuids()
                    for net in self.sim.networks:
                        if isinstance(net, ss.SexualNetwork):
                            self.concurrency[f_uids] = net.partners[f_uids] > 1
                self.rel_sus[self.concurrency.uids] *= self.pars.rr_concurrency
            self.rel_sus[self.stable_cst1.uids] *= self.pars.rr_stable_cst1
            self.rel_sus[self.stable_cst4.uids] *= self.pars.rr_stable_cst4

        return

    def init_post(self):
        """Initialize with sim properties"""
        ss.Module.init_post(self)  # Skip the disease init_post() since we create infections in a different way
        for state in self.state_list:
            if not state.initialized:
                state.init_vals()
        self.initialized = True

        # Set cst, hygiene states and circumcision
        self.set_cst()
        self.set_hygiene_states()
        self.set_circumcision()
        return

    def spontaneous(self, uids, cst="cst1"):
        """
        Determine the probability of transitioning to worse CST states
        for agents in CST 1, we determine probability of transitioning to CST 3 (assuming no transition directly to CST 4)
        for agents in CST 3, we determine probability of transitioning to CST 4
        agents in CST 4 have no further transitions

        """
        p_cst_change = np.full_like(uids, fill_value=0, dtype=float)
        p_cst_change = self.pars.p_cst_change[cst] * self.rel_sus[uids]
        return p_cst_change

    def change_cst(self):
        # Transition CST states
        for cst, new_cst in zip(["cst1", "cst3"], [3, 4]):
            cst_uids = getattr(self, cst).uids
            p_cst = self.spontaneous(cst_uids, cst)
            self.pars.p_spontaneous.set(p_cst)
            cst_cases = self.pars.p_spontaneous.filter(cst_uids)
            self.set_prognoses(cst_cases, new_cst)
        return

    def step(self):
        upper_age = self.t.dt.years # TODO this looks wrong!
        self.set_cst(upper_age=upper_age)
        self.set_hygiene_states(upper_age=upper_age)
        self.set_circumcision(upper_age=upper_age)
        self.set_rel_sus(spontaneous=True)

        # First, spontaneous transitions
        self.change_cst()

        if self.pars.sexual_transmission:
            # Now do transmission-based infection
            # Reset susceptibility
            self.rel_sus[:] = 1
            # Update susceptibility for sexual transmission (circumcision, douching)
            self.set_rel_sus(spontaneous=False)
            new_cases = self.sexual_transmission()
            # pull out male uids
            m_uids = new_cases & self.sim.people.male
            self.set_male_prognoses(m_uids)
            f_uids = new_cases & self.sim.people.female
            self.set_prognoses(f_uids, new_cst=4)
            self.results.new_sexually_acquired_infections[self.ti] += len(f_uids)

    def sexual_transmission(self):
        """Determine who gets infected on this timestep via transmission on the network"""
        new_cases, _, _ = super().infect()
        return new_cases

    def clear_infection(self, uids):
        """Clear infection"""
        # First we need to determine where they are going, which depends upon stable CST state
        stable_cst1 = uids.intersect(self.stable_cst1)
        stable_cst3 = uids.intersect(self.stable_cst3)
        stable_cst4 = uids.intersect(self.stable_cst4)
        self.cst[stable_cst1] = self.pars.spontaneous_clearance.stable_cst1
        self.cst[stable_cst3] = self.pars.spontaneous_clearance.stable_cst3
        self.cst[stable_cst4] = self.pars.spontaneous_clearance.stable_cst4
        self.infected[uids] = False
        self.susceptible[uids] = True
        self.ti_clearance[uids] = self.ti
        self.symptomatic[uids] = False
        self.asymptomatic[uids] = False

    def step_state(self):
        """Updates for this timestep"""
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

        # Return to stable CST states
        new_cleared = (self.ti_clearance <= ti).uids
        self.clear_infection(new_cleared)

        return

    def wipe_dates(self, uids):
        """Clear all previous dates"""
        self.ti_cst_change[uids] = np.nan
        self.ti_clearance[uids] = np.nan
        self.ti_symptomatic[uids] = np.nan
        return

    def set_symptoms(self, uids):
        p = self.pars
        self.asymptomatic[uids] = True
        stable_cst1 = uids.intersect(self.stable_cst1)
        stable_cst3 = uids.intersect(self.stable_cst3)
        stable_cst4 = uids.intersect(self.stable_cst4)
        for cst_uids, cst in zip(
            [stable_cst1, stable_cst3, stable_cst4],
            ["stable_cst1", "stable_cst3", "stable_cst4"],
        ):
            if len(cst_uids):
                symp, asymp = p.p_symp[cst].split(cst_uids)
                dur_presymp = self.pars.dur_presymp.rvs(symp)
                self.ti_symptomatic[symp] = self.ti_cst_change[symp] + dur_presymp
                dur = self.pars.dur2clear["cst4"].rvs(symp)
                self.ti_clearance[symp] = dur + dur_presymp + self.ti_cst_change[symp]
                dur = self.pars.dur2clear["cst4"].rvs(asymp)
                self.ti_clearance[asymp] = dur + self.ti_cst_change[asymp]
                # self.results.new_symptomatic[self.ti] += len(symp)

        return

    def set_duration(self, uids, cst=3):
        dur = self.pars.dur2clear[f"cst{cst}"].rvs(uids)
        self.ti_clearance[uids] = dur + self.ti_cst_change[uids]
        return

    def update_pregnancy(self, uids, cleared=False):
        pregnancy = self.sim.demographics.pregnancy
        if cleared:  # Update PTB risk for those who have cleared infection
            pregnancy.rel_sus_ptb[uids] = 1
            pregnant = uids.intersect(pregnancy.pregnant.uids)
            pregnancy.set_prognoses(pregnant, bv_update=True)

        else:  # Update PTB risk for those who are newly infected
            bv_pregnant = uids.intersect(pregnancy.pregnant.uids)
            bv_not_pregnant = np.setdiff1d(uids, bv_pregnant)
            pregnancy.rel_sus_ptb[bv_not_pregnant] = self.pars.or_ptb[1]
            # calculate trimester of pregnancy
            if len(bv_pregnant):
                trimester = pregnancy.trimester[bv_pregnant]
                trimester_or_ptb = [self.pars.or_ptb[i] for i in trimester]
                pregnancy.rel_sus_ptb[bv_pregnant] = trimester_or_ptb
                pregnancy.set_prognoses(bv_pregnant, bv_update=True)

    def set_male_prognoses(self, uids):
        """
        Set initial prognoses for newly infected males
        """
        self.ti_infected[uids] = self.ti
        self.infected[uids] = True
        self.susceptible[uids] = False
        self.ti_clearance[uids] = self.ti + self.pars.dur2clear["male"].rvs(uids)
        return

    def set_prognoses(self, uids, new_cst=3):
        """
        Set initial prognoses for newly infected females
        """
        # Check that all the UIDS are female
        if ut.count(uids & self.sim.people.male) > 0:
            errormsg = 'set_prognoses should only be used for female agents'
            raise ValueError(errormsg)

        self.wipe_dates(uids)  # Clear prior dates
        self.cst[uids] = new_cst
        self.ti_cst_change[uids] = self.ti
        if new_cst == 4:
            self.ti_infected[uids] = self.ti
            self.infected[uids] = True
            self.susceptible[uids] = False
            self.set_symptoms(uids)  # Set symptoms and duration
            if "pregnancy" in self.sim.demographics:
                self.update_pregnancy(uids)
        else:
            self.set_duration(uids, new_cst)  # Set duration of CST transition

        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        women = (self.sim.people.age >= 15) & self.sim.people.female
        men = (self.sim.people.age >= 15) & self.sim.people.male
        self.results["cst1_prevalence"][ti] = ut.cond_prob(self.cst1, women)
        self.results["cst3_prevalence"][ti] = ut.cond_prob(self.cst3, women)
        self.results["cst4_prevalence"][ti] = ut.cond_prob(self.cst4, women)
        self.results["symp_prevalence"][ti] = ut.cond_prob(self.symptomatic, women)
        return
