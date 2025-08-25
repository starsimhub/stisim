"""
Template class for SEIS-type STIs
Used for chlamydia, gonorrhea, and trich
"""

import numpy as np
import starsim as ss
import sciris as sc
import stisim.utils as ut
import pandas as pd

ss_int_ = ss.dtypes.int
ss_float_ = ss.dtypes.float

__all__ = ['BaseSTI', 'SEIS', 'BaseSTIPars', 'STIPars']


class BaseSTIPars(ss.Pars):
    def __init__(self, **kwargs):
        super().__init__()

        # Settings
        self.include_care = True  # Determines whether testing results are included

        # Time
        self.dt = 'month'

        # Transmission
        self.beta = 0  # Placeholder: no transmission. This will be set in validate_beta
        self.beta_m2f = None
        self.rel_beta_f2m = 0.5
        self.beta_m2c = None
        self.beta_m2m = None
        self.eff_condom = 1
        self.rel_init_prev = 1
        self.update(kwargs)
        return


class STIPars(BaseSTIPars):
    def __init__(self, **kwargs):
        super().__init__()

        # Natural history
        self.dur_exp = ss.constant(ss.weeks(1))  # How long after exposure before you can infect others

        # Symptoms and symptomatic testing
        self.p_symp_dist = ss.bernoulli(p=0.5)  # Distribution of symptomatic vs asymptomatic
        self.p_symp = [0.375, 0.375]
        self.dur_presymp_dist = ss.lognorm_ex()
        self.dur_presymp = [  # For those who develop symptoms, how long before symptoms appear
            [ss.weeks(1), ss.weeks(12)],    # Women
            [ss.weeks(0.25), ss.weeks(3)],  # Men
        ]
        self.p_symp_clear_dist = ss.bernoulli(p=0)
        self.p_symp_clear = [0.0, 0.0]
        self.dur_symp = [
            ss.lognorm_ex(ss.weeks(1), ss.weeks(26)),  # Duration of symptoms
            ss.lognorm_ex(ss.weeks(1), ss.weeks(26)),  # Duration of symptoms
        ]
        self.p_symp_care_dist=ss.bernoulli(p=0)
        self.p_symp_care=[0.3, 0.2]
        self.dur_symp2care_dist=ss.lognorm_ex()
        self.dur_symp2care = [  # For those who test, how long before they seek care - reset for each individual STI
            [ss.weeks(4), ss.weeks(4)],  # Women
            [ss.weeks(6), ss.weeks(4)],  # Men
        ]

        # PID and PID care-seeking
        self.p_pid = ss.bernoulli(p=0.2)
        self.dur_prepid = ss.lognorm_ex(ss.weeks(6), ss.weeks(4))
        self.p_pid_care = ss.bernoulli(p=0.1)  # Women
        self.dur_pid2care = ss.lognorm_ex(ss.weeks(2), ss.weeks(4))  # Women

        # Clearance
        self.dur_asymp2clear_dist = ss.lognorm_ex()
        self.dur_asymp2clear = [  # Duration of untreated asymptomatic infection (excl initial latent)
            [ss.weeks(52), ss.weeks(5)],  # Women
            [ss.weeks(52), ss.weeks(5)],  # Men
        ]
        self.dur_symp2clear_dist = ss.lognorm_ex()
        self.dur_symp2clear = [  # Duration of untreated symptomatic infection (excl initial latent)
            [ss.weeks(52), ss.weeks(5)],  # Women
            [ss.weeks(52), ss.weeks(5)],  # Men
        ]
        self.dur_pid2clear=ss.lognorm_ex(ss.weeks(52), ss.weeks(5))

        # Initial conditions
        self.init_prev=ss.bernoulli(p=0.01)
        self.update(kwargs)
        return


# Main class
class BaseSTI(ss.Infection):
    """
    Base class for sexually transmitted infections.
    Modifies make_new_cases to account for barrier protection.
    """
    def __init__(self, name=None, pars=None, init_prev_data=None, **kwargs):
        super().__init__(name=name)

        # Handle parameters
        default_pars = BaseSTIPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        self.define_states(
            ss.FloatArr('ti_transmitted_sex'),
            ss.FloatArr('ti_transmitted_mtc'),
            ss.FloatArr('ti_transmitted'),
            ss.FloatArr('new_transmissions'),
            ss.FloatArr('new_transmissions_sex'),
        )

        # Set initial prevalence
        self.init_prev_data = init_prev_data

        # Results
        self.age_range = [15, 50]  # Age range for main results e.g. prevalence
        self.age_bins = np.array([0, 15, 20, 25, 30, 35, 50, 65, 100])  # Age bins for results
        self.sex_keys = {'': 'alive', 'f': 'female', 'm': 'male'}

        return

    def make_init_prev(self, uids=None, active=True):
        """ Initialize prevalence by sex and risk group """
        sim = self.sim
        data = self.init_prev_data
        if uids is None: uids = sim.people.auids  # Everyone

        if sc.isnumber(data):
            init_prev = data

        elif isinstance(data, pd.DataFrame):

            init_prev = np.zeros(len(uids), dtype=ss_float_)
            df = data

            nw = sim.networks.structuredsexual
            n_risk_groups = nw.pars.n_risk_groups
            for rg in range(n_risk_groups):
                for sex in ['female', 'male']:
                    for sw in [0, 1]:
                        thisdf = df.loc[(df.risk_group==rg) & (df.sex==sex) & (df.sw==sw)]
                        conditions = sim.people[sex] & (nw.risk_group==rg)
                        if active:
                            conditions = conditions & nw.active(sim.people)
                        if sw:
                            if sex == 'female': conditions = conditions & sim.networks.structuredsexual.fsw
                            if sex == 'male':   conditions = conditions & sim.networks.structuredsexual.client
                        init_prev[conditions[uids]] = thisdf.init_prev.values[0]

        else:
            errormsg = 'Format of init_prev_data must be float or dataframe.'
            raise ValueError(errormsg)

        # Scale and validate
        init_prev = init_prev * self.pars.rel_init_prev
        init_prev = np.clip(init_prev, a_min=0, a_max=1)

        return init_prev

    def validate_beta(self):
        betamap = super().validate_beta()
        if self.pars.beta_m2f is not None and betamap and 'structuredsexual' in betamap.keys():
            betamap['structuredsexual'][0] = self.pars.beta_m2f
            betamap['structuredsexual'][1] = self.pars.beta_m2f * self.pars.rel_beta_f2m
        if self.pars.beta_m2c is not None and betamap and 'maternal' in betamap.keys():
            betamap['maternal'][0] = ss.permonth(self.pars.beta_m2c)
        if self.pars.beta_m2m is not None and betamap and 'msm' in betamap.keys():
            betamap['msm'][0] = self.pars.beta_m2m
            betamap['msm'][1] = self.pars.beta_m2m
        return betamap

    def agehist(self, a):
        """ Return an age histogram """
        aa = self.sim.people.age[a]
        return np.histogram(aa, bins=self.age_bins)[0]

    @property
    def treatable(self):
        """ Assume infected people are treatable, can be overwritten in subclasses """
        return self.infected

    def init_results(self):
        """ Initialize results """
        super().init_results()
        results = sc.autolist()

        # Most results are stored by age and sex
        for sk in self.sex_keys.keys():
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' ({sk.upper()})'
            results += [
                ss.Result(f'incidence{skk}', scale=False, label=f"Incidence{skl}"),
            ]
            if skk != '':
                results += [
                    ss.Result(f'new_infections{skk}', dtype=int, label=f"New infections{skl}"),
                    ss.Result(f'prevalence{skk}', scale=False, label=f"Prevalence{skl}"),
                    ss.Result(f'n_infected{skk}', dtype=int, label=f"Number infected{skl}"),
                ]

            for ab1,ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                ask = f'{skk}_{ab1}_{ab2}'
                asl = f' ({skl}, {ab2}-{ab2})'
                results += [
                    ss.Result(f'new_infections{ask}', dtype=int, label=f"New infections{asl}"),
                    ss.Result(f'n_infected{ask}', dtype=int, label=f"Number infected{asl}"),
                    ss.Result(f'incidence{ask}', scale=False, label=f"Incidence{asl}"),
                    ss.Result(f'prevalence{ask}', scale=False, label=f"Prevalence{asl}"),
                ]

        if self.pars.include_care:
            for sk in self.sex_keys.keys():
                skk = '' if sk == '' else f'_{sk}'
                skl = '' if sk == '' else f' ({sk.upper()})'
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

    def infect(self):
        """ Determine who gets infected on this timestep via transmission on the network """

        new_cases = []
        sources = []
        networks = []
        betamap = self.validate_beta()

        rel_trans = self.rel_trans.asnew(self.infectious * self.rel_trans)
        rel_sus   = self.rel_sus.asnew(self.susceptible * self.rel_sus)

        for i, (nkey,net) in enumerate(self.sim.networks.items()):
            nk = ss.standardize_netkey(nkey)
            if len(net): # Skip networks with no edges
                edges = net.edges
                p1p2b0 = [edges.p1, edges.p2, betamap[nk][0]] # Person 1, person 2, beta 0
                p2p1b1 = [edges.p2, edges.p1, betamap[nk][1]] # Person 2, person 1, beta 1
                for src, trg, beta in [p1p2b0, p2p1b1]:
                    if beta: # Skip networks with no transmission
                        disease_beta = beta.to_prob(self.t.dt) if isinstance(beta, ss.Rate) else beta
                        beta_per_dt = net.net_beta(disease_beta=disease_beta, disease=self) # Compute beta for this network and timestep
                        randvals = self.trans_rng.rvs(src, trg) # Generate a new random number based on the two other random numbers
                        args = (src, trg, rel_trans, rel_sus, beta_per_dt, randvals) # Set up the arguments to calculate transmission
                        target_uids, source_uids = self.compute_transmission(*args) # Actually calculate it
                        new_cases.append(target_uids)
                        sources.append(source_uids)
                        networks.append(np.full(len(target_uids), dtype=ss_int_, fill_value=i))

        # Finalize
        if len(new_cases) and len(sources):
            new_cases = ss.uids.cat(new_cases)
            new_cases, inds = new_cases.unique(return_index=True)
            sources = ss.uids.cat(sources)[inds]
            networks = np.concatenate(networks)[inds]
        else:
            new_cases = ss.uids()
            sources = ss.uids()
            networks = np.empty(0, dtype=ss_int_)

        return new_cases, sources, networks

    def set_outcomes(self, uids, sources=None):
        super().set_outcomes(uids, sources=sources)
        self.new_transmissions[:] = 0  # Total
        self.new_transmissions_sex[:] = 0  # Sexual transmissions only

        if sources is not None:
            # Separate sexual and congenital transmission
            congenital = self.sim.people.age[uids] <= 0
            sex_transmission = sources[~congenital]
            mtc_transmission = sources[congenital]

            # Get the unique sources and counts
            # Only needed for sexual transmission, not MTC.
            unique_sources, counts = np.unique(sources, return_counts=True)
            self.ti_transmitted[unique_sources] = self.ti
            unique_sources_sex, counts_sex = np.unique(sex_transmission, return_counts=True)
            self.ti_transmitted_sex[sex_transmission] = self.ti
            self.ti_transmitted_mtc[mtc_transmission] = self.ti
            self.new_transmissions[unique_sources] = counts
            self.new_transmissions_sex[unique_sources_sex] = counts_sex

        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        ppl = self.sim.people

        adults = (self.sim.people.age >= self.age_range[0]) & (self.sim.people.age <= self.age_range[1])

        # Main results, looping over people keys and attributes
        for pkey, pattr in self.sex_keys.items():
            skk = '' if pkey == '' else f'_{pkey}'

            # Collate results
            new_inf = (self.ti_infected == ti) & ppl[pattr]
            n_sus = self.susceptible & ppl[pattr]
            n_inf = self.infected & ppl[pattr]

            # Store main results
            self.results[f'prevalence{skk}'][ti] = ut.cond_prob(self.infected, adults & ppl[pattr])
            self.results[f'new_infections{skk}'][ti] = ut.count(new_inf)
            self.results[f'n_infected{skk}'][ti] = ut.count(n_inf)
            self.results[f'incidence{skk}'][ti] = ut.countdiv(new_inf, n_sus)

            # Compute age results
            age_results = dict(
                new_infections  = self.agehist(new_inf),
                n_infected      = self.agehist(n_inf),
                incidence       = ut.div(self.agehist(new_inf), self.agehist(n_sus)),
                prevalence      = ut.div(self.agehist(n_inf), self.agehist(ppl[pattr])),
            )

            # Store age results
            for akey, ares in age_results.items():
                ai = 0
                for ab1, ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                    ask = f'{skk}_{ab1}_{ab2}'
                    self.results[f'{akey}{ask}'][ti] = ares[ai]
                    ai += 1
        return


class SEIS(BaseSTI):

    def __init__(self, pars=None, name=None, init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)

        # Handle parameters
        default_pars = STIPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        self.define_states(
            # Natural history
            ss.BoolState('exposed'),
            # ss.BoolState('infected'), # Inherited from ss.Infection
            ss.BoolState('asymptomatic'),
            ss.BoolState('symptomatic'),
            ss.BoolState('pid'),
            ss.BoolState('seeking_care'),
            ss.FloatArr('dur_inf'),
            ss.FloatArr('ti_exposed'),
            ss.FloatArr('ti_symptomatic'),
            ss.FloatArr('ti_symp_clear'),
            ss.FloatArr('ti_seeks_care'),
            ss.FloatArr('ti_pid'),
            ss.FloatArr('ti_clearance'),
        )

        self.age_range = [15, 65]  # Age range for main results e.g. prevalence
        self.age_bins = np.array([0, 15, 20, 25, 30, 35, 50, 65, 100])  # Age bins for results

        return

    @property
    def treatable(self):
        """ Active bacterial presence -- includes exposed and infected, and responds to treatment """
        return self.exposed | self.infected

    def init_results(self):
        """ Initialize results """
        super().init_results()
        results = sc.autolist()

        # Most results are stored by age and sex
        for sk in self.sex_keys.keys():
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' ({sk.upper()})'
            results += [
                ss.Result(f'new_symptomatic{skk}', dtype=int, label=f"New symptomatic{skl}"),
                ss.Result(f'symp_prevalence{skk}', scale=False, label=f"Symptomatic prevalence{skl}"),
            ]
            if skk != '':
                results += [
                    ss.Result(f'n_symptomatic{skk}', dtype=int, label=f"Number symptomatic{skl}"),
                ]

            for ab1,ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                ask = f'{skk}_{ab1}_{ab2}'
                asl = f' ({skl}, {ab2}-{ab2})'
                results += [
                    ss.Result(f'new_symptomatic{ask}', dtype=int, label=f"New symptomatic{asl}"),
                    ss.Result(f'symp_prevalence{ask}', scale=False, label=f"Symptomatic prevalence{asl}"),
                ]

        self.define_results(*results)

        return

    def clear_infection(self, uids):
        self.exposed[uids] = False
        self.infected[uids] = False
        self.symptomatic[uids] = False
        self.asymptomatic[uids] = False
        self.pid[uids] = False
        self.seeking_care[uids] = False
        self.susceptible[uids] = True
        past_care_seekers = uids[(self.ti_seeks_care[uids] < self.ti).nonzero()[-1]]
        self.ti_seeks_care[past_care_seekers] = np.nan
        self.ti_clearance[uids] = self.ti
        self.dur_inf[uids] = self.ti - self.ti_infected[uids]

    def step_state(self):
        """ Updates for this timestep """
        ti = self.ti

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Exposed -> infected
        new_infected = (self.exposed & (self.ti_infected <= ti)).uids
        if len(new_infected):
            self.exposed[new_infected] = False
            self.infected[new_infected] = True
            self.ti_infected[new_infected] = ti

        # Presymptomatic -> symptomatic
        new_symptomatic = (self.asymptomatic & (self.ti_symptomatic <= ti)).uids
        if len(new_symptomatic):
            self.asymptomatic[new_symptomatic] = False
            self.symptomatic[new_symptomatic] = True
            self.ti_symptomatic[new_symptomatic] = ti

        # Clear infections
        new_cleared = (self.infected & (self.ti_clearance <= ti)).uids
        self.clear_infection(new_cleared)

        # Progress PID
        new_pid = (~self.pid & (self.ti_pid <= ti)).uids
        self.pid[new_pid] = True
        self.ti_pid[new_pid] = ti

        # Symptomatic/PID care seeking
        old_seekers = (self.seeking_care).uids
        self.seeking_care[old_seekers] = False
        self.ti_seeks_care[old_seekers] = np.nan  # Remove the old
        new_seekers = (self.infected & (self.ti_seeks_care <= ti)).uids
        self.seeking_care[new_seekers] = True
        self.ti_seeks_care[new_seekers] = ti
        # If we don't update this results here, it's possible that someone could seek care,
        # then get reinfected, and the reinfection would wipe all their dates so we'd miss
        # counting them in the results.
        if self.pars.include_care:
            for sk, sl in self.sex_keys.items():
                skk = '' if sk == '' else f'_{sk}'
                self.results['new_care_seekers'+skk][ti] = ut.count((self.ti_seeks_care == ti) & self.sim.people[sl])

        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        ppl = self.sim.people

        adults = (self.sim.people.age >= self.age_range[0]) & (self.sim.people.age <= self.age_range[1])

        # Main results, looping over people keys and attributes
        for pkey, pattr in self.sex_keys.items():
            skk = '' if pkey == '' else f'_{pkey}'

            # Collate results
            new_sym = (self.ti_symptomatic == ti) & ppl[pattr]
            n_sym = self.symptomatic & ppl[pattr]

            # Store main results
            self.results[f'symp_prevalence{skk}'][ti] = ut.count(self.symptomatic & adults & ppl[pattr]) / ut.count(adults & ppl[pattr])
            self.results[f'new_symptomatic{skk}'][ti] = ut.count(new_sym)
            self.results[f'n_symptomatic{skk}'][ti] = ut.count(n_sym)

            # Compute age results
            age_results = dict(
                new_symptomatic=self.agehist(new_sym),
                symp_prevalence=ut.div(self.agehist(n_sym), self.agehist(ppl[pattr]))
            )

            # Store age results
            for akey, ares in age_results.items():
                ai = 0
                for ab1, ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                    ask = f'{skk}_{ab1}_{ab2}'
                    self.results[f'{akey}{ask}'][ti] = ares[ai]
                    ai += 1
        return

    def set_exposure(self, uids):
        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.asymptomatic[uids] = True
        self.ti_exposed[uids] = self.ti
        dur_exp = self.pars.dur_exp.rvs(uids)
        self.ti_infected[uids] = self.ti + dur_exp
        return

    def set_pars(self, par, uids):
        if sc.isnumber(par[0]):  # Single par, e.g. for bernoulli
            arr = np.full(len(uids), np.nan)
            arr[self.sim.people.female[uids]] = par[0]
            arr[self.sim.people.male[uids]] = par[1]
        elif len(par[0]) > 1:
            arr = np.full((len(uids), 2), np.nan)
            for idx in range(2):
                arr[self.sim.people.female[uids], idx] = par[0][idx]
                arr[self.sim.people.male[uids], idx] = par[1][idx]
        return arr

    def set_symptoms(self, p, uids):
        # Probability of symptoms
        p_symp_vals = self.set_pars(p.p_symp, uids)
        p.p_symp_dist.set(p_symp_vals)
        symp, asymp = p.p_symp_dist.split(uids)

        # Duration of presymptomatic period
        dur_presymp_vals = self.set_pars(p.dur_presymp, symp)
        p.dur_presymp_dist.set(mean=dur_presymp_vals[:,0], std=dur_presymp_vals[:,1])
        dur_presymp = self.pars.dur_presymp_dist.rvs(symp)
        self.ti_symptomatic[symp] = self.ti_infected[symp] + dur_presymp
        return symp, asymp

    def set_care_seeking(self, p, symp):
        # Probability of symptoms
        p_care_vals = self.set_pars(p.p_symp_care, symp)
        p.p_symp_care_dist.set(p_care_vals)
        care, no_care = p.p_symp_care_dist.split(symp)

        # Duration before seeking care
        dur_precare_vals = self.set_pars(p.dur_symp2care, care)
        p.dur_symp2care_dist.set(mean=dur_precare_vals[:, 0], std=dur_precare_vals[:, 1])
        dur_symp2care = p.dur_symp2care_dist.rvs(care)
        self.ti_seeks_care[care] = self.ti_symptomatic[care] + dur_symp2care
        return

    def set_pid(self, p, f_uids):
        pid = p.p_pid.filter(f_uids)
        dur_prepid = p.dur_prepid.rvs(pid)
        self.ti_pid[pid] = self.ti_infected[pid] + dur_prepid
        return pid

    def set_pid_care_seeking(self, p, pid):
        pid_care = p.p_pid_care.filter(pid)
        dur_pid2care = p.dur_pid2care.rvs(pid_care)
        self.ti_seeks_care[pid_care] = np.minimum(self.ti_seeks_care[pid_care], self.ti_infected[pid_care] + dur_pid2care)
        return

    def set_duration(self, p, symp, asymp, pid):
        dur_symp2clear_vals = self.set_pars(p.dur_symp2clear, symp)
        p.dur_symp2clear_dist.set(mean=dur_symp2clear_vals[:, 0], std=dur_symp2clear_vals[:, 1])
        dur_symp2clear = p.dur_symp2clear_dist.rvs(symp)
        self.ti_clearance[symp] = self.ti_symptomatic[symp] + dur_symp2clear

        dur_asymp2clear_vals = self.set_pars(p.dur_asymp2clear, asymp)
        p.dur_asymp2clear_dist.set(mean=dur_asymp2clear_vals[:, 0], std=dur_asymp2clear_vals[:, 1])
        dur_asymp2clear = p.dur_asymp2clear_dist.rvs(asymp)
        self.ti_clearance[asymp] = self.ti_infected[asymp] + dur_asymp2clear

        dur_inf_pid = p.dur_pid2clear.rvs(pid)
        self.ti_clearance[pid] = np.maximum(self.ti_clearance[pid], dur_inf_pid + self.ti_pid[pid])
        return

    def wipe_dates(self, uids):
        """ Clear all previous dates """
        self.ti_exposed[uids] = np.nan
        self.ti_infected[uids] = np.nan
        self.ti_symptomatic[uids] = np.nan
        self.ti_symp_clear[uids] = np.nan
        self.ti_pid[uids] = np.nan
        # self.ti_seeks_care[uids] = np.nan
        self.ti_clearance[uids] = np.nan
        # self.dur_inf[uids] = np.nan
        return

    def set_prognoses(self, uids, sources=None):
        """
        Set initial prognoses for adults newly infected
        """
        super().set_prognoses(uids, sources)
        self.wipe_dates(uids)
        self.dur_inf[uids] = np.nan  # Overwrite. Not done in wipe_dates because that's also called during treatment

        ppl = self.sim.people
        p = self.pars
        f_uids = ppl.female.uids.intersect(uids)  # For PID

        self.set_exposure(uids)  # Set exposure
        symp, asymp = self.set_symptoms(p, uids)  # Set symptoms & presymptomatic duration
        self.set_care_seeking(p, symp)  # Determine who seeks care and when
        pid = self.set_pid(p, f_uids)  # Determine who developes PID and when
        self.set_pid_care_seeking(p, pid)
        self.set_duration(p, symp, asymp, pid)

        # Determine overall duration of infection, but don't set until clearance
        dur_inf = self.ti_clearance[uids] - self.ti_infected[uids]

        if (dur_inf < 0).any():
            errormsg = 'Invalid durations of infection'
            raise ValueError(errormsg)

        return

