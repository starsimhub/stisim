"""
What STIsim's diseases share: results by sex and age (`STIResults`), initial prevalence by group, and `SEIS`,
the declared natural history of the curable STIs (chlamydia, gonorrhea, and trichomoniasis)
"""

import numpy as np
import starsim as ss
import sciris as sc
import stisim.utils as ut
import pandas as pd

ss_float_ = ss.dtypes.float

__all__ = ['STIResults', 'SEIS', 'init_prev_by_group', 'sexual_beta', 'hetero_networks', 'msm_networks', 'by_sex', 'lognorm_by_sex', 'default_age_bins', 'default_sex_keys']

# Default stratifications for disease results. Pass age_bins=None or sex_keys=None to a
# disease module to skip creating the corresponding stratified results.
default_age_bins = [0, 15, 20, 25, 30, 35, 50, 65, 100]
default_sex_keys = {'': 'alive', 'f': 'female', 'm': 'male'}


hetero_networks = ['structuredsexual', 'mfnetwork', 'swnetwork', 'sexualnet', 'sexualnetwork'] # The sexual networks between men and women, by name: STIsim's, Starsim's (ss.SexualNet), and others named like it (e.g. HPVsim's)
msm_networks = ['msm', 'msmscalefreenetwork'] # ... and between men


def sexual_beta(m2f, f2m, m2m=0, hetero=None, msm=None):
    """
    A beta for each of STIsim's sexual networks, by name: a probability per act from a man to a woman and the reverse
    on the networks between men and women, and one between men on the others; 0 on any other network (e.g. maternal)

    Args:
        m2f (float/par): from a man to a woman
        f2m (float/par): from a woman to a man
        m2m (float/par): between men
        hetero (list): the networks between men and women (default `sti.hetero_networks`)
        msm (list): the networks between men (default `sti.msm_networks`)
    """
    beta = {name:ss.peract(m2f=m2f, f2m=f2m) for name in sc.ifelse(hetero, hetero_networks)}
    beta |= {name:ss.peract(m2m) for name in sc.ifelse(msm, msm_networks)}
    return beta | dict(default=0)


def by_sex(female, male):
    """
    A value for women and one for men, read from each person's sex

    Shorthand for `ss.by(sex=dict(female=female, male=male))`.

    Examples:
        ```python
        ng = sti.Gonorrhea(p_symp=sti.by_sex(0.35, 0.65))
        ```
    """
    return ss.by(sex=dict(female=female, male=male))


def lognorm_by_sex(female, male):
    """
    A lognormal duration whose mean and standard deviation differ by sex, each given as `[mean, std]`

    Examples:
        ```python
        dur = sti.lognorm_by_sex(female=[ss.months(8), ss.months(2)], male=[ss.months(6), ss.months(3)])
        ```
    """
    return ss.lognorm_ex(mean=by_sex(female[0], male[0]), std=by_sex(female[1], male[1]))


class STIResults:
    """
    Results by sex and age group, shared by STIsim's diseases

    A disease that uses it calls `set_strata()` when it is made; in a sim, it provides `infected` and
    `susceptible` (boolean arrays) and `new_infected()`, the people infected in the current step.
    """
    def set_strata(self, age_range, age_bins, sex_keys):
        """ The adult ages of the main results (e.g. prevalence), the age groups (or None), and the sex strata (or None, for the whole population only) """
        self.age_range = age_range
        self.age_bins = np.array(age_bins) if age_bins is not None else None
        self.sex_keys = dict(sex_keys) if sex_keys is not None else {'': 'alive'}
        return

    def sex_strata(self):
        """ For each sex stratum: the suffix of its results ('', '_f', or '_m'), its label, and the attribute of people it selects """
        for sk, attr in self.sex_keys.items():
            yield ('' if sk == '' else f'_{sk}'), ('' if sk == '' else f' ({sk.upper()})'), attr

    def age_strata(self):
        """ The age groups, as pairs of bin edges """
        return list(zip(self.age_bins[:-1], self.age_bins[1:])) if self.age_bins is not None else []

    def agehist(self, a):
        """ Return an age histogram """
        aa = self.sim.people.age[a]
        w = ut.weights(a)
        return np.histogram(aa, bins=self.age_bins, weights=None if w is None else w[np.asarray(a.values, dtype=bool)])[0] # Each agent counted by its weight, if they differ

    def adults(self):
        """ Who is in the age range of the main results """
        age = self.sim.people.age
        return (age >= self.age_range[0]) & (age <= self.age_range[1])

    def transmitters(self):
        """ Who infected others sexually in this step, and how many each: by the sexual transmission, `infection` """
        sources, _ = self.engine.sources(self.infection.uids)
        uids, counts = np.unique(sources[~np.isnan(sources)].astype(int), return_counts=True)
        return ss.uids(uids), counts

    def define_sti_results(self, include_care=True):
        """ Define the results of infection by sex and age, of care, and of the routes of transmission (those already defined, e.g. by `ss.Infection`, are kept) """
        results = sc.autolist()
        for skk, skl, attr in self.sex_strata():
            is_subpop = skk != ''  # Sex-stratified results are not auto-plotted
            results += [
                ss.Result(f'incidence{skk}', scale=False, label=f"Incidence{skl}", auto_plot=not is_subpop),
                ss.Result(f'new_infections{skk}', dtype=int, label=f"New infections{skl}", auto_plot=False),
                ss.Result(f'prevalence{skk}', scale=False, label=f"Prevalence{skl}", auto_plot=not is_subpop),
                ss.Result(f'infected{skk}', dtype=int, label=f"Number infected{skl}", auto_plot=False),
            ]
            for ab1, ab2 in self.age_strata():
                ask = f'{skk}_{ab1}_{ab2}'
                asl = f' ({skl}, {ab2}-{ab2})'
                results += [
                    ss.Result(f'new_infections{ask}', dtype=int, label=f"New infections{asl}", auto_plot=False),
                    ss.Result(f'infected{ask}', dtype=int, label=f"Number infected{asl}", auto_plot=False),
                    ss.Result(f'incidence{ask}', scale=False, label=f"Incidence{asl}", auto_plot=False),
                    ss.Result(f'prevalence{ask}', scale=False, label=f"Prevalence{asl}", auto_plot=False),
                ]

        # Results of care, which the testing and treatment interventions write
        if include_care:
            for skk, skl, attr in self.sex_strata():
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

        # Transmission route results
        results += [
            ss.Result('new_infections_sex', dtype=int, label='New infections (sexual)', auto_plot=False),
            ss.Result('new_infections_mtct', dtype=int, label='New infections (MTCT)', auto_plot=False),
            ss.Result('new_infections_prenatal', dtype=int, label='New infections (prenatal MTCT)', auto_plot=False),
            ss.Result('new_infections_postnatal', dtype=int, label='New infections (postnatal MTCT)', auto_plot=False),
        ]

        declared = set(self.graph.exprs) | set(self.graph.outcomes) if self.graph is not None else set() # Those declared as results of the states and flows, made by the engine
        self.define_results(*[res for res in results if res.name not in self.results and res.name not in declared])
        return

    def update_sti_results(self):
        """ Count infections by sex and age; incidence is counted at acquisition, including the infections seeded at the start """
        ti = self.ti
        ppl = self.sim.people
        adults = self.adults()
        new_inf = self.new_infected()

        for skk, skl, attr in self.sex_strata():
            new_inf_s = new_inf & ppl[attr]
            sus_s = self.susceptible & ppl[attr]
            inf_s = self.infected & ppl[attr]

            self.results[f'prevalence{skk}'][ti] = ut.cond_prob(self.infected, adults & ppl[attr])
            self.results[f'new_infections{skk}'][ti] = ut.count(new_inf_s)
            self.results[f'infected{skk}'][ti] = ut.count(inf_s)
            self.results[f'incidence{skk}'][ti] = ut.countdiv(new_inf_s, sus_s)

            if self.age_bins is not None:
                age_results = dict(
                    new_infections = self.agehist(new_inf_s),
                    infected       = self.agehist(inf_s),
                    incidence      = ut.div(self.agehist(new_inf_s), self.agehist(sus_s)),
                    prevalence     = ut.div(self.agehist(inf_s), self.agehist(ppl[attr])),
                )
                for akey, ares in age_results.items():
                    for ai, (ab1, ab2) in enumerate(self.age_strata()):
                        self.results[f'{akey}{skk}_{ab1}_{ab2}'][ti] = ares[ai]
        return


def init_prev_by_group(sim, uids, data, rel_init_prev=1, active=True):
    """
    The initial prevalence of each of some agents, by risk group, sex, and sex work, from data

    Args:
        sim (Sim): the sim
        uids (uids): the agents
        data (float/DataFrame): a prevalence, or a DataFrame with columns risk_group, sex, sw, and init_prev
        rel_init_prev (float): a multiplier on it
        active (bool): only for the sexually active
    """
    if sc.isnumber(data):
        init_prev = np.full(len(uids), data, dtype=ss_float_)
    elif isinstance(data, pd.DataFrame):
        init_prev = np.zeros(len(uids), dtype=ss_float_)
        nw = sim.networks.structuredsexual
        for rg in range(nw.pars.n_risk_groups):
            for sex in ['female', 'male']:
                for sw in [0, 1]:
                    thisdf = data.loc[(data.risk_group==rg) & (data.sex==sex) & (data.sw==sw)]
                    conditions = sim.people[sex] & (nw.risk_group==rg)
                    if active:
                        conditions = conditions & nw.active(sim.people)
                    if sw:
                        conditions = conditions & (nw.fsw if sex == 'female' else nw.client)
                    init_prev[conditions[uids]] = thisdf.init_prev.values[0]
    else:
        errormsg = 'Format of init_prev_data must be float or dataframe.'
        raise ValueError(errormsg)
    return np.clip(init_prev * rel_init_prev, a_min=0, a_max=1)


class SEIS(ss.Disease, STIResults):
    """
    The natural history of a curable STI, declared: susceptible, exposed, infectious, and susceptible again

    People infected go through a latent period (`exposed`), then become infectious, with symptoms
    to come (`presymptomatic`) or with none (`asymptomatic`). Those with symptoms may seek care,
    and women may develop pelvic inflammatory disease (`pid`). Infections clear on their own, and
    everyone who clears is susceptible again. Used for chlamydia, gonorrhea, and trichomoniasis,
    which differ only in their parameters.

    States:

    - `susceptible`
    - `exposed`: infected, not yet infectious
    - `presymptomatic`: infectious, with symptoms to come
    - `asymptomatic`: infectious, with no symptoms to come
    - `symptomatic`: infectious, with symptoms
    - `pid`: pelvic inflammatory disease (women)

    `infectious` and `infected` (`exposed` too) are conditions over these.

    Transitions (each is an event that interventions can be triggered by, e.g. `who=ng.seek_care`):

    - `infection`: sexual transmission, from someone infectious
    - `mtct`: transmission from mother to child in pregnancy (if `beta_m2c` is given)
    - `onset`: the end of the latent period, `dur_exp` after infection, to `presymptomatic` (with probability `p_symp`) or `asymptomatic`
    - `symptoms`: `dur_presymp` after becoming infectious
    - `seek_care`: a share `p_symp_care` of the symptomatic seek care, `dur_symp2care` after their symptoms start, if they still have them; they stay symptomatic
    - `clear_symp`: `dur_symp2clear` after symptoms start, whether or not they sought care
    - `clear_asymp`: `dur_asymp2clear` after becoming infectious
    - `persist`: a share `1 - p_clear` of asymptomatic women do not clear, for `dur_persist` (100 years by default)
    - `pid_onset`: a share `p_pid` of infected women develop PID, `dur_prepid` after becoming infectious
    - `seek_care_pid`: a share `p_pid_care` of women with PID seek care, `dur_pid2care` after PID starts
    - `clear_pid`: `dur_pid2clear` after PID starts

    Values that differ by sex are given with `sti.by_sex()` (e.g. `p_symp=sti.by_sex(0.35, 0.65)`), and
    lognormal durations that do with `sti.lognorm_by_sex()`.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module
        age_range (list): age range of the main (adult) results, e.g. prevalence (default [15, 65])
        age_bins (list): bin edges for age-stratified results; ``None`` skips them
        sex_keys (dict): sex strata for results; ``None`` skips sex-stratified results
        include_care (bool): whether to define the results of care (care seekers, test outcomes, treatments)
        kwargs (dict): parameter values (e.g. `beta_m2f=0.06`), `init_prev`, or `label`

    Examples:
        ```python
        ng = sti.Gonorrhea(beta_m2f=0.08, init_prev=0.02)
        sim = sti.Sim(diseases=ng, n_agents=5000)
        ```
    """
    def __init__(self, pars=None, name=None, age_range=None, age_bins=default_age_bins,
                 sex_keys=default_sex_keys, include_care=True, **kwargs):

        # Transmission: a probability per act on the sexual networks, and a rate per month from mother to child
        beta_m2f = ss.par(0.0, doc='Probability of transmission per act from a man to a woman')
        rel_beta_f2m = ss.par(0.5, doc='Transmission from a woman to a man, relative to from a man to a woman')
        beta_m2m = ss.par(0.0, doc='Probability of transmission per act between men')
        beta_m2c = ss.par(0.0, doc='Rate of transmission from mother to child in pregnancy, per month')

        # Natural history
        dur_exp = ss.par(ss.constant(ss.weeks(1)), doc='Latent period: time from infection to being infectious')
        p_symp = ss.par(by_sex(0.375, 0.375), doc='Probability that an infection becomes symptomatic')
        dur_presymp = ss.par(lognorm_by_sex([ss.weeks(1), ss.weeks(12)], [ss.weeks(0.25), ss.weeks(3)]), doc='Time from being infectious to symptoms')
        p_symp_care = ss.par(by_sex(0.3, 0.2), doc='Probability that someone with symptoms seeks care')
        dur_symp2care = ss.par(lognorm_by_sex([ss.weeks(4), ss.weeks(4)], [ss.weeks(6), ss.weeks(4)]), name='dur_symp2care', doc='Time from symptoms to seeking care')
        dur_symp2clear = ss.par(lognorm_by_sex([ss.weeks(52), ss.weeks(5)], [ss.weeks(52), ss.weeks(5)]), name='dur_symp2clear', doc='Time from symptoms to clearance')
        dur_asymp2clear = ss.par(lognorm_by_sex([ss.weeks(52), ss.weeks(5)], [ss.weeks(52), ss.weeks(5)]), doc='Time from being infectious to clearance, without symptoms')
        p_clear = ss.par(1.0, doc='Probability that an asymptomatic infection in a woman clears (the rest persist)')
        dur_persist = ss.par(ss.constant(ss.years(100)), doc='Duration of a persistent infection')
        p_pid = ss.par(0.2, doc='Probability that an infected woman develops pelvic inflammatory disease (PID)')
        dur_prepid = ss.par(ss.lognorm_ex(ss.weeks(6), ss.weeks(4)), name='dur_prepid', doc='Time from being infectious to PID')
        p_pid_care = ss.par(0.1, doc='Probability that a woman with PID seeks care')
        dur_pid2care = ss.par(ss.lognorm_ex(ss.weeks(2), ss.weeks(4)), name='dur_pid2care', doc='Time from PID to seeking care')
        dur_pid2clear = ss.par(ss.lognorm_ex(ss.weeks(52), ss.weeks(5)), doc='Time from PID to clearance')

        # States
        susceptible, exposed, presymptomatic, asymptomatic, symptomatic, pid = ss.states('susceptible exposed presymptomatic asymptomatic symptomatic pid')
        infectious = presymptomatic + asymptomatic + symptomatic + pid
        female = ss.ppl.female

        # Transmission
        sexual = sexual_beta(beta_m2f, beta_m2f*rel_beta_f2m, beta_m2m)
        infection = susceptible >> exposed @ sexual*infectious
        mtct = susceptible >> exposed @ beta_m2c*dict(maternalnet=ss.permonth(1), default=0)*infectious

        # Natural history
        onset = ss.Transition(exposed, [presymptomatic, asymptomatic], p=[p_symp, None], dur=dur_exp)
        symptoms = presymptomatic >> symptomatic @ dur_presymp
        seek_care = ss.Transition(symptomatic, [symptomatic, ss.stay], p=[p_symp_care, None], on=symptoms, dur=dur_symp2care) # A share seek care, and stay symptomatic
        clear_symp = ss.Transition(symptomatic, susceptible, on=symptoms, dur=dur_symp2clear) # Timed from the symptoms, so seeking care does not restart it
        clear_asymp = asymptomatic >> susceptible @ dur_asymp2clear
        persist = ss.Transition(asymptomatic, susceptible, p=1 - p_clear, dur=dur_persist, where=female) # The rest of the asymptomatic clear by clear_asymp
        pid_onset = ss.Transition(presymptomatic + asymptomatic + symptomatic, [pid, ss.stay], p=[p_pid, None], on=onset, dur=dur_prepid, where=female)
        seek_care_pid = ss.Transition(pid, [pid, ss.stay], p=[p_pid_care, None], on=pid_onset, dur=dur_pid2care)
        clear_pid = pid >> susceptible @ dur_pid2clear

        # The parameters that the transitions do not read: the condoms' efficacy, which the sexual networks read (see sti.BaseNetwork.net_beta());
        # and, as the graph does not yet register them [PROPOSAL ST-D], the times after an event (on=), named above, whose values are set once it is built
        values = sc.mergedicts(dict(init_prev=0.01), pars, kwargs)
        delays = dict(dur_symp2care=dur_symp2care, dur_symp2clear=dur_symp2clear, dur_prepid=dur_prepid, dur_pid2care=dur_pid2care)
        given = {key:values.pop(key) for key in list(delays) + ['eff_condom'] if key in values}

        # [PROPOSAL ST-A] Declared here rather than in the class body, since on= a sibling transition does not yet work there
        super().__init__(
            infection=infection, mtct=mtct, onset=onset, symptoms=symptoms, seek_care=seek_care,
            clear_symp=clear_symp, clear_asymp=clear_asymp, persist=persist,
            pid_onset=pid_onset, seek_care_pid=seek_care_pid, clear_pid=clear_pid,
            new_infections_sex=1*infection, new_infections_prenatal=1*mtct, new_infections_mtct=1*mtct, # Results of the flows, by route
            name=name, **sc.mergedicts(dict(dt='month'), values), # Monthly unless dt is given
        )
        self.define_pars(eff_condom=1.0, **delays)
        self.pars.update(given)
        self.set_strata(age_range if age_range is not None else [15, 65], age_bins, sex_keys)
        self.include_care = include_care
        return

    @property
    def infectious(self):
        """ Who is infectious: all but the latent period """
        return self.presymptomatic | self.asymptomatic | self.symptomatic | self.pid

    @property
    def infected(self):
        """ Who is infected, including the latent period """
        return self.exposed | self.infectious

    @property
    def treatable(self):
        """ Who a treatment can cure: everyone infected """
        return self.infected

    def new_infected(self):
        """ Who was infected in this step: by transmission, or seeded at the start """
        new = self.infection.uids | self.mtct.uids
        if self.ti == 0:
            new = new | self.infected.uids
        return self.sim.people.uid.isin(new) # As a boolean array

    def clear_infection(self, uids):
        """ Cure people (e.g. treated): they are susceptible again """
        self.move(uids, self.susceptible)
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        self.define_results(ss.Result('infectious', dtype=int, label='Number infectious', auto_plot=False))
        self.define_sti_results(include_care=self.include_care)
        self.define_results(ss.Result('cum_infections', dtype=int, label='Cumulative infections'))

        # Symptoms, by sex and age
        results = sc.autolist()
        for skk, skl, attr in self.sex_strata():
            is_subpop = skk != ''
            results += [
                ss.Result(f'new_symptomatic{skk}', dtype=int, label=f"New symptomatic{skl}", auto_plot=not is_subpop),
                ss.Result(f'symp_prevalence{skk}', scale=False, label=f"Symptomatic prevalence{skl}", auto_plot=not is_subpop),
            ]
            if is_subpop:
                results += ss.Result(f'symptomatic{skk}', dtype=int, label=f"Number symptomatic{skl}", auto_plot=False)
            for ab1, ab2 in self.age_strata():
                ask = f'{skk}_{ab1}_{ab2}'
                asl = f' ({skl}, {ab2}-{ab2})'
                results += [
                    ss.Result(f'new_symptomatic{ask}', dtype=int, label=f"New symptomatic{asl}", auto_plot=False),
                    ss.Result(f'symp_prevalence{ask}', scale=False, label=f"Symptomatic prevalence{asl}", auto_plot=False),
                ]
        self.define_results(*results)
        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        ppl = self.sim.people
        res = self.results
        self.update_sti_results()
        res.infectious[ti] = ut.count(self.infectious)

        # Symptoms and care seeking, by sex and age
        adults = self.adults()
        new_sym = ppl.uid.isin(self.symptoms.uids)
        seekers = ppl.uid.isin(self.seek_care.uids | self.seek_care_pid.uids)
        for skk, skl, attr in self.sex_strata():
            new_sym_s = new_sym & ppl[attr]
            sym_s = self.symptomatic & ppl[attr]
            res[f'symp_prevalence{skk}'][ti] = ut.cond_prob(self.symptomatic, adults & ppl[attr])
            res[f'new_symptomatic{skk}'][ti] = ut.count(new_sym_s)
            if skk:
                res[f'symptomatic{skk}'][ti] = ut.count(sym_s)
            if self.include_care:
                res[f'new_care_seekers{skk}'][ti] = ut.count(seekers & ppl[attr])
            if self.age_bins is not None:
                age_results = dict(
                    new_symptomatic = self.agehist(new_sym_s),
                    symp_prevalence = ut.div(self.agehist(sym_s), self.agehist(ppl[attr])),
                )
                for akey, ares in age_results.items():
                    for ai, (ab1, ab2) in enumerate(self.age_strata()):
                        res[f'{akey}{skk}_{ab1}_{ab2}'][ti] = ares[ai]
        return

    def finalize_results(self):
        super().finalize_results()
        self.results.cum_infections[:] = np.cumsum(self.results.new_infections[:])
        return
