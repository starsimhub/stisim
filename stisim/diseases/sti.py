"""
What STIsim's diseases share: results by sex and age (`STIResults`), initial prevalence by group, and `SEIS`,
the natural history of the curable STIs (chlamydia, gonorrhea, and trichomoniasis)
"""

import numpy as np
import starsim as ss
import sciris as sc
import stisim.utils as ut
import pandas as pd
import scipy.special as sps

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
    sexual_sources = ss.uids() # Who infected someone sexually in this step, as the disease's step() records them: see transmitters()

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
        """ Who infected others sexually in this step, and how many each (from `sexual_sources`, which the disease's `step()` records) """
        uids, counts = np.unique(self.sexual_sources, return_counts=True)
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

        self.define_results(*[res for res in results if res.name not in self.results])
        return

    def tally(self):
        """
        How this step's results count the living agents, each count in one pass: `vals(arr)`, their values as booleans; `count(mask)`, the
        number for whom it holds (by weight, if agents' weights differ); `hist(mask)`, the same by age group (as `np.histogram()` gives it,
        the last group including its upper edge); `div(a, b)` and `hdiv(a, b)`, a ratio of counts or of histograms (0 where b is); and
        `adults`, who is in the age range of the main results
        """
        ppl = self.sim.people
        t = sc.objdict()
        t.vals = lambda arr: np.asarray(arr.values, dtype=bool)
        age = np.asarray(ppl.age.values)
        w = ut.weights(ppl.alive)
        t.adults = (age >= self.age_range[0]) & (age <= self.age_range[1])
        t.count = lambda mask: np.count_nonzero(mask) if w is None else float(w[mask].sum())
        t.div = lambda a, b: a/b if b else 0.0
        bins = self.age_bins
        if bins is not None:
            n_bins = len(bins) - 1
            band = np.searchsorted(bins, age, side='right') - 1
            band[age == bins[-1]] = n_bins - 1
            band[(band < 0) | (band >= n_bins) | np.isnan(age)] = n_bins # Outside the groups
            t.hist = lambda mask: np.bincount(band[mask], weights=None if w is None else w[mask], minlength=n_bins + 1)[:n_bins]
            t.hdiv = lambda a, b: np.divide(a, b, out=np.zeros(n_bins), where=b != 0)
        return t

    def set_age_results(self, skk, **values):
        """ Store results by age group (one value per group) for a sex stratum """
        for akey, ares in values.items():
            for ai, (ab1, ab2) in enumerate(self.age_strata()):
                self.results[f'{akey}{skk}_{ab1}_{ab2}'][self.ti] = ares[ai]
        return

    def update_sti_results(self, t=None):
        """ Count infections by sex and age (with `t` from `tally()`); incidence is counted at acquisition, including the infections seeded at the start """
        ti = self.ti
        ppl = self.sim.people
        t = t or self.tally()
        infected, new, susceptible = t.vals(self.infected), t.vals(self.new_infected()), t.vals(self.susceptible)
        for skk, skl, attr in self.sex_strata():
            sex = t.vals(ppl[attr])
            new_s = new & sex
            sus_s = susceptible & sex
            inf_s = infected & sex
            n_new = t.count(new_s)
            self.results[f'prevalence{skk}'][ti] = t.div(t.count(inf_s & t.adults), t.count(t.adults & sex))
            if skk or ti == 0: # For everyone, ss.Infection counts the infections as they happen (including those of people who die in the step), but not those seeded at the start
                self.results[f'new_infections{skk}'][ti] = n_new
            self.results[f'infected{skk}'][ti] = t.count(inf_s)
            self.results[f'incidence{skk}'][ti] = t.div(n_new, t.count(sus_s))
            if self.age_bins is not None:
                h_new, h_inf = t.hist(new_s), t.hist(inf_s)
                self.set_age_results(skk, new_infections=h_new, infected=h_inf, incidence=t.hdiv(h_new, t.hist(sus_s)), prevalence=t.hdiv(h_inf, t.hist(sex)))
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


class SEIS(ss.Infection, STIResults):
    """
    The natural history of a curable STI: susceptible, exposed, infectious, and susceptible again

    People infected go through a latent period (`exposed`), then become infectious, with symptoms
    to come (`presymptomatic`) or with none (`asymptomatic`). Those with symptoms may seek care,
    and women may develop pelvic inflammatory disease (`pid`). Infections clear on their own, and
    everyone who clears is susceptible again. Used for chlamydia, gonorrhea, and trichomoniasis,
    which differ only in their parameters.

    Each infection's course is drawn when it starts (`set_prognoses()`): when it becomes infectious,
    whether and when symptoms start, care is sought, and PID starts, and when it clears. Each step
    then takes the people whose time has come (`step_state()`), in the order of the course, so
    someone can pass through several stages in one step.

    States: `susceptible`, `exposed`, `presymptomatic`, `asymptomatic`, `symptomatic`, `pid`, and
    `infected` and `infectious` (`exposed` and the rest, and all but `exposed`). The people who
    sought care in the step are `seek_care` (with symptoms) and `seek_care_pid` (with PID), and
    those whose symptoms started, `symptoms`: e.g. `who=ng.seek_care` for an `ss.treat`, or
    `eligibility=lambda sim: sim.diseases.ng.seek_care.uids` for an STIsim treatment.

    The course, by its parameters:

    - `dur_exp` after infection, infectious, with symptoms to come with probability `p_symp`
    - `dur_presymp` after that, symptoms; a share `p_symp_care` seek care `dur_symp2care` after their symptoms start, if they still have them
    - clearance: `dur_symp2clear` after symptoms start, or `dur_asymp2clear` after becoming infectious without symptoms; a share `1 - p_clear` of asymptomatic women do not clear, for `dur_persist` (100 years by default)
    - PID: a share `p_pid` of infected women, `dur_prepid` after becoming infectious, if still infected; a share `p_pid_care` of them seek care `dur_pid2care` after it starts; it clears `dur_pid2clear` after it starts

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
    durs = ['dur_exp', 'dur_presymp', 'dur_symp2care', 'dur_symp2clear', 'dur_asymp2clear', 'dur_persist', 'dur_prepid', 'dur_pid2care', 'dur_pid2clear']
    chances = ['p_symp', 'p_symp_care', 'p_clear', 'p_pid', 'p_pid_care']

    def __init__(self, pars=None, name=None, age_range=None, age_bins=default_age_bins,
                 sex_keys=default_sex_keys, include_care=True, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.define_pars(
            dt = 'month',
            init_prev = 0.01,

            # Transmission: a probability per act on the sexual networks, and a rate per month from mother to child
            beta = 0.0, # On each network, from the ones below (see validate_beta())
            beta_m2f = 0.0, # From a man to a woman
            rel_beta_f2m = 0.5, # From a woman to a man, relative to from a man to a woman
            beta_m2m = 0.0, # Between men
            beta_m2c = 0.0, # From mother to child in pregnancy, per month
            eff_condom = 1.0, # Read by the sexual networks (see sti.BaseNetwork.net_beta())

            # Natural history
            dur_exp = ss.constant(ss.weeks(1)), # Latent period: time from infection to being infectious
            p_symp = by_sex(0.375, 0.375), # Probability that an infection becomes symptomatic
            dur_presymp = lognorm_by_sex([ss.weeks(1), ss.weeks(12)], [ss.weeks(0.25), ss.weeks(3)]), # Time from being infectious to symptoms
            p_symp_care = by_sex(0.3, 0.2), # Probability that someone with symptoms seeks care
            dur_symp2care = lognorm_by_sex([ss.weeks(4), ss.weeks(4)], [ss.weeks(6), ss.weeks(4)]), # Time from symptoms to seeking care
            dur_symp2clear = lognorm_by_sex([ss.weeks(52), ss.weeks(5)], [ss.weeks(52), ss.weeks(5)]), # Time from symptoms to clearance
            dur_asymp2clear = lognorm_by_sex([ss.weeks(52), ss.weeks(5)], [ss.weeks(52), ss.weeks(5)]), # Time from being infectious to clearance, without symptoms
            p_clear = 1.0, # Probability that an asymptomatic infection in a woman clears (the rest persist)
            dur_persist = ss.constant(ss.years(100)), # Duration of a persistent infection
            p_pid = 0.2, # Probability that an infected woman develops pelvic inflammatory disease (PID)
            dur_prepid = ss.lognorm_ex(ss.weeks(6), ss.weeks(4)), # Time from being infectious to PID
            p_pid_care = 0.1, # Probability that a woman with PID seeks care
            dur_pid2care = ss.lognorm_ex(ss.weeks(2), ss.weeks(4)), # Time from PID to seeking care
            dur_pid2clear = ss.lognorm_ex(ss.weeks(52), ss.weeks(5)), # Time from PID to clearance
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.BoolState('exposed', label='Exposed'),
            ss.BoolState('infectious', label='Infectious'),
            ss.BoolState('presymptomatic', label='Presymptomatic'),
            ss.BoolState('asymptomatic', label='Asymptomatic'),
            ss.BoolState('symptomatic', label='Symptomatic'),
            ss.BoolState('pid', label='PID'),
            ss.BoolArr('will_symp', label='Symptoms to come'),
            ss.BoolArr('symptoms', label='Symptoms started this step'),
            ss.BoolArr('seek_care', label='Sought care this step, with symptoms'),
            ss.BoolArr('seek_care_pid', label='Sought care this step, with PID'),
            *[ss.TimeArr(f'ti_{key}') for key in ['infectious', 'symptoms', 'seek_care', 'clear', 'pid', 'seek_care_pid', 'clear_pid']], # When each happens, drawn at infection (NaN for never), due at the nearest step
            overwrite = True, # In place of ss.Infection's alias, infectious='infected'
        )
        self.u = sc.objdict({key:ss.random(name=f'u_{key}') for key in self.durs + self.chances}) # A uniform for each draw of the course, so each person's is their own
        self.set_strata(age_range if age_range is not None else [15, 65], age_bins, sex_keys)
        self.include_care = include_care
        return

    @property
    def treatable(self):
        """ Who a treatment can cure: everyone infected """
        return self.infected

    def validate_beta(self):
        """ The betas of the sexual networks and the maternal one, from the parameters (so a fitted beta_m2f is read as it is now) """
        p = self.pars
        self.pars.beta = sexual_beta(p.beta_m2f, p.beta_m2f*p.rel_beta_f2m, p.beta_m2m) | dict(maternalnet=ss.permonth(p.beta_m2c))
        return super().validate_beta()

    def sexed(self, val, uids):
        """ A value for some agents: for each one's sex if it is by sex (`sti.by_sex()`), or as it is; a duration in steps """
        step = lambda v: v/self.t.dt if isinstance(v, ss.dur) else v
        if isinstance(val, ss.by):
            sex = val.given['sex']
            return np.where(self.sim.people.female.raw[uids], step(sex['female']), step(sex['male']))
        return step(val)

    def draw(self, key, uids):
        """ A draw of the course for some agents: a chance (True with the probability) or a duration in steps (a lognormal by sex, or any distribution) """
        val = self.pars[key]
        if key in self.chances:
            return self.u[key].rvs(uids) < self.sexed(val, uids)
        if isinstance(val, ss.lognorm_ex) and isinstance(val.pars.mean, ss.by): # By sex: as the lognormal draws, from each one's mean and standard deviation
            mean, std = self.sexed(val.pars.mean, uids), self.sexed(val.pars.std, uids)
            sigma = np.sqrt(np.log1p((std/mean)**2))
            return np.exp(np.log(mean) - sigma**2/2 + sigma*sps.ndtri(self.u[key].rvs(uids)))
        return val.rvs(uids) if isinstance(val, ss.Dist) else np.full(len(uids), self.sexed(val, uids), dtype=float)

    def set_prognoses(self, uids, sources=None):
        """ Infect some people, and draw their course: when each stage starts, and whether it does """
        super().set_prognoses(uids, sources)
        ti = self.ti
        female = self.sim.people.female.raw
        self.infected[uids] = True
        self.exposed[uids] = True
        self.ti_infected[uids] = ti if self.sim.initialized else 0 # The infections seeded at the start are of its first point, though their course runs from the first step, as ss.Infection.init_post() dates it
        start = self.ti + 1 if self.sim.initialized else self.ti # A course runs from when the infection happens: at the end of the step, in transmission (as a declared disease dates it), or as ss.Infection.init_post() dates those seeded at the start
        infectious = start + self.draw('dur_exp', uids)
        self.ti_infectious[uids] = infectious

        # Symptoms, and seeking care for them; or none, and clearing (or, for some women, persisting)
        symp = self.draw('p_symp', uids)
        self.will_symp[uids] = symp
        s, a = uids[symp], uids[~symp]
        symptoms = infectious[symp] + self.draw('dur_presymp', s)
        self.ti_symptoms[s] = symptoms
        self.ti_seek_care[s] = np.where(self.draw('p_symp_care', s), symptoms + self.draw('dur_symp2care', s), np.nan)
        self.ti_clear[s] = symptoms + self.draw('dur_symp2clear', s)
        persist = female[a] & ~self.draw('p_clear', a)
        self.ti_clear[a] = infectious[~symp] + np.where(persist, self.draw('dur_persist', a), self.draw('dur_asymp2clear', a))

        # PID, for some women, and seeking care for it
        women = female[uids]
        w = uids[women]
        pid = self.draw('p_pid', w)
        w, pid_start = w[pid], infectious[women][pid] + self.draw('dur_prepid', w[pid])
        self.ti_pid[w] = pid_start
        self.ti_seek_care_pid[w] = np.where(self.draw('p_pid_care', w), pid_start + self.draw('dur_pid2care', w), np.nan)
        self.ti_clear_pid[w] = pid_start + self.draw('dur_pid2clear', w)
        return

    def step_state(self):
        """ Take the people whose time has come through their course, in its order """
        due = lambda state, when: (state & when.due()).uids
        for arr in [self.symptoms, self.seek_care, self.seek_care_pid]: # Events of the step only
            arr[:] = False

        uids = due(self.exposed, self.ti_infectious) # Infectious
        self.exposed[uids] = False
        self.infectious[uids] = True
        symp = self.will_symp[uids]
        self.presymptomatic[uids[symp]] = True
        self.asymptomatic[uids[~symp]] = True

        uids = due(self.presymptomatic, self.ti_symptoms) # Symptoms
        self.presymptomatic[uids] = False
        self.symptomatic[uids] = True
        self.symptoms[uids] = True

        uids = due(self.infectious & ~self.pid, self.ti_pid) # PID, from any other infectious stage
        for state in [self.presymptomatic, self.asymptomatic, self.symptomatic]:
            state[uids] = False
        self.pid[uids] = True

        for state, when, event in [(self.symptomatic, self.ti_seek_care, self.seek_care), (self.pid, self.ti_seek_care_pid, self.seek_care_pid)]: # Seeking care, once
            uids = due(state, when)
            event[uids] = True
            when[uids] = np.nan

        cleared = due(self.infectious & ~self.pid, self.ti_clear) | due(self.pid, self.ti_clear_pid)
        self.clear_infection(cleared)
        return

    def clear_infection(self, uids):
        """ Cure people (on their own, or treated): they are susceptible again """
        for state in [self.infected, self.exposed, self.infectious, self.presymptomatic, self.asymptomatic, self.symptomatic, self.pid, self.will_symp]:
            state[uids] = False
        self.susceptible[uids] = True
        for key in ['infectious', 'symptoms', 'seek_care', 'clear', 'pid', 'seek_care_pid', 'clear_pid']:
            self[f'ti_{key}'][uids] = np.nan
        return

    def step(self):
        """ Transmission (see `ss.Infection.step()`), and the new infections by route """
        new, sources, networks = super().step()
        maternal = np.isin(networks, [i for i,net in enumerate(self.sim.networks.values()) if isinstance(net, ss.PrenatalNet)])
        self.sexual_sources = sources[~maternal]
        ti = self.ti
        self.results.new_infections_sex[ti] = np.count_nonzero(~maternal)
        self.results.new_infections_prenatal[ti] = self.results.new_infections_mtct[ti] = np.count_nonzero(maternal)
        return new, sources, networks

    def new_infected(self):
        """ Who was infected in this step, including those seeded at the start, as a boolean array """
        return self.ti_infected == self.ti

    def init_results(self):
        """ Initialize results """
        super().init_results()
        self.define_sti_results(include_care=self.include_care)

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
        self.define_results(*[res for res in results if res.name not in self.results])
        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        ppl = self.sim.people
        res = self.results
        t = self.tally()
        self.update_sti_results(t)

        # Symptoms and care seeking, by sex and age
        new_sym = t.vals(self.symptoms)
        seekers = t.vals(self.seek_care) | t.vals(self.seek_care_pid)
        symptomatic = t.vals(self.symptomatic)
        for skk, skl, attr in self.sex_strata():
            sex = t.vals(ppl[attr])
            new_sym_s = new_sym & sex
            sym_s = symptomatic & sex
            res[f'symp_prevalence{skk}'][ti] = t.div(t.count(sym_s & t.adults), t.count(t.adults & sex))
            res[f'new_symptomatic{skk}'][ti] = t.count(new_sym_s)
            if skk:
                res[f'symptomatic{skk}'][ti] = t.count(sym_s)
            if self.include_care:
                res[f'new_care_seekers{skk}'][ti] = t.count(seekers & sex)
            if self.age_bins is not None:
                self.set_age_results(skk, new_symptomatic=t.hist(new_sym_s), symp_prevalence=t.hdiv(t.hist(sym_s), t.hist(sex)))
        return
