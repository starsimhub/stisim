"""
Bacterial STIs as declared Starsim diseases: the base class, and the SEIS template
used for chlamydia, gonorrhea, and trichomoniasis
"""

import numpy as np
import sciris as sc
import starsim as ss

__all__ = ['BaseSTI', 'SEIS', 'sexed', 'default_age_bins', 'default_sex_keys']

# Default stratifications of results: sti.Sim stratifies by sex and by these age groups
default_age_bins = [0, 15, 20, 25, 30, 35, 50, 65, 100]
default_sex_keys = {'': 'alive', 'f': 'female', 'm': 'male'}
female = ss.people.female


def sexed(val):
    """
    Convert a value by sex as STIsim has written it to a Starsim value

    `[women, men]` becomes `ss.by(sex=[women, men])`, `[[mean, std], [mean, std]]` a lognormal
    with a mean and standard deviation by sex, and `ss.bernoulli(p=x)` its probability `x`.
    Anything else is returned unchanged.

    Examples:
        ```python
        sti.sexed([0.35, 0.65]) # ss.by(sex=[0.35, 0.65])
        sti.sexed([[ss.months(8), ss.months(2)], [ss.months(6), ss.months(3)]]) # ss.lognorm_ex by sex
        ```
    """
    if isinstance(val, ss.bernoulli):
        return val.pars.p
    if isinstance(val, (list, tuple)) and len(val) == 2:
        f, m = val
        if all(isinstance(v, (list, tuple)) for v in val):
            return ss.lognorm_ex(mean=ss.by(sex=[f[0], m[0]]), std=ss.by(sex=[f[1], m[1]]))
        return ss.by(sex=list(val))
    return val


def init_prev_fn(data):
    """
    Initial prevalence by risk group, sex, and sex work, as a function of `(sim, uids)` for `init_prev=` (agents only)

    Args:
        data (dataframe): columns `risk_group`, `sex` (`'female'` or `'male'`), `sw` (0 or 1), and `init_prev`; or a number
    """
    if sc.isnumber(data):
        return data
    def init_prev(sim, uids):
        nw = sim.networks.structuredsexual
        ppl = sim.people
        out = np.zeros(len(uids))
        active = nw.active(ppl)[uids]
        for row in data.sort_values('sw').itertuples(): # As in STIsim 3: by group, then sex workers and clients overwrite theirs
            cond = ppl[row.sex][uids] & (nw.risk_group[uids] == row.risk_group) & active
            if row.sw:
                cond &= (nw.fsw if row.sex == 'female' else nw.client)[uids]
            out[np.asarray(cond)] = row.init_prev
        return np.clip(out, 0, 1)
    return init_prev


class BaseSTI(ss.Disease):
    """
    Base class for sexually transmitted infections declared as states and transitions

    Transmission is given per network, as in STIsim 3: `beta_m2f` (and `rel_beta_f2m`) on
    `structuredsexual`, `beta_m2c` on `maternal`, and `beta_m2m` on `msm`; other networks
    do not transmit. These are turned into the beta of each network when the sim is
    initialized. Give `beta` directly (a number, a rate, or a dict by network) to bypass
    this, e.g. `beta=ss.peryear(2)` for homogeneous mixing with no networks.

    Subclasses declare their states and transitions in the class body; `defaults` gives
    the values of the parameters that differ from the parent's, in STIsim's forms (see
    `sti.sexed()`).

    Args:
        pars (dict): parameter values
        name (str): the name of the module
        init_prev (float): the initial prevalence (or `ss.bernoulli(p=...)`)
        init_prev_data (dataframe): initial prevalence by risk group, sex, and sex work (columns `risk_group`, `sex`, `sw`, `init_prev`), among the sexually active in `structuredsexual`; instead of init_prev
        age_range (list): the ages of the adults that `prevalence` and `symp_prevalence` count (default [15, 65])
        age_bins (list): the age groups `sti.Sim` stratifies by, if this disease is in it; None for none
        sex_keys (dict): the sexes of the flat results (`prevalence_f`, ...); None for none
        kwargs (dict): parameter values by name, or `label`, `dt`
    """
    beta = ss.par(0.0, doc='Beta by network; worked out from beta_m2f, beta_m2c, and beta_m2m unless given')
    beta_m2f = ss.par(None, doc='Per-act probability of transmission, male to female')
    rel_beta_f2m = ss.par(0.5, doc='Relative per-act probability of transmission, female to male')
    beta_m2c = ss.par(None, doc='Monthly probability of mother-to-child transmission')
    beta_m2m = ss.par(None, doc='Per-act probability of transmission, male to male')
    eff_condom = ss.par(1.0, doc='Efficacy of condoms')
    defaults = {}
    transmits = True # False for a condition that is not transmitted (BV): it needs no networks

    def __init__(self, pars=None, name=None, init_prev=None, init_prev_data=None, age_range=None,
                 age_bins=default_age_bins, sex_keys=default_sex_keys, **kwargs):
        names = ss.transitions.class_members(type(self)).keys()
        if init_prev_data is not None:
            init_prev = init_prev_fn(init_prev_data)
        given = sc.mergedicts(pars, {k:kwargs.pop(k) for k in list(kwargs) if k in names}, None if init_prev is None else dict(init_prev=init_prev))
        if not self.transmits:
            if given.get('init_prev') is not None:
                errormsg = f'{type(self).__name__} is not transmitted, so has no init_prev; set the shares in its states instead'
                raise ValueError(errormsg)
            given.pop('init_prev', None)
            pars = sc.mergedicts(self.defaults, given)
            super().__init__(name=name, **kwargs)
        else:
            pars = sc.mergedicts(dict(init_prev=0.01), self.defaults, given)
            init_prev = sexed(pars.pop('init_prev'))
            super().__init__(name=name, init_prev=ss.par(init_prev, name='init_prev'), **kwargs)
        self.set_pars({k:sexed(v) for k,v in pars.items()})
        self.beta_given = 'beta' in given
        self.age_range = age_range if age_range is not None else [15, 65]
        self.age_bins = None if age_bins is None else list(age_bins)
        self.sex_keys = {'': 'alive'} if sex_keys is None else dict(sex_keys)
        return

    def value(self, key):
        """ The plain value of a parameter """
        val = self.pars[key]
        return val.value if isinstance(val, ss.par) else val

    def betamap(self, sim):
        """ The beta of each network in the sim, from beta_m2f, rel_beta_f2m, beta_m2c, and beta_m2m """
        m2f, rel, m2c, m2m = [self.value(k) for k in ['beta_m2f', 'rel_beta_f2m', 'beta_m2c', 'beta_m2m']]
        out = {}
        for name in sim.networks.keys():
            key = ss.standardize_netkey(name)
            if key == 'structuredsexual' and m2f is not None:
                out[name] = [m2f, m2f*rel]
            elif key == 'msm' and m2m is not None:
                out[name] = [m2m, m2m]
            elif key == 'maternal' and m2c:
                out[name] = [ss.permonth(m2c), 0]
            else:
                out[name] = 0
        return out

    def modifiers(self, sim):
        """ Modifiers that are part of the disease (e.g. syphilis latent infectiousness, BV risk factors from the network): added to the sim's connectors """
        return []

    def add_modifiers(self, sim):
        """ Add the disease's own modifiers to the sim, if they are not there already. TODO(v4-rs): a way for a module to declare its own modifiers """
        for mod in self.modifiers(sim):
            if mod.name not in sim.connectors:
                sim.connectors[mod.name] = mod
                mod.init_pre(sim)
        return

    def init_pre(self, sim):
        if self.transmits and not self.beta_given:
            if not len(sim.networks):
                errormsg = f'"{self.name}" has no networks to transmit over; give beta= (e.g. ss.peryear(2)) for homogeneous mixing'
                raise ValueError(errormsg)
            self.pars.beta.value = self.betamap(sim)
        super().init_pre(sim)
        self.add_modifiers(sim)
        return

    @property
    def strats(self):
        """ The sim's stratifications, by name """
        return {s.name:s for s in self.sim.strats}

    def adult_levels(self):
        """ The age groups within age_range, if the sim is stratified by age """
        age = self.strats.get('age')
        if age is None or age.edges is None:
            return None
        lo, hi = self.age_range
        edges = list(age.edges) + [np.inf]
        return [label for label,a,b in zip(age.labels, edges[:-1], edges[1:]) if a >= lo and b <= hi]

    def flat_results(self):
        """ Results by sex as STIsim 3 named them, among adults where they are prevalences: (name, numerator, denominator or None, label) """
        return [
            ('new_infections', 'new_infections', None, 'New infections'),
            ('n_infected', 'n_infected', None, 'Number infected'),
            ('incidence', 'incidence', None, 'Incidence'),
        ]

    def finalize_results(self):
        """ Add results by sex (`prevalence_f`, ...), and make `prevalence` that of adults, as in STIsim 3 """
        super().finalize_results()
        res = self.results
        strats = self.strats
        adults = self.adult_levels()
        sexes = [(k, v) for k,v in self.sex_keys.items() if k and 'sex' in strats]
        for key,num,den,label in self.flat_results():
            if getattr(res[num], 'strata', None) is None: # E.g. a ratio of results with different dimensions
                continue
            if adults is not None and den is not None: # A prevalence among adults
                res[key].values[:] = res[key].select(age=adults).values
            for sk,attr in sexes:
                level = 'female' if attr == 'female' else 'male'
                sel = dict(sex=level) if den is None or adults is None else dict(sex=level, age=adults)
                vals = res[num].select(**sel).values
                self.results += ss.Result(f'{key}_{sk}', values=vals, scale=False, label=f'{label} ({sk.upper()})', timevec=res[num].timevec, auto_plot=False, module=self.name)
        return


class SEIS(BaseSTI):
    """
    Template for bacterial STIs: an SEIS infection with symptoms, care seeking, and PID

    States (S susceptible, E exposed, P presymptomatic, A asymptomatic, Y symptomatic,
    C sought care, PID pelvic inflammatory disease):

        S -> E -> P (p_symp) | A
        P -> Y                                   after dur_presymp
        Y -> C (p_symp_care) | PID (women, p_pid), or -> S after dur_symp2clear
        A -> PID (women, p_pid) | S persistent (women, 1 - p_clear) | S after dur_asymp2clear
        PID -> C (p_pid_care), or -> S after dur_pid2clear
        C -> S                                   after dur_symp2clear
        (E | P | A | Y | C | PID) -> S           cure: moved only by treatment

    Values by sex are `ss.by(sex=[women, men])`; STIsim's forms (`[women, men]`, and
    `[[mean, std], [mean, std]]` for durations) are converted by `sti.sexed()`. Seeking care
    is entering C (`dis.C.enter`), the event that care-seeking interventions are triggered
    by. Treatment moves people along `cure`; failures stay where they are.

    Results: the counts and flows of the states and transitions; `new_infections`,
    `n_infected`, `prevalence` (among adults in `age_range` when the sim is stratified by
    age), `incidence`, `new_symptomatic`, `n_symptomatic`, `symp_prevalence`, and
    `new_care_seekers`, with `.select(sex=..., age=...)` and `.by()` when the sim is
    stratified, and flat `_f` and `_m` versions of each.
    """
    S, E, P, A, Y, C, PID = ss.states('S E P A Y C PID')

    # Natural history
    dur_exp = ss.par(ss.constant(ss.weeks(1)), doc='Time from exposure to being infectious')
    p_symp = ss.par(sexed([0.375, 0.375]), doc='Probability of symptoms')
    dur_presymp = ss.par(sexed([[ss.weeks(1), ss.weeks(12)], [ss.weeks(0.25), ss.weeks(3)]]), doc='Time to symptoms, for those who develop them')
    p_symp_care = ss.par(sexed([0.3, 0.2]), doc='Probability that the symptomatic seek care')
    dur_symp2care = ss.par(sexed([[ss.weeks(4), ss.weeks(4)], [ss.weeks(6), ss.weeks(4)]]), doc='Time from symptoms to seeking care')
    p_pid = ss.par(0.2, doc='Probability of PID (women)')
    dur_prepid = ss.par(ss.lognorm_ex(ss.weeks(6), ss.weeks(4)), doc='Time to PID')
    p_pid_care = ss.par(0.1, doc='Probability of seeking care for PID')
    dur_pid2care = ss.par(ss.lognorm_ex(ss.weeks(2), ss.weeks(4)), doc='Time from PID to seeking care')
    dur_asymp2clear = ss.par(sexed([[ss.weeks(52), ss.weeks(5)], [ss.weeks(52), ss.weeks(5)]]), doc='Duration of untreated asymptomatic infection')
    dur_symp2clear = ss.par(sexed([[ss.weeks(52), ss.weeks(5)], [ss.weeks(52), ss.weeks(5)]]), doc='Duration of untreated symptomatic infection')
    dur_pid2clear = ss.par(ss.lognorm_ex(ss.weeks(52), ss.weeks(5)), doc='Duration of untreated PID')
    p_clear = ss.par(1.0, doc='Probability that an asymptomatic infection in a woman clears spontaneously (the rest persist)')
    dur_persist = ss.par(ss.constant(ss.years(100)), doc='Duration of a persistent infection')

    # Transitions
    infection = S >> E @ BaseSTI.beta*(P + A + Y + C + PID)
    latency   = ss.Transition(E, [P, A], p=[p_symp, None], dur=dur_exp)
    onset     = ss.Transition(P, Y, dur=dur_presymp)
    seek      = ss.Transition(Y, C, p=p_symp_care, dur=dur_symp2care)
    pid_y     = ss.Transition(Y, PID, p=p_pid, dur=dur_prepid, where=female)
    clear_y   = ss.Transition(Y, S, dur=dur_symp2clear)
    pid_a     = ss.Transition(A, PID, p=p_pid, dur=dur_prepid, where=female)
    persist   = ss.Transition(A, S, p=(1 - p_clear)*(1 - p_pid), dur=dur_persist, where=female)
    clear_a   = ss.Transition(A, S, p=None, dur=dur_asymp2clear)
    seek_pid  = ss.Transition(PID, C, p=p_pid_care, dur=dur_pid2care)
    clear_pid = ss.Transition(PID, S, dur=dur_pid2clear)
    clear_c   = ss.Transition(C, S, dur=dur_symp2clear)
    cure      = ss.Transition(E | P | A | Y | C | PID, S, undriven=True)

    # Results
    n_infected       = E + P + A + Y + C + PID
    n_symptomatic    = Y + C
    prevalence       = (E + P + A + Y + C + PID)/(S + E + P + A + Y + C + PID)
    symp_prevalence  = (Y + C)/(S + E + P + A + Y + C + PID)
    incidence        = infection/S
    new_symptomatic  = ss.outcome(onset)
    new_care_seekers = seek + seek_pid

    def flat_results(self):
        return super().flat_results() + [
            ('prevalence', 'prevalence', True, 'Prevalence'),
            ('new_symptomatic', 'new_symptomatic', None, 'New symptomatic'),
            ('n_symptomatic', 'n_symptomatic', None, 'Number symptomatic'),
            ('symp_prevalence', 'symp_prevalence', True, 'Symptomatic prevalence'),
            ('new_care_seekers', 'new_care_seekers', None, 'New care seekers'),
        ]

    # States as STIsim 3 named them, as agents, for connectors and custom code
    @property
    def susceptible(self):  return self.S
    @property
    def exposed(self):      return self.E
    @property
    def infectious(self):   return self.P | self.A | self.Y | self.C | self.PID
    @property
    def infected(self):     return self.E | self.infectious
    @property
    def symptomatic(self):  return self.Y | self.C
    @property
    def asymptomatic(self): return self.E | self.P | self.A
    @property
    def pid(self):          return self.PID
    @property
    def seeking_care(self): # Entered C since this point in the last step
        out = self.C.asnew(np.zeros(len(self.C.raw), dtype=bool))
        out[self.graph.states.C.enter.uids] = True
        return out
    @property
    def treatable(self):    return self.infected
