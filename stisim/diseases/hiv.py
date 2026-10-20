"""
HIV: natural history and transmission, as a declared Starsim disease

Treatment is a module of its own, `sti.ART`, which changes HIV through modifiers; HIV knows nothing about it.
"""

import numpy as np
import sciris as sc
import starsim as ss
from stisim.diseases.sti import BaseSTI, default_age_bins, default_sex_keys

__all__ = ['HIV']

# Parameters that belong to ART, PrEP, or VMMC, and the module that has each
moved_pars = dict(
    art_cd4_growth='ART', effective_art_efficacy='ART', nonsupp_art_efficacy='ART', time_to_art_efficacy='ART',
    p_effective_art='ART', art_cd4_pars='ART', dur_on_art='ART', rel_dur_on_art='ART', dur_on_art_trend='ART',
    dur_post_art='ART', dur_post_art_scale_factor='ART', rel_art_mortality_effective='ART',
    rel_art_mortality_unsupp_m='ART', rel_art_mortality_unsupp_f='ART', art_death_age='ART', art_death_dur='ART',
    care_seeking='ART', maternal_care_scale='ART', eff_circ='VMMC', eff_circ_traditional='VMMC',
)
dropped_pars = ['p_hiv_death', 'dist_ti_init_infected', 'beta', 'include_care']

# CD4 bands of the death rate, as STIsim 3 gave them: descending edges, and the annual rate in each
default_cd4_death_bins = [1000, 500, 350, 200, 50, 0]
default_cd4_death_rates = [0.003, 0.003, 0.005, 0.01, 0.05, 0.30]


def cd4_bands(bins, rates):
    """ STIsim 3's CD4 death bands (descending edges, read with np.digitize) as ascending edges and a rate for each band """
    bins = np.asarray(bins, dtype=float)
    rates = np.asarray(rates, dtype=float)
    if len(bins) != len(rates) or np.any(np.diff(bins) >= 0):
        errormsg = f'cd4_death_bins must be descending, with one rate in cd4_death_rates for each (CD4 >= the first bin uses the first rate); got {bins} and {rates}'
        raise ValueError(errormsg)
    edges = list(bins[::-1]) # E.g. [0, 50, 200, 350, 500, 1000]: the band [0, 50) has the last rate, and >= 1000 the first
    return edges, list(rates[::-1])


class HIV(BaseSTI):
    """
    HIV: the natural history without treatment, and transmission

    States:

        S -> acute               sexual transmission (infection), and from mother to child in pregnancy (prenatal) and breastfeeding (postnatal)
        acute -> latent          latency, after dur_acute
        latent -> falling        progression, after dur_latent (CD4 falls)
        falling -> dead          aids_death, after dur_falling (CD4 reaches 0)
        acute | latent | falling -> dead   cd4_death: a rate by CD4 band (cd4_death_bins, cd4_death_rates), x rel_death, x rel_death_f for women

    The CD4 count is an attribute (`hiv.cd4`): drawn at infection (`cd4_start`, then `cd4_latent`), falling through
    acute infection to `cd4_latent`, constant while latent, and falling linearly to 0 over the falling stage.
    Infectiousness is `rel_trans_acute` in acute infection and `rel_trans_falling` with CD4 below 200, each drawn
    for each person at infection. `rel_sus_age` gives relative susceptibility by age and sex.

    Treatment, prevention, and testing are modules of their own: `sti.ART` (which changes CD4, infectiousness,
    progression, mortality, and mother-to-child transmission through modifiers), `sti.Prep`, `sti.VMMC`, and
    `sti.HIVTest`. HIV keeps who has been diagnosed (`diagnosed`, `ti_diagnosed`), which the tests set.

    Results: the counts and flows of the states and transitions; `new_infections` (with `new_infections_sex`,
    `new_infections_mtct`, `new_infections_prenatal`, `new_infections_postnatal`), `n_infected`, `prevalence`
    (of everyone), `prevalence_15_49`, `incidence`, `new_deaths` and `cum_deaths` (with ART's deaths after
    stopping, if `sti.ART` is in the sim), `new_diagnoses`, `cum_diagnoses`, `n_diagnosed`; by sex (`_f`,
    `_m`, among adults for prevalence) and age; and with `sti.StructuredSexual`, by sex work and risk group.

    Args:
        init_prev (float): initial prevalence, back-dated over `init_window` so people start part way through their natural history
        init_prev_data (dataframe): initial prevalence by risk group, sex, and sex work (see `sti.BaseSTI`), also back-dated
        rel_init_prev (float): a multiplier on the initial prevalence
        init_window (tuple): how long before the start the initial infections happened (default 5 months to 10 years, as in STIsim 3)
        cd4_death_bins (list): descending CD4 edges of the death rate (STIsim 3's form)
        cd4_death_rates (list): the annual death rate in each band
        kwargs (dict): parameter values by name, and those of `sti.BaseSTI`

    Examples:
        ```python
        hiv = sti.HIV(beta_m2f=0.05, init_prev=0.1)
        sim = sti.Sim(diseases=hiv, interventions=[sti.HIVTest(test_prob_data=0.2), sti.ART(coverage=0.8)])
        ```
    """
    S, acute, latent, falling = ss.states('S acute latent falling')

    # Natural history without treatment
    dur_acute = ss.par(ss.lognorm_ex(ss.months(1.7), ss.months(1)), doc='Duration of acute infection (Bellan 2015: central estimate 1.7 months)')
    dur_latent = ss.par(ss.lognorm_ex(ss.years(10), ss.years(3)), doc='Duration of latent infection, untreated')
    dur_falling = ss.par(ss.lognorm_ex(ss.years(3), ss.years(1)), doc='Duration of late-stage infection, as CD4 falls to 0 (then AIDS death)')
    cd4_start = ss.par(ss.normal(800, 50), doc='CD4 count before infection')
    cd4_latent = ss.par(ss.normal(500, 50), doc='CD4 count during latent infection')
    rel_death = ss.par(1.0, doc='Multiplier on the CD4-banded death rate')
    rel_death_f = ss.par(0.74, doc='Multiplier on the CD4-banded death rate, for women')
    cd4_death_rate = ss.par(ss.by(cd4=ss.peryear(default_cd4_death_rates[::-1])), doc='Death rate by CD4 band (from cd4_death_bins and cd4_death_rates)')
    sex_death = ss.par(ss.by(sex=dict(female=rel_death_f, male=1.0)), doc='Death rate multiplier by sex (from rel_death_f)')

    # Transmission
    beta_breastfeed = ss.par(0.005, doc='Monthly probability of transmission to a breastfed child (~14% over 12 months without ART)')
    rel_trans_acute = ss.par(ss.normal(5.3, 0.5), doc='Relative infectiousness in acute infection (Bellan 2015: RH=5.3), drawn at infection')
    rel_trans_falling = ss.par(ss.normal(8, 0.5), doc='Relative infectiousness with CD4 below 200, drawn at infection')
    sexual_beta = ss.par(dict(
        structuredsexual = ss.peract(m2f=BaseSTI.beta_m2f, f2m=BaseSTI.beta_m2f*BaseSTI.rel_beta_f2m),
        mfnetwork = ss.peract(m2f=BaseSTI.beta_m2f, f2m=BaseSTI.beta_m2f*BaseSTI.rel_beta_f2m),
        msm = ss.peract(BaseSTI.beta_m2m),
        default = 0,
    ), doc='Sexual transmission by network, from beta_m2f, rel_beta_f2m, and beta_m2m')
    prenatal_beta = ss.par(dict(maternalnet=BaseSTI.beta_m2c*ss.permonth(1), default=0), doc='Transmission in pregnancy, from beta_m2c')
    postnatal_beta = ss.par(dict(breastfeedingnet=beta_breastfeed*ss.permonth(1), default=0), doc='Transmission in breastfeeding, from beta_breastfeed')

    # Transitions
    infection = S >> acute @ sexual_beta*(acute + latent + falling)
    prenatal = S >> acute @ prenatal_beta*(acute + latent + falling)
    postnatal = S >> acute @ postnatal_beta*(acute + latent + falling)
    latency = ss.Transition(acute, latent, dur=dur_acute)
    progression = ss.Transition(latent, falling, dur=dur_latent)
    aids_death = ss.Transition(falling, ss.null, dur=dur_falling)

    # CD4 and infectiousness, for each person
    acquired = infection | prenatal | postnatal
    cd4_at_infection = ss.attr(cd4_start, on=acquired)
    cd4_set_point = ss.attr(cd4_latent, on=acquired)
    rta = ss.attr(rel_trans_acute, on=acquired)
    rtf = ss.attr(rel_trans_falling, on=acquired)
    cd4 = ss.attr({
        acute:   cd4_at_infection - (cd4_at_infection - cd4_set_point)*acute.progress,
        latent:  cd4_set_point,
        falling: np.maximum(1, cd4_set_point*(1 - falling.progress)),
    })
    cd4_death = (acute + latent + falling) >> ss.null @ cd4_death_rate*sex_death*rel_death

    # Results
    n_infected = acute + latent + falling
    prevalence = (acute + latent + falling)/(S + acute + latent + falling)
    new_infections = infection + prenatal + postnatal
    new_infections_sex = ss.outcome(infection)
    new_infections_prenatal = ss.outcome(prenatal)
    new_infections_postnatal = ss.outcome(postnatal)
    new_infections_mtct = prenatal + postnatal

    defaults = dict(beta_m2f=0.05, beta_m2c=0.025, eff_condom=0.9, init_prev=0.05)
    options = ['cd4_death_bins', 'cd4_death_rates', 'include_aids_deaths', 'init_diagnosed', 'init_window', 'rel_init_prev', 'rel_sus_age'] # Constructor arguments documented with the parameters

    def __init__(self, pars=None, name='hiv', init_prev_data=None, rel_init_prev=1.0, init_window=None,
                 age_bins=default_age_bins, sex_keys=default_sex_keys, **kwargs):
        pars = sc.mergedicts(pars)
        for src in [pars, kwargs]:
            for key in list(src):
                if key in moved_pars:
                    errormsg = f'"{key}" is a parameter of sti.{moved_pars[key]}, not sti.HIV (see docs/design/hiv_art_split.md); give it there, e.g. sti.{moved_pars[key]}({key}=...)'
                    raise ValueError(errormsg)
                if key in dropped_pars:
                    src.pop(key) # Starsim 3 parameters with no STIsim 4 equivalent
        if 'cd4_death_bins' in pars or 'cd4_death_rates' in pars or 'cd4_death_bins' in kwargs or 'cd4_death_rates' in kwargs:
            bins = pars.pop('cd4_death_bins', kwargs.pop('cd4_death_bins', default_cd4_death_bins))
            rates = pars.pop('cd4_death_rates', kwargs.pop('cd4_death_rates', default_cd4_death_rates))
        else:
            bins, rates = default_cd4_death_bins, default_cd4_death_rates
        self.cd4_edges, cd4_rates = cd4_bands(bins, rates)
        pars['cd4_death_rate'] = ss.by(cd4=ss.peryear(cd4_rates))
        for key in ['init_diagnosed', 'rel_sus_age', 'rel_init_prev', 'include_aids_deaths']: # Options rather than parameters of the graph
            for src in [pars, kwargs]:
                if key in src:
                    if key == 'rel_init_prev':
                        rel_init_prev = src.pop(key)
                    else:
                        setattr(self, f'_{key}', src.pop(key))
        init_diagnosed = getattr(self, '_init_diagnosed', 0)
        self.init_diagnosed = init_diagnosed.pars.p if isinstance(init_diagnosed, ss.bernoulli) else init_diagnosed
        self.rel_sus_age = getattr(self, '_rel_sus_age', None)
        self.include_aids_deaths = getattr(self, '_include_aids_deaths', True)
        self.rel_init_prev = rel_init_prev
        self.init_window = init_window if init_window is not None else (ss.months(5), ss.years(10))

        # Initial prevalence: back-dated, so people start part way through their natural history
        init_prev = pars.pop('init_prev', kwargs.pop('init_prev', None))
        if isinstance(init_prev, ss.bernoulli):
            init_prev = init_prev.pars.p
        self.init_data = None
        if init_prev_data is not None:
            if sc.isnumber(init_prev_data):
                init_prev = init_prev_data
            else:
                self.init_data = init_prev_data.copy()
                self.init_data['init_prev'] = np.clip(self.init_data['init_prev']*rel_init_prev, 0, 1)
                init_prev = 0 # Seeded by burnin_spread()
        if init_prev is None:
            init_prev = self.defaults['init_prev']
        if self.init_data is None and init_prev:
            init_prev = ss.spread(np.clip(init_prev*rel_init_prev, 0, 1), over=self.init_window)
        super().__init__(name=name, pars=pars, init_prev=init_prev, age_bins=age_bins, sex_keys=sex_keys, **kwargs)
        self.age_range = [15, 50]
        self.stratify(cd4=self.cd4_edges, states=self.graph.states.acute + self.graph.states.latent + self.graph.states.falling) # The CD4 bands of the death rate (derived: no groups)
        self.other_deaths = [] # (module, result) of the deaths that other modules count as HIV deaths (see finalize_results())
        self.define_states(
            ss.BoolArr('diagnosed', label='Diagnosed'), # Not a BoolState: its automatic result is not counted in a declared disease (see n_diagnosed)
            ss.FloatArr('ti_diagnosed', label='Time of diagnosis'),
        )
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if not self.include_aids_deaths:
            self.add_modifiers_from([ss.modify(self.aids_death, multiply=0, name=f'{self.name}_no_aids_deaths')], sim)
        return

    def add_modifiers_from(self, mods, sim):
        for mod in mods:
            if mod.name not in sim.connectors:
                sim.connectors[mod.name] = mod
                mod.init_pre(sim)
        return

    def modifiers(self, sim):
        """ Infectiousness by stage and CD4, and susceptibility by age and sex """
        st, name = self.graph.states, self.name
        rt = f'{name}.rel_trans'
        mods = [
            ss.modify(rt, multiply=self.graph.attrs.rta, where=st.acute, name=f'{name}_rel_trans_acute'),
            ss.modify(rt, multiply=self.graph.attrs.rtf, where=(st.acute | st.latent | st.falling) & (self.graph.attrs.cd4 < 200), name=f'{name}_rel_trans_aids'),
        ]
        for i,(lo, hi, sex, mult) in enumerate(self.rel_sus_age or []):
            who = (ss.ppl.age >= lo) & (ss.ppl.age < hi)
            if sex == 'f':
                who = who & ss.ppl.female
            elif sex == 'm':
                who = who & ss.ppl.male
            mods.append(ss.modify(f'{name}.rel_sus', multiply=mult, where=who, name=f'{name}_rel_sus_age{i}'))
        return mods

    def burnin_spread(self, engine):
        """ With init_prev_data: each person's chance of being infected at the start, back-dated over init_window (see ss.spread) """
        if self.init_data is None:
            return None
        from stisim.diseases.sti import init_prev_fn
        fn = init_prev_fn(self.init_data)
        sim = engine.sim
        lo, hi = [ss.years(v).value if isinstance(v, ss.TimePar) else v for v in self.init_window] # In years
        steps = max(1, int(round(hi/sim.t.dt.years))) # UPSTREAM: a per-person ss.spread is only possible through this hook, which back-dates over [0, steps), not over a window
        def share():
            out = np.zeros(len(self.rel_sus.raw))
            uids = sim.people.auids
            out[uids] = fn(sim, uids)
            return out
        note = f'{self.name}: initial infections by risk group, sex, and sex work, back-dated over the {steps} steps before the start'
        return sc.objdict(spread=None, share=share, steps=steps, entry=None, tr=self.graph.transitions.infection, target=['acute', 'latent', 'falling'], note=note, after=None)

    def init_post(self):
        super().init_post()
        if self.init_diagnosed:
            infected = self.infected.uids
            dx = infected[np.random.default_rng(int(self.sim.pars.rand_seed or 0) + 11).random(len(infected)) < self.init_diagnosed]
            self.diagnosed[dx] = True
            self.ti_diagnosed[dx] = 0
        return

    def init_results(self):
        super().init_results()
        results = [
            ss.Result('new_deaths', dtype=float, label='HIV deaths', scale=False),
            ss.Result('cum_deaths', dtype=float, label='Cumulative HIV deaths', scale=False, auto_plot=False),
            ss.Result('new_diagnoses', dtype=float, label='Diagnoses', scale=False),
            ss.Result('cum_diagnoses', dtype=float, label='Cumulative diagnoses', scale=False, auto_plot=False),
            ss.Result('n_diagnosed', dtype=float, label='Number diagnosed', scale=False, auto_plot=False),
            ss.Result('prevalence_15_49', dtype=float, label='Prevalence 15-49', scale=False),
            ss.Result('p_diagnosed_pregnant', dtype=float, label='Proportion of HIV+ pregnant women diagnosed', scale=False, auto_plot=False),
        ]
        if 'structuredsexual' in self.sim.networks:
            results += [
                ss.Result('prevalence_sw', dtype=float, label='FSW prevalence', scale=False, auto_plot=False),
                ss.Result('new_infections_sw', dtype=float, scale=False, label='New infections - FSW', auto_plot=False),
                ss.Result('new_infections_not_sw', dtype=float, scale=False, label='New infections - Other F', auto_plot=False),
                ss.Result('prevalence_client', dtype=float, label='Client prevalence', scale=False, auto_plot=False),
                ss.Result('new_infections_client', dtype=float, scale=False, label='New infections - Clients', auto_plot=False),
                ss.Result('new_infections_not_client', dtype=float, scale=False, label='New infections - Other M', auto_plot=False),
            ]
            for rg in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    results += [
                        ss.Result(f'prevalence_risk_group_{rg}_{sex}', scale=False, auto_plot=False),
                        ss.Result(f'new_infections_risk_group_{rg}_{sex}', dtype=float, scale=False, auto_plot=False),
                    ]
        self.define_results(*results)
        return

    def update_results(self):
        """ The results that need agents: diagnoses, prevalence 15-49, and by sex work and risk group """
        super().update_results()
        ti = self.ti
        res = self.results
        ppl = self.sim.people
        w = ss.results.people_scale(self.sim) # A declared module's results are written in people
        alive = ppl.alive.raw
        infected = self.infected.raw & alive
        newly = (self.ti_acute.raw == ti) & alive # Infected in the step just taken (dated at its end)
        dx = self.diagnosed.raw & alive
        res.new_diagnoses[ti] = w*np.count_nonzero((self.ti_diagnosed.raw == ti - 1) & alive) # Dated at the step's start by the tests
        res.cum_diagnoses[ti] = res.new_diagnoses[:ti+1].sum()
        res.n_diagnosed[ti] = w*np.count_nonzero(dx)
        age = ppl.age.raw
        adults = alive & (age >= 15) & (age < 50)
        res.prevalence_15_49[ti] = sc.safedivide(np.count_nonzero(infected & adults), np.count_nonzero(adults))
        preg = self.sim.demographics.get('pregnancy') if hasattr(self.sim, 'demographics') else None
        if preg is not None:
            pregnant = preg.pregnant.raw & alive
            res.p_diagnosed_pregnant[ti] = sc.safedivide(np.count_nonzero(dx & infected & pregnant), np.count_nonzero(infected & pregnant))
        if 'structuredsexual' in self.sim.networks:
            nw = self.sim.networks.structuredsexual
            fsw, client = nw.fsw.raw & alive, nw.client.raw & alive
            female = ppl.female.raw
            if fsw.any():
                res.prevalence_sw[ti] = np.count_nonzero(infected & fsw)/np.count_nonzero(fsw)
                res.new_infections_sw[ti] = w*np.count_nonzero(newly & fsw)
                res.new_infections_not_sw[ti] = w*np.count_nonzero(newly & ~fsw & female)
            if client.any():
                res.prevalence_client[ti] = np.count_nonzero(infected & client)/np.count_nonzero(client)
                res.new_infections_client[ti] = w*np.count_nonzero(newly & client)
                res.new_infections_not_client[ti] = w*np.count_nonzero(newly & ~client & ~female)
            rg = nw.risk_group.raw
            for g in range(nw.pars.n_risk_groups):
                for sex,sel in [('female', female), ('male', ~female)]:
                    grp = alive & sel & (rg == g)
                    res[f'prevalence_risk_group_{g}_{sex}'][ti] = sc.safedivide(np.count_nonzero(infected & grp), np.count_nonzero(grp))
                    res[f'new_infections_risk_group_{g}_{sex}'][ti] = w*np.count_nonzero(newly & grp)
        return

    def flat_results(self):
        return super().flat_results() + [
            ('prevalence', 'prevalence', True, 'Prevalence'),
        ]

    def finalize_results(self):
        """ As BaseSTI, with prevalence of everyone (not only adults), as in STIsim 3; and HIV deaths by cause """
        res = self.results
        n_sus = np.maximum(np.concatenate([res.n_S.values[:1], res.n_S.values[:-1]]), 1e-12) # Susceptible at the start of each step (no infections at the start)
        vals = res.new_infections.values/n_sus # UPSTREAM: a sum of flows over a state (new_infections/S) cannot be declared as a result
        self.results += ss.Result('incidence', values=vals, scale=False, label='Incidence', timevec=res.n_S.timevec, module=self.name)
        super().finalize_results()
        res.prevalence.values[:] = (res.n_infected.values/np.maximum(res.n_S.values + res.n_infected.values, 1e-12))
        res.new_deaths.values[:] = res.aids_death.values + res.cd4_death.values + sum(self.sim.results[mod][key].values for mod,key in self.other_deaths) # And deaths that other modules count as HIV deaths (e.g. sti.ART's, after stopping)
        res.cum_deaths.values[:] = np.cumsum(res.new_deaths.values)
        return

    def plot(self):
        """ Plot key HIV results """
        import matplotlib.pyplot as plt
        res = self.results
        keys = [k for k in ['new_infections', 'new_deaths', 'prevalence', 'prevalence_15_49', 'new_diagnoses'] if k in res]
        with sc.options.with_style('fancy'):
            fig, axs = sc.getrowscols(len(keys), make=True)
            for ax, k in zip(axs.flatten(), keys):
                ax.plot(res.timevec, res[k], lw=2, alpha=0.8)
                ax.set_title(res[k].label)
            sc.figlayout()
        return ss.return_fig(fig)

    # States as STIsim 3 named them, as agents, for connectors and custom code
    @property
    def susceptible(self): return self.S
    @property
    def infected(self):    return self.acute | self.latent | self.falling
    @property
    def infectious(self):  return self.infected
    @property
    def symptomatic(self): return self.infected
    @property
    def treatable(self):   return self.infected
    @property
    def ti_exposed(self):
        """ The time of acquisition (entering acute), as STIsim 3's BaseSTI named it """
        return self.ti_acute
    @property
    def ti_infected(self):
        """ The time of infection (entering acute), for code written for STIsim 3 """
        return self.ti_acute
    @property
    def aids(self):
        """ People with CD4 below 200 """
        out = self.S.asnew(np.zeros(len(self.S.raw), dtype=bool))
        uids = self.infected.uids
        out[uids] = self.cd4_values(uids) < 200
        return out

    def cd4_values(self, uids=None):
        """ The CD4 count of some people (default: everyone infected), with the modifiers on it (e.g. ART's) """
        if uids is None:
            uids = self.infected.uids
        return self.engine.attr_values(self.graph.attrs.cd4, ss.uids(uids))
