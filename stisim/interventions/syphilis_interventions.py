"""
Syphilis testing and treatment, as Starsim products and deliveries
"""

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from stisim.interventions.sti_interventions import STIDx, STITest, STITreatment, add_aliases

__all__ = ['SyphDx', 'SyphTx', 'FetalSyphTx', 'NewbornTreatment', 'SyphTest', 'ANCSyphTest', 'NewbornSyphTest']

# The states of sti.Syphilis that STIsim 3's state names cover
syph_states = dict(
    naive         = ['S'],
    sus_not_naive = ['Sp', 'Spl', 'Sr'],
    exposed       = ['E'],
    primary       = ['Pv', 'Ph'],
    secondary     = ['S2v', 'S2h'],
    latent        = ['EL', 'LL'],
    early         = ['EL'],
    late          = ['LL'],
    tertiary      = ['T'],
    symptomatic   = ['Pv', 'Ph', 'S2v', 'S2h'],
    ulcerative    = ['Pv'],
    visibly_symptomatic = ['Pv', 'S2v'],
    mat_active    = ['E', 'Pv', 'Ph', 'S2v', 'S2h'],
    sexually_transmissible = ['Pv', 'Ph', 'S2v', 'S2h', 'EL'],
    infectious    = ['Pv', 'Ph', 'S2v', 'S2h', 'EL', 'LL'],
    infected      = ['E', 'Pv', 'Ph', 'S2v', 'S2h', 'EL', 'LL', 'T'],
    congenital    = ['CS'],
)


def p_of(val):
    """ The probability of a Starsim 3 bernoulli, or the value """
    return val.pars.p if isinstance(val, ss.bernoulli) else val


def add_module(sim, mod):
    """ Add a delivery that another one made to the sim's interventions. TODO(v4-rs): a module that adds modules it needs """
    if mod.name not in sim.interventions:
        sim.interventions[mod.name] = mod
        mod.init_pre(sim)
    return


class SyphDx(STIDx):
    """
    A syphilis test, with the probability of a positive result by state

    As `sti.STIDx`, for `sti.Syphilis`: the dataframe's states can be STIsim 3's names
    (`primary`, `secondary`, `latent`, `early`, `late`, `tertiary`, `exposed`,
    `sus_not_naive`, ...) or the states of `sti.Syphilis`; the row for `susceptible` gives
    1 - specificity. A non-treponemal (RPR) test is `ss.Test(detects=syph.rpr, ...)`.

    Args:
        df (dataframe): columns `state` and `p_positive`
        kwargs (dict): passed to `sti.STIDx` (e.g. `sensitivity=`, `name=`)

    Examples:
        ```python
        dual = sti.SyphDx(sensitivity=dict(Sp=0.95, Spl=0.95, E=0.2, Pv=0.2, Ph=0.2, S2v=0.95, S2h=0.95, EL=0.95, LL=0.95, T=0.95), specificity=0.99)
        ```
    """
    state_map = syph_states

    def __init__(self, df=None, *args, disease='syph', **kwargs):
        super().__init__(df, disease, *args, **kwargs)
        return


class FetalSyphTx(ss.treat):
    """
    Treatment of the fetuses of mothers cured of syphilis: a fetus infected in utero is
    cured (`syph.fetal_cure`) with `fetus_treat_eff`, or `treat_eff_reduced` within
    `-fetus_age_cutoff_treat_eff` years of birth; made by `sti.SyphTx`

    Results: those of `ss.treat`, and `new_fetus_treated_success`, `new_fetus_treated_failure`, and `new_fetus_treated_unnecessary` as in STIsim 3.
    """
    def __init__(self, fetus_treat_eff=0.98, treat_eff_reduced=0.75, fetus_age_cutoff_treat_eff=-0.25, disease='syph', network='maternalnet', **kwargs):
        self.effs = (fetus_treat_eff, treat_eff_reduced, fetus_age_cutoff_treat_eff)
        product = ss.Treatment(moves=f'{disease}.fetal_cure', efficacy=self.fetal_eff, name=kwargs.get('name', 'fetal_tx') + '_drug')
        on = ss.state(f'{disease}.cure_early') | ss.state(f'{disease}.cure_late')
        super().__init__(product, on=on, trace=network, disease=disease, **kwargs)
        return

    def fetal_eff(self, sim):
        """ Lower in the last trimester """
        full, reduced, cutoff = self.effs
        return np.where(sim.people.age.raw > cutoff, reduced, full)

    def finalize_results(self):
        super().finalize_results()
        add_aliases(self, dict(new_fetus_treated_success='cured', new_fetus_treated_failure='failed', new_fetus_treated_unnecessary='overtreated'))
        return


class SyphTx(STITreatment):
    """
    Syphilis treatment (benzathine penicillin), from any stage

    Cures with probability `treat_eff`: people in an early stage (incubation to early
    latent) move along `syph.cure_early` to Sp, and their RPR titre falls within months;
    people in late latent or tertiary move along `syph.cure_late` to Spl, where it stays.
    Fetuses of mothers cured while pregnant are treated too (`sti.FetalSyphTx`, added to
    the sim when it has a maternal network). As for `sti.STITreatment`, people are treated
    when referred (e.g. by a test's `on_positive=`), or at each step they meet `eligibility`.

    Args:
        pars (dict): `treat_eff`, `treat_prob`, `rel_treat_prob`, `fetus_treat_eff`, `treat_eff_reduced`, `fetus_age_cutoff_treat_eff`
        max_capacity (float): the most that can be treated (see `ss.Delivery` `capacity=`)
        years (list): the first and last years of treatment
        eligibility (Cond/func): who is treated at each step; if None, only those referred
        fetus (bool): whether to treat the fetuses of cured mothers
        kwargs (dict): parameter values, or passed to `sti.STITreatment`

    Examples:
        ```python
        tx = sti.SyphTx(name='syph_tx')
        test = sti.SyphTest(product=dual, test_prob_data=0.1, on_positive=tx)
        ```
    """
    def __init__(self, pars=None, max_capacity=None, years=None, eligibility=None, name=None, fetus=True, **kwargs):
        p = sc.mergedicts(dict(rel_treat_prob=1, treat_prob=1, treat_eff=0.98, fetus_treat_eff=0.98, treat_eff_reduced=0.75, fetus_age_cutoff_treat_eff=-0.25), pars)
        for key in list(kwargs):
            if key in p or key == 'nontrep_revert_months':
                p[key] = kwargs.pop(key)
        if 'nontrep_revert_months' in p:
            errormsg = 'nontrep_revert_months is now dur_nontrep_revert of sti.Syphilis'
            raise ValueError(errormsg)
        p = sc.objdict({k:p_of(v) for k,v in p.items()})
        self.syph_pars = p
        eff = p.treat_eff
        name = sc.ifelse(name, 'syphtx')
        product = ss.Treatment(moves={'syph.cure_early':eff, 'syph.cure_late':eff}, name=f'{name}_bpg')
        super().__init__(product=product, treat_prob=p.treat_prob*p.rel_treat_prob, eligibility=eligibility, max_capacity=max_capacity,
                         years=years, name=name, **kwargs)
        self.treats = ['syph']
        self.fetus = fetus
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        nets = [name for name,net in sim.networks.items() if isinstance(net, ss.MaternalNet)]
        if self.fetus and nets and 'syph' in sim.diseases:
            p = self.syph_pars
            fetal = FetalSyphTx(fetus_treat_eff=p.fetus_treat_eff, treat_eff_reduced=p.treat_eff_reduced, fetus_age_cutoff_treat_eff=p.fetus_age_cutoff_treat_eff,
                                network=nets[0], name=f'{self.name}_fetus')
            add_module(sim, fetal)
        return


class NewbornTreatment(STITreatment):
    """
    Treatment of newborns with congenital syphilis (`syph.cs_cure`), with efficacy `treat_eff`

    Treated when referred (e.g. by `sti.NewbornSyphTest`'s `on_positive=`), or at each step they meet `eligibility`.

    Args:
        pars (dict): `treat_eff`, `treat_prob`
        eligibility (Cond/func): who is treated at each step; if None, only those referred
        max_capacity (float): the most that can be treated
        years (list): the first and last years of treatment
        kwargs (dict): passed to `sti.STITreatment`
    """
    def __init__(self, pars=None, eligibility=None, max_capacity=None, years=None, name=None, **kwargs):
        p = sc.mergedicts(dict(treat_eff=0.98, treat_prob=1), pars, {k:kwargs.pop(k) for k in ['treat_eff', 'treat_prob'] if k in kwargs})
        kwargs.pop('by_sex', None)
        name = sc.ifelse(name, 'newborntreatment')
        product = ss.Treatment(moves='syph.cs_cure', efficacy=p_of(p['treat_eff']), name=f'{name}_drug')
        super().__init__(product=product, treat_prob=p_of(p['treat_prob']), eligibility=eligibility, max_capacity=max_capacity, years=years, name=name, **kwargs)
        self.treats = ['syph']
        return


class SyphTest(STITest):
    """
    Syphilis testing, with the probability of testing as for `sti.STITest`, or by risk group

    `test_prob_data` can also be a dataframe with columns `year`, `risk_group`, `sex`,
    `sw`, and `symp_test_prob`: an annual probability for each risk group and sex of
    `structuredsexual`, with the rows with `sw=1` for female sex workers and clients
    (interpolated between years). Agents only.

    Args:
        test_prob_data (float/array/dataframe): the probability of testing
        product (`sti.SyphDx`): the test
        kwargs (dict): passed to `sti.STITest` (e.g. `years=`, `start=`, `eligibility=`, `on_positive=`)
    """
    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None,
                 name=None, label=None, newborn_test=None, dt_scale=True, **kwargs):
        pars = sc.mergedicts(pars)
        dt_scale = pars.pop('dt_scale', dt_scale)
        if newborn_test is not None:
            errormsg = 'newborn_test= is for sti.ANCSyphTest'
            raise ValueError(errormsg)
        if isinstance(test_prob_data, pd.DataFrame):
            self.by_group = test_prob_data
            kwargs['repeat'] = True # A probability per agent each step
            test_prob_data = 1.0
        else:
            self.by_group = None
        name = sc.ifelse(name, type(self).__name__.lower()) # As in Starsim 3
        super().__init__(product=product if product is not None else SyphDx(sensitivity=dict.fromkeys(syph_states['infected'], 1.0)),
                         test_prob_data=test_prob_data, years=years, start=start, stop=stop, eligibility=eligibility, dt_scale=dt_scale,
                         disease='syph', name=name, label=label, **kwargs, **pars)
        if self.by_group is not None:
            self.coverage = self.group_prob
        return

    def group_prob(self, sim):
        """ The probability of testing in this step for each agent, from a dataframe by risk group, sex, and sex work """
        df = self.by_group
        nw = sim.networks.structuredsexual
        ppl = sim.people
        year = sim.t.now('year')
        out = np.zeros(len(ppl.uid.raw))
        for (rg, sex, sw), rows in df.sort_values('sw').groupby(['risk_group', 'sex', 'sw'], sort=False):
            rows = rows.sort_values('year')
            p = np.interp(year, rows.year.values, rows.symp_test_prob.values)
            cond = ppl[sex].raw & (nw.risk_group.raw == rg)
            if sw:
                cond &= (nw.fsw if sex == 'female' else nw.client).raw
            out[cond] = p
        out = np.clip(out*self.rel_test, 0, 1)
        return 1 - (1 - out)**sim.t.dt.years if self.dt_scale else out


class ANCSyphTest(SyphTest):
    """
    Syphilis screening at antenatal care, once a pregnancy

    Each pregnancy has a visit `visit_ga` steps after conception (drawn for each
    pregnancy), and a pregnant woman is tested there with probability
    `test_prob_data` (per pregnancy). Those who test positive go to `on_positive=` (e.g. an
    `sti.SyphTx`), and their newborns can be tested at birth (`newborn_test=`, an
    `sti.NewbornSyphTest`). Needs `ss.Pregnancy`.

    Args:
        test_prob_data (float/array): the probability of testing at the visit
        visit_ga (Dist): the number of steps from conception to the visit (default 1 to 7)
        newborn_test (`sti.NewbornSyphTest`): testing of the newborns of mothers who test positive
        kwargs (dict): passed to `sti.SyphTest`

    Examples:
        ```python
        tx = sti.SyphTx(name='syph_tx')
        nb_tx = sti.NewbornTreatment(name='newborn_tx')
        nb = sti.NewbornSyphTest(product=dual, test_prob_data=0.8, on_positive=nb_tx)
        anc = sti.ANCSyphTest(product=dual, test_prob_data=0.9, start=2005, on_positive=tx, newborn_test=nb)
        sim = sti.Sim(diseases='syph', interventions=[anc, nb])
        ```
    """
    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None,
                 name=None, label=None, newborn_test=None, visit_ga=None, on_positive=None, **kwargs):
        self.visit_ga = sc.ifelse(visit_ga, ss.randint(1, 8)) # Steps after conception, as in STIsim 3
        self.newborn_test = newborn_test
        on_positive = sc.tolist(on_positive) + [self.record_positive]
        super().__init__(test_prob_data=test_prob_data, years=years, start=start, stop=stop, pars=pars, product=product,
                         eligibility=None, name=name, label=label, dt_scale=False, on_positive=on_positive, **kwargs)
        self.user_eligible = eligibility
        self.eligible = ss.conditions.as_cond(self.visiting) # At the visit, once a pregnancy
        self.define_states(ss.FloatArr('ti_positive', label='Time of the last positive test'), ss.FloatArr('ti_visit', label='Time of the ANC visit'))
        if newborn_test is not None:
            newborn_test.anc = self
        return

    def init_pre(self, sim):
        preg = sim.demographics.get('pregnancy')
        if preg is None:
            errormsg = 'ANCSyphTest needs ss.Pregnancy'
            raise ValueError(errormsg)
        super().init_pre(sim)
        return

    def visiting(self, sim):
        """ Pregnant women at their visit. TODO(v4-rs): preg.attr() added by a delivery, so eligible=preg.P & (preg.P.time >= ga) """
        preg = sim.demographics.pregnancy
        new = (preg.ti_pregnant == sim.ti).uids
        if len(new):
            self.ti_visit[new] = sim.ti + self.visit_ga.rvs(new)
        out = np.asarray(preg.pregnant.raw, dtype=bool) & (self.ti_visit.raw == sim.ti)
        if self.user_eligible is not None:
            out &= np.asarray(self.user_eligible(sim), dtype=bool)
        return out

    def record_positive(self, sim, uids):
        self.ti_positive[uids] = sim.ti
        return


class NewbornSyphTest(SyphTest):
    """
    Syphilis testing of newborns at birth, if their mothers tested positive at antenatal care in this pregnancy

    Each such newborn is tested with probability `test_prob_data`; those who test positive
    go to `on_positive=` (e.g. an `sti.NewbornTreatment`). Pass it to an `sti.ANCSyphTest`
    as `newborn_test=` (or give `anc=`), and to the sim. Agents only.

    Args:
        test_prob_data (float): the probability of testing
        anc (`sti.ANCSyphTest`): the antenatal testing whose positives' newborns are tested
        kwargs (dict): passed to `sti.SyphTest`
    """
    def __init__(self, test_prob_data=None, anc=None, **kwargs):
        self.anc = anc
        kwargs.setdefault('product', SyphDx(sensitivity=dict(CS=1.0, Ccs=1.0, Cnnd=1.0)))
        if kwargs.get('eligibility') is not None:
            errormsg = 'NewbornSyphTest tests the newborns of mothers who tested positive at ANC; it takes no eligibility='
            raise ValueError(errormsg)
        super().__init__(test_prob_data=test_prob_data, dt_scale=False, on=ss.state('pregnancy.birth'), **kwargs)
        self.eligible = ss.conditions.as_cond(self.mother_positive)
        return

    def init_pre(self, sim):
        if self.anc is None:
            ancs = [i for i in sim.interventions.values() if isinstance(i, ANCSyphTest)]
            if not ancs:
                errormsg = 'NewbornSyphTest needs an sti.ANCSyphTest in the sim'
                raise ValueError(errormsg)
            self.anc = ancs[0]
        super().init_pre(sim)
        return

    def mother_positive(self, sim):
        """ Whose mother tested positive at antenatal care within the last year """
        ppl = sim.people
        anc = sim.interventions[self.anc.name]
        out = np.zeros(len(ppl.uid.raw), dtype=bool)
        kids = ppl.auids
        moms = ppl.parent.raw[kids]
        has = moms >= 0
        kids, moms = kids[has], ss.uids(moms[has])
        recent = (sim.ti - anc.ti_positive.raw[moms]) <= ss.years(1)/sim.t.dt
        out[kids[recent]] = True
        return out
