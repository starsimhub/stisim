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


class SyphDx(STIDx):
    """
    A syphilis test, with the probability of a positive result by state

    As `sti.STIDx`, for `sti.Syphilis`: the dataframe's states can be STIsim 3's names
    (`primary`, `secondary`, `latent`, `early`, `late`, `tertiary`, `exposed`,
    `sus_not_naive`, ...) or the states of `sti.Syphilis`; the row for `susceptible` gives
    1 - specificity. A non-treponemal (RPR) test is `ss.Diagnostic(detects=syph.rpr, ...)`.

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
    `-fetus_age_cutoff_treat_eff` years of birth. Delivered over the maternal network to
    the fetuses of the people `who=` selects (e.g. `who=syph_tx.cured`); see `sti.SyphTx`.

    Results: those of `ss.treat`, and `new_fetus_treated_success`, `new_fetus_treated_failure`, and `new_fetus_treated_unnecessary` as in STIsim 3.
    """
    def __init__(self, fetus_treat_eff=0.98, treat_eff_reduced=0.75, fetus_age_cutoff_treat_eff=-0.25, disease='syph', network='maternalnet', name='fetal_tx', **kwargs):
        eff = lambda sim: np.where(sim.people.age.raw > fetus_age_cutoff_treat_eff, treat_eff_reduced, fetus_treat_eff) # Lower in the last trimester
        product = ss.Treatment(moves={f'{disease}.fetal_cure': eff}, name=f'{name}_drug')
        super().__init__(product, trace=network, disease=disease, name=name, **kwargs)
        return

    def finalize_results(self):
        super().finalize_results()
        add_aliases(self, dict(new_fetus_treated_success='cured', new_fetus_treated_failure='failed', new_fetus_treated_unnecessary='unneeded'))
        return


class SyphTx(STITreatment):
    """
    Syphilis treatment (benzathine penicillin), from any stage

    Cures with probability `treat_eff`: people in an early stage (incubation to early
    latent) move along `syph.cure_early` to Sp, and their RPR titre falls within months;
    people in late latent or tertiary move along `syph.cure_late` to Spl, where it stays.
    The fetuses of mothers it cures are treated by a step of their own, `fetal`
    (`sti.FetalSyphTx`, over the maternal network): pass it to the sim too. Who is treated
    is `who=` (e.g. `who=test.positive`) or `eligibility=`, as for `sti.STITreatment`.

    Args:
        pars (dict): `treat_eff`, `treat_prob`, `rel_treat_prob`, `fetus_treat_eff`, `treat_eff_reduced`, `fetus_age_cutoff_treat_eff`
        max_capacity (float): the most that can be treated (see `ss.Delivery` `capacity=`)
        years (list): the first and last years of treatment
        eligibility (Cond/func): who is treated, as a condition
        kwargs (dict): parameter values, or passed to `sti.STITreatment` (e.g. `who=`)

    Examples:
        ```python
        test = sti.SyphTest(product=dual, test_prob_data=0.1, name='syph_test')
        tx = sti.SyphTx(who=test.positive, name='syph_tx')
        sim = sti.Sim(diseases='syph', interventions=[test, tx, tx.fetal])
        ```
    """
    def __init__(self, pars=None, max_capacity=None, years=None, eligibility=None, name=None, **kwargs):
        p = sc.mergedicts(dict(rel_treat_prob=1, treat_prob=1, treat_eff=0.98, fetus_treat_eff=0.98, treat_eff_reduced=0.75, fetus_age_cutoff_treat_eff=-0.25), pars)
        for key in list(kwargs):
            if key in p or key == 'nontrep_revert_months':
                p[key] = kwargs.pop(key)
        if 'nontrep_revert_months' in p:
            errormsg = 'nontrep_revert_months is now dur_nontrep_revert of sti.Syphilis'
            raise ValueError(errormsg)
        if 'fetus' in kwargs:
            errormsg = 'The treatment of fetuses is a step of its own, tx.fetal: pass it to the sim, or leave it out'
            raise TypeError(errormsg)
        p = sc.objdict({k:p_of(v) for k,v in p.items()})
        self.syph_pars = p
        eff = p.treat_eff
        name = sc.ifelse(name, 'syphtx')
        product = ss.Treatment(moves={'syph.cure_early':eff, 'syph.cure_late':eff}, name=f'{name}_bpg')
        super().__init__(product=product, treat_prob=p.treat_prob*p.rel_treat_prob, eligibility=eligibility, max_capacity=max_capacity,
                         years=years, name=name, **kwargs)
        self.treats = ['syph']
        self.fetal = FetalSyphTx(fetus_treat_eff=p.fetus_treat_eff, treat_eff_reduced=p.treat_eff_reduced, fetus_age_cutoff_treat_eff=p.fetus_age_cutoff_treat_eff,
                                 who=self.cured, name=f'{name}_fetus')
        return


class NewbornTreatment(STITreatment):
    """
    Treatment of newborns with congenital syphilis (`syph.cs_cure`), with efficacy `treat_eff`

    Who is treated is `who=` (e.g. `who=newborn_test.positive`) or `eligibility=`.

    Args:
        pars (dict): `treat_eff`, `treat_prob`
        eligibility (Cond/func): who is treated, as a condition
        max_capacity (float): the most that can be treated
        years (list): the first and last years of treatment
        kwargs (dict): passed to `sti.STITreatment` (e.g. `who=`)
    """
    def __init__(self, pars=None, eligibility=None, max_capacity=None, years=None, name=None, **kwargs):
        p = sc.mergedicts(dict(treat_eff=0.98, treat_prob=1), pars, {k:kwargs.pop(k) for k in ['treat_eff', 'treat_prob'] if k in kwargs})
        kwargs.pop('by_sex', None)
        name = sc.ifelse(name, 'newborntreatment')
        product = ss.Treatment(moves={'syph.cs_cure': p_of(p['treat_eff'])}, name=f'{name}_drug')
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
        kwargs (dict): passed to `sti.STITest` (e.g. `years=`, `start=`, `eligibility=`)
    """
    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None,
                 name=None, label=None, newborn_test=None, dt_scale=True, **kwargs):
        pars = sc.mergedicts(pars)
        dt_scale = pars.pop('dt_scale', dt_scale)
        if newborn_test is not None:
            errormsg = 'newborn_test= is gone: sti.NewbornSyphTest(anc=anc_test) is a step of its own'
            raise ValueError(errormsg)
        if isinstance(test_prob_data, pd.DataFrame):
            self.by_group = test_prob_data
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
        """ The probability of testing for each agent, from a dataframe by risk group, sex, and sex work, as a typed probability (so it is offered at every step) """
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
        return ss.probperyear(out) if self.dt_scale else ss.prob(out, sim.t.dt)


class ANCSyphTest(SyphTest):
    """
    Syphilis screening at antenatal care, once a pregnancy

    Each pregnant woman is tested with probability `test_prob_data` at a visit `visit_ga`
    after conception (drawn for each pregnancy), if she is still pregnant then. The
    treatment of those who test positive is a step of its own (`who=anc.positive`), as is
    the testing of their newborns (`sti.NewbornSyphTest(anc=anc)`). Needs `ss.Pregnancy`;
    agents only.

    Args:
        test_prob_data (float/array): the probability of testing at the visit
        visit_ga (dur/Dist): the time from conception to the visit (default 1 to 8 months)
        kwargs (dict): passed to `sti.SyphTest`

    Examples:
        ```python
        anc = sti.ANCSyphTest(product=dual, test_prob_data=0.9, start=2005, name='anc')
        tx = sti.SyphTx(who=anc.positive, name='syph_tx')
        nb = sti.NewbornSyphTest(anc=anc, product=dual, test_prob_data=0.8, name='nb_test')
        nb_tx = sti.NewbornTreatment(who=nb.positive, name='newborn_tx')
        sim = sti.Sim(diseases='syph', interventions=[anc, tx, tx.fetal, nb, nb_tx])
        ```
    """
    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None,
                 name=None, label=None, visit_ga=None, **kwargs):
        for key in ['newborn_test', 'on_positive']:
            if kwargs.get(key) is not None:
                errormsg = f'{key}= is gone: what follows a positive test is a step of its own (who=anc.positive, or sti.NewbornSyphTest(anc=anc))'
                raise TypeError(errormsg)
        visit_ga = sc.ifelse(visit_ga, ss.uniform(ss.months(1), ss.months(8)))
        pregnant = ss.state('pregnancy.P')
        who = pregnant.enter & pregnant # At conception; tested at the visit if still pregnant
        super().__init__(test_prob_data=test_prob_data, years=years, start=start, stop=stop, pars=pars, product=product,
                         eligibility=eligibility, name=name, label=label, dt_scale=None, who=who, delay=visit_ga, **kwargs)
        return


class NewbornSyphTest(SyphTest):
    """
    Syphilis testing of newborns at birth, if their mothers tested positive at antenatal care within the last year

    Each such newborn is tested with probability `test_prob_data`; their treatment is a step
    of its own (e.g. `sti.NewbornTreatment(who=nb.positive)`). Agents only.

    Args:
        test_prob_data (float): the probability of testing
        anc (`sti.ANCSyphTest`): the antenatal testing whose positives' newborns are tested
        kwargs (dict): passed to `sti.SyphTest`
    """
    def __init__(self, test_prob_data=None, anc=None, **kwargs):
        if anc is None:
            errormsg = 'NewbornSyphTest tests the newborns of mothers who tested positive at an sti.ANCSyphTest: give it as anc='
            raise ValueError(errormsg)
        if kwargs.get('eligibility') is not None:
            errormsg = 'NewbornSyphTest tests the newborns of mothers who tested positive at ANC; it takes no eligibility='
            raise ValueError(errormsg)
        kwargs.setdefault('product', SyphDx(sensitivity=dict(CS=1.0, Ccs=1.0, Cnnd=1.0)))
        on = ss.state('pregnancy.birth') & ss.mother(anc.positive.since <= ss.years(1), born=True) # Born now, to a mother who tested positive at ANC in the last year
        super().__init__(test_prob_data=test_prob_data, dt_scale=None, on=on, **kwargs)
        return
