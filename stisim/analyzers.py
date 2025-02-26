"""
Common analyzers for STI analyses
"""

# %% Imports and settings
import numpy as np
import sciris as sc
import starsim as ss
import pandas as pd

import stisim as sti
import pylab as pl

__all__ = ["result_grouper", "coinfection_stats", "sw_stats"]


class result_grouper(ss.Analyzer):
    @staticmethod
    def cond_prob(numerator, denominator):
        numer = len((denominator & numerator).uids)
        denom = len(denominator.uids)
        out = sc.safedivide(numer, denom)
        return out


# TODO: generalize this for any 2 diseases
class coinfection_stats(result_grouper):
    def __init__(self, diseases=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'coinfection_stats'
        self.diseases = diseases
        return

    def init_results(self):
        results = [
            ss.Result('syph_prev_no_hiv', dtype=float, scale=False),
            ss.Result('syph_prev_has_hiv', dtype=float, scale=False),
            ss.Result('syph_prev_no_hiv_f', dtype=float, scale=False),
            ss.Result('syph_prev_has_hiv_f', dtype=float, scale=False),
            ss.Result('syph_prev_no_hiv_m', dtype=float, scale=False),
            ss.Result('syph_prev_has_hiv_m', dtype=float, scale=False),
            ss.Result('syph_inf_fsw', dtype=float, scale=True),
            ss.Result('syph_inci_fsw', dtype=float, scale=False),
        ]
        self.define_results(*results)
        return

    def step(self):
        sim = self.sim
        ti = self.ti
        hiv = sim.diseases.hiv
        syph = sim.diseases.syphilis
        ppl = sim.people

        denom = (ppl.age >= 15) & (ppl.age < 50)  # Adults
        has_hiv = denom & hiv.infected  # Adults with HIV
        has_syph = denom & syph.infected  # Adults with syphilis

        has_syph_f = has_syph & ppl.female  # Women with syphilis
        has_syph_m = has_syph & ppl.male  # Men with syphilis
        has_hiv_f = has_hiv & ppl.female  # Women with HIV
        has_hiv_m = has_hiv & ppl.male  # Men with HIV
        no_hiv = denom & hiv.susceptible  # Adults without HIV
        no_hiv_f = no_hiv & ppl.female  # Women without HIV
        no_hiv_m = no_hiv & ppl.male  # Men without HIV

        self.results['syph_prev_no_hiv'][ti] = self.cond_prob(has_syph, no_hiv)
        self.results['syph_prev_has_hiv'][ti] = self.cond_prob(has_syph, has_hiv)
        self.results['syph_prev_no_hiv_f'][ti] = self.cond_prob(has_syph_f, no_hiv_f)
        self.results['syph_prev_has_hiv_f'][ti] = self.cond_prob(has_syph_f, has_hiv_f)
        self.results['syph_prev_no_hiv_m'][ti] = self.cond_prob(has_syph_m, no_hiv_m)
        self.results['syph_prev_has_hiv_m'][ti] = self.cond_prob(has_syph_m, has_hiv_m)

        return


class sw_stats(result_grouper):
    def __init__(self, diseases=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'sw_stats'
        self.diseases = diseases
        return

    def init_results(self):
        results = sc.autolist()
        for d in self.diseases:
            results += [
                ss.Result('share_new_infections_fsw_'+d, scale=False),
                ss.Result('share_new_infections_client_'+d,scale=False),
                ss.Result('new_infections_fsw_'+d, dtype=int),
                ss.Result('new_infections_client_'+d, dtype=int),
                ss.Result('new_infections_non_fsw_'+d, dtype=int),
                ss.Result('new_infections_non_client_'+d, dtype=int),
                ss.Result('new_transmissions_fsw_'+d, dtype=int),
                ss.Result('new_transmissions_client_'+d, dtype=int),
                ss.Result('new_transmissions_non_fsw_'+d, dtype=int),
                ss.Result('new_transmissions_non_client_'+d, dtype=int),
            ]
        self.define_results(*results)
        return

    def step(self):
        sim = self.sim
        ti = self.ti

        if ti > 0:

            for d in self.diseases:
                dis = sim.diseases[d]
                nw = sim.networks.structuredsexual
                fsw = nw.fsw
                client = nw.client
                non_fsw = sim.people.female & ~nw.fsw
                non_client = sim.people.male & ~nw.client
                newly_infected = dis.ti_infected == ti
                total_acq = len(newly_infected.uids)

                newly_transmitting_fsw = (dis.ti_transmitted == ti) & fsw
                newly_transmitting_clients = (dis.ti_transmitted == ti) & client
                newly_transmitting_non_fsw = (dis.ti_transmitted == ti) & non_fsw
                newly_transmitting_non_client = (dis.ti_transmitted == ti) & non_client

                new_transmissions_fsw = dis.new_transmissions[newly_transmitting_fsw]
                new_transmissions_client = dis.new_transmissions[newly_transmitting_clients]
                new_transmissions_non_fsw = dis.new_transmissions[newly_transmitting_non_fsw]
                new_transmissions_non_client = dis.new_transmissions[newly_transmitting_non_client]

                self.results['share_new_infections_fsw_'+d][ti] = self.cond_prob(fsw, newly_infected)
                self.results['share_new_infections_client_'+d][ti] = self.cond_prob(client, newly_infected)

                self.results['new_infections_fsw_'+d][ti] = len((fsw & newly_infected).uids)
                self.results['new_infections_client_'+d][ti] = len((client & newly_infected).uids)
                self.results['new_infections_non_fsw_'+d][ti] = len((non_fsw & newly_infected).uids)
                self.results['new_infections_non_client_'+d][ti] = len((non_client & newly_infected).uids)

                self.results['new_transmissions_fsw_'+d][ti] = sum(new_transmissions_fsw)
                self.results['new_transmissions_client_'+d][ti] = sum(new_transmissions_client)
                self.results['new_transmissions_non_fsw_'+d][ti] = sum(new_transmissions_non_fsw)
                self.results['new_transmissions_non_client_'+d][ti] = sum(new_transmissions_non_client)

                total_trans = sum(new_transmissions_fsw) + sum(new_transmissions_client) + sum(new_transmissions_non_fsw) + sum(new_transmissions_non_client)
                if total_trans != len(newly_infected.uids):
                    errormsg = f'Infections acquired should equal number transmitted: {total_acq} vs {total_trans}'
                    raise ValueError(errormsg)

        return



