from typing import Callable, Iterable, Dict

import sciris as sc
import starsim as ss


class TimeToAIDSTracker(ss.Analyzer):
    """
    Records the time to AIDS (falling) for each infected model agent. Results are obtainable from the analyzer by key 'hiv.ti_to_aids' .
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # This is a toggle that allows the first time step to count pre-simulation infections (once)
        self.has_results = False

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        # results are a list of times to AIDS for agents infected at each timestep
        self.define_results(ss.Result('hiv.ti_to_aids', dtype=list, scale=False))

    def update_results(self):
        hiv = self.sim.diseases.hiv
        infected = hiv.infected

        if self.has_results:
            infected_this_step = infected & (hiv.ti_infected == self.ti)
        else:
            infected_this_step = infected & (hiv.ti_infected <= self.ti)
            self.has_results = True

        times_to_aids = hiv.ti_falling[infected_this_step] - hiv.ti_infected[infected_this_step]
        self.results['hiv.ti_to_aids'][self.ti] = times_to_aids


class CD4ByUIDTracker(ss.Analyzer):
    """
    Records an agent-uid-keyed dict of timeseries of CD4 count for specified subpopulation. Results obtainable by
    analyzer key 'hiv.ts_cd4' .
    """
    result_name = 'hiv.ts_cd4'

    INFECTED = 'infected'
    ONART = 'on_art'
    SUBPOPS = [INFECTED, ONART]

    def __init__(self, subpop: str = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.subpop = self.INFECTED if subpop is None else subpop
        if self.subpop not in self.SUBPOPS:
            raise Exception(f"Unknown sub-population specified for analyzer CD4ByUIDTracker: {self.subpop}")

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.results[self.result_name] = {}

    def update_results(self):
        hiv = self.sim.diseases.hiv
        result_name = self.result_name

        # determine which agent uids are in the sub-population of interest
        if self.subpop == self.INFECTED:
            uids = hiv.infected.uids
        elif self.subpop == self.ONART:
            uids = hiv.on_art.uids
        else:
            raise Exception(f"Unknown sub-population specified for analyzer CD4ByUIDTracker: {self.subpop}")

        # record cd4 count for agents of interest
        for uid in uids:
            if uid not in self.results[result_name]:
                self.results[result_name][uid] = []
            cd4 = hiv.cd4[uid]
            self.results[result_name][uid].append(cd4)


class RelativeInfectivityByUIDTracker(ss.Analyzer):
    """
    Records an agent-uid-keyed dict of timeseries of CD4 count for specified subpopulation. Results obtainable by
    analyzer key 'hiv.ts_cd4' .
    """
    result_name = 'hiv.ts_rel_trans'

    INFECTED = 'infected'
    ONART = 'on_art'
    SUBPOPS = [INFECTED, ONART]

    def __init__(self, subpop: str = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.subpop = self.INFECTED if subpop is None else subpop
        if self.subpop not in self.SUBPOPS:
            raise Exception(f"Unknown sub-population specified for analyzer RelativeInfectivityByUIDTracker: {self.subpop}")

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.results[self.result_name] = {}

    def update_results(self):
        hiv = self.sim.diseases.hiv
        result_name = self.result_name

        # determine which agent uids are in the sub-population of interest
        if self.subpop == self.INFECTED:
            uids = hiv.infected.uids
        elif self.subpop == self.ONART:
            uids = hiv.on_art.uids
        else:
            raise Exception(f"Unknown sub-population specified for analyzer RelativeInfectivityByUIDTracker: {self.subpop}")

        # record rel_trans for agents of interest
        for uid in uids:
            if uid not in self.results[result_name]:
                self.results[result_name][uid] = []
            rel_trans = hiv.rel_trans[uid]
            self.results[result_name][uid].append(rel_trans)


# TODO: consider replacing usage of this analyzer with RelativeInfectivityByUIDTracker usage & updates
class RelativeInfectivityTracker(ss.Analyzer):
    """
    Records the rel_trans (infectivity ratio) of agents in the specified states (acute, falling, and/or latent).
    Results are obtainable by the analyzer keys below.
    """
    STATES = {
        'acute':   {'result_name': 'hiv.acute_rel_trans',   'filter': lambda hiv: hiv.acute},
        'falling': {'result_name': 'hiv.falling_rel_trans', 'filter': lambda hiv: hiv.falling},
        'latent':  {'result_name': 'hiv.latent_rel_trans',  'filter': lambda hiv: hiv.latent},
        'aids':    {'result_name': 'hiv.aids_rel_trans',    'filter': lambda hiv: hiv.aids}

    }

    def __init__(self, states: list[str], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.states_to_track = states

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        for state in self.states_to_track:
            if state not in self.STATES:
                raise Exception(f"Unknown infectivity state: {state}")
            state_dict = self.STATES[state]
            self.define_results(ss.Result(state_dict['result_name'], dtype=list, scale=False))

    def update_results(self):
        hiv = self.sim.diseases.hiv
        for state in self.states_to_track:
            state_dict = self.STATES[state]
            # we ignore ti_infected == self.ti because relative transmission is updated BEFORE infection in the model
            # timestep and all just-infected agents have rel_trans == 1 (but this value is not used in any infection
            # events, as rel_trans will be properly updated for the NEXT transmission step)
            ratios = hiv.rel_trans[(state_dict['filter'](hiv=hiv) & (hiv.ti_infected < self.ti))]
            self.results[state_dict['result_name']][self.ti] = ratios


class SexualTransmissionCountTracker(ss.Analyzer):
    """
    Records the number of sexual HIV transmissions per timestep, results obtainable by analyzer key
    'hiv.n_sexual_transmissions' .
    """

    result_name = 'hiv.n_sexual_transmissions'

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result(self.result_name, dtype=list, scale=False))

    def update_results(self):
        hiv = self.sim.diseases.hiv
        transmissions = sum(hiv.new_transmissions_sex.notnanvals)
        self.results[self.result_name][self.ti] = transmissions


class MTCTransmissionCountTracker(ss.Analyzer):
    """
    Records the number of mother-to-child HIV transmissions per timestep, results obtainable by analyzer key
    'hiv.n_mtc_transmissions'. Uses the built-in new_infections_mtct result from BaseSTI.
    """

    result_name = 'hiv.n_mtc_transmissions'

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result(self.result_name, dtype=list, scale=False))

    def update_results(self):
        self.results[self.result_name][self.ti] = self.sim.diseases.hiv.results['new_infections_mtct'][self.ti]


class PrepCoverageAnalyzer(ss.Analyzer):
    """
    Records the PrEP coverage per 1+ eligibiliity groups/functions per timestep, results obtainable by analyzer key
    'hiv.prep_coverage' . consider_new_infections=True yields results equivalent to "end of intervention time" instead
    of "end of model timestep time".
    """

    result_name = 'hiv.prep_coverage'

    def __init__(self, eligibilities: Dict[str, Callable], consider_new_infections=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.eligibilities = eligibilities
        self.consider_new_infections = consider_new_infections

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.results[self.result_name] = {group_name: [] for group_name, _ in self.eligibilities.items()}

    def update_results(self):
        sim = self.sim
        hiv = sim.diseases.hiv

        for group_name, eligibility in self.eligibilities.items():
            # The reason for "OR on prep + infected this ti" is because the model infection step occurs between:
            # - when the intervention uses the eligibility functions used here
            # - when analyzer update_results() is called
            # Including this additional OR clause allows for checking PrEP coverage at the time of the intervention (as
            # opposed to AFTER intervention + infection) which is useful for test code. Ignoring the additional clause
            # yields "end of timestep" results.
            if self.consider_new_infections:
                eligible_uids =  (eligibility(sim) | self.infected_this_ti_on_prep(sim)).uids
            else:
                eligible_uids = eligibility(sim).uids
            n_eligible = len(eligible_uids)
            n_on_prep = len(hiv.on_prep.uids & eligible_uids)
            coverage = n_on_prep / n_eligible
            self.results[self.result_name][group_name].append(coverage)

    @staticmethod
    def infected_this_ti_on_prep(sim):
        hiv = sim.diseases.hiv
        agents = hiv.on_prep & (hiv.ti_infected == sim.ti)
        return agents

class PrepEfficacyAnalyzer(ss.Analyzer):
    """
    Records the PrEP efficacy of agents on PrEP, per agent uid.
    """

    result_name = 'hiv.prep_efficacy'
    ti_name = 'ti'

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.results[self.result_name] = {}
        self.results[self.ti_name] = {}

    def update_results(self):
        sim = self.sim
        ti = sim.ti
        hiv = sim.diseases.hiv
        result_name = self.result_name
        ti_name = self.ti_name

        uids = hiv.on_prep.uids

        # record prep efficacy for agents of interest
        for uid in uids:
            if uid not in self.results[result_name]:
                self.results[result_name][uid] = []
                self.results[ti_name][uid] = []
            prep_eff = hiv.prep_eff[uid]
            self.results[result_name][uid].append(prep_eff)
            self.results[ti_name][uid].append(ti)


class PrepDurationAnalyzer(ss.Analyzer):
    """
    Records the time on PrEP of agents, per agent uid.
    """

    result_name = 'hiv.prep_duration'
    ti_name = 'ti'

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.results[self.result_name] = {}
        self.results[self.ti_name] = {}

    def update_results(self):
        sim = self.sim
        ti = sim.ti
        hiv = sim.diseases.hiv
        result_name = self.result_name
        ti_name = self.ti_name

        uids = hiv.on_prep.uids

        # record prep efficacy for agents of interest
        for uid in uids:
            if uid not in self.results[result_name]:
                self.results[result_name][uid] = []
                self.results[ti_name][uid] = []
            ti_on_prep = ti - hiv.ti_prep_start[uid]
            self.results[result_name][uid].append(ti_on_prep)
            self.results[ti_name][uid].append(ti)

# perinatal infection progression not currently implemented in hivsim, so leaving this untested analyzer out for
# future work
# class BirthTracker(ss.Analyzer):
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.uids_seen = {}
#
#     def step(self):
#         pass
#
#     def init_results(self):
#         super().init_results()
#         # self.define_results(ss.Result('hiv.perinatally_infected_uids', dtype=list, scale=False))
#         self.results['hiv.perinatally_infected_uids'] = []
#         self.results['hiv.tis_acute'] = {}
#         self.results['hiv.tis_latent'] = {}
#         self.results['hiv.tis_falling'] = {}
#
#     def update_results(self):
#         hiv = self.sim.diseases.hiv
#
#         tis_acute = hiv.ti_acute - hiv.ti_infected
#         tis_latent = hiv.ti_latent - hiv.ti_acute
#         tis_falling = hiv.ti_falling - hiv.ti_latent
#
#         # checking for any new infections and their computed times to different stages
#         infected_uids = hiv.infected.uids
#         for uid in infected_uids:
#             if uid not in self.uids_seen:
#                 self.results['hiv.tis_acute'][uid] = tis_acute[uid]
#                 self.results['hiv.tis_latent'][uid] = tis_latent[uid]
#                 self.results['hiv.tis_falling'][uid] = tis_falling[uid]
#                 self.uids_seen[uid] = None  # record to skip next time we see this uid
#
#         # checking birth stats for this timestep, looking for perinatally infected agents
#         people = self.sim.people
#         born = (people.age >= 0) & (people.age < float(self.sim.dt))  # agents where this is their first non-negative age
#         born_infected = born & (hiv.ti_infected <= self.ti)
#         born_infected_uids = born_infected.uids
#         self.results['hiv.perinatally_infected_uids'].extend(born_infected_uids)