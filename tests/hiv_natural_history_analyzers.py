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
    ONART = 'onart'
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
    'hiv.n_mtc_transmissions' .
    """

    result_name = 'hiv.n_mtc_transmissions'

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result(self.result_name, dtype=list, scale=False))

    def update_results(self):
        hiv = self.sim.diseases.hiv
        transmitting_mothers = hiv.ti_transmitted_mtc == self.ti

        # now back-out mtc transmissions by removing any potential sexual transmissions
        # TODO: Update how mtc transmissions are counted once this issue is fixed:
        #  https://github.com/starsimhub/stisim/issues/325
        transmissions =  sum(hiv.new_transmissions[transmitting_mothers] - hiv.new_transmissions_sex[transmitting_mothers])
        self.results[self.result_name][self.ti] = transmissions


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