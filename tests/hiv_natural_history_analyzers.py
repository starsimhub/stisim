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
    Records an agent-uid-keyed dict of timeseries of CD4 count. Results obtainable by analyzer key 'hiv.ts_cd4'
    """

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.results['hiv.ts_cd4'] = {}

    def update_results(self):
        hiv = self.sim.diseases.hiv
        infected_uids = hiv.infected.uids
        for uid in infected_uids:
            if uid not in self.results['hiv.ts_cd4']:
                self.results['hiv.ts_cd4'][uid] = []
            cd4 = hiv.cd4[uid]
            self.results['hiv.ts_cd4'][uid].append(cd4)


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


class MTCTTracker(ss.Analyzer):
    """
    Records the number of mother-to-child hiv transmissions per timestep, accessible by analyzer key
    'hiv.n_mtc_transmissions'. MTCT definition used is "unborn child becomes infected".
    """

    result_name = 'hiv.n_mtc_transmissions'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.has_results = False

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result(self.result_name, dtype=list, scale=False))

    def update_results(self):
        hiv = self.sim.diseases.hiv
        people = self.sim.people
        infected = hiv.infected

        # we count infections this step for unborn children (transmitted by mother)
        if self.has_results:
            infected_this_step = infected & (hiv.ti_infected == self.ti)
        else:
            infected_this_step = (hiv.ti_infected <= self.ti)
            self.has_results = True
        mtct_this_step = (infected_this_step & (people.age < 0)).uids

        self.results[self.result_name][self.ti] = len(mtct_this_step)


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
