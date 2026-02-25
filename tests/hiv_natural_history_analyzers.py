import starsim as ss


class TimeToAIDSTracker(ss.Analyzer):
    # Records the time to AIDS for each infected model agent. Results are obtainable from the analyzer by key
    # 'hiv.ti_to_aids' .

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # This is a toggle that allows the first time step to count pre-simulation infections (once)
        self.has_results = False

    def step(self):
        pass

    def init_results(self):
        super().init_results()
        # results are a list of times to AIDS for agents infected at each timestep
        self.define_results(
            ss.Result('hiv.ti_to_aids', dtype=list, scale=False),
        )

    def update_results(self):
        infected = self.sim.people.filter('hiv.infected')  # self.sim.people.hiv.infected.uids
        if self.has_results:
            infected_this_step = infected('hiv.ti_infected') == self.ti
        else:
            infected_this_step = infected('hiv.ti_infected') <= self.ti
            self.has_results = True

        times_to_aids = infected_this_step.states['hiv.ti_falling'] - infected_this_step.states['hiv.ti_infected']
        self.results['hiv.ti_to_aids'][self.ti] = times_to_aids


class CD4ByUIDTracker(ss.Analyzer):
    # Records an agent-uid-keyed dict of timeseries of CD4 count. Results obtainable by analyzer key 'hiv.ts_cd4'

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
    # records the rel_trans (infectivity ratio) of agents in the specified states (acute, falling, and/or latent).
    # Results are obtainable by the analyzer keys below.
    STATES = {
        'acute':   {'result_name': 'hiv.acute_infectivity_ratios',   'filter': lambda hiv: hiv.acute},
        'falling': {'result_name': 'hiv.falling_infectivity_ratios', 'filter': lambda hiv: hiv.falling},
        'latent':  {'result_name': 'hiv.latent_infectivity_ratios',  'filter': lambda hiv: hiv.latent}
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
            # we ignore ti_infected < self.ti because relative transmission is updated BEFORE infection in the model
            # timestep and all just-infected agents have rel_trans == 1 (but this value is not used in any infection
            # events, as rel_trans will be properly updated for the NEXT transmission step)
            ratios = hiv.rel_trans[(state_dict['filter'](hiv=hiv) & (hiv.ti_infected < self.ti))]
            self.results[state_dict['result_name']][self.ti] = ratios


class TransmissionTracker(ss.Analyzer):
    # records the number of hiv transmissions per timestep, accessible by analyzer key 'hiv.n_transmissions' .
    def step(self):
        pass

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('hiv.n_transmissions', dtype=list, scale=False))

    def update_results(self):
        hiv = self.sim.diseases.hiv
        transmissions = len((hiv.ti_infected == self.ti).uids)
        self.results['hiv.n_transmissions'][self.ti] = transmissions
