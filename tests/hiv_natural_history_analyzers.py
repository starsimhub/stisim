import starsim as ss


class TimeToAIDSTracker(ss.Analyzer):
    # Records the time to AIDS for each infected model agent. Results are obtainable from the analyzer by key
    # 'hiv.ti_to_aids' .

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # This is a toggle that allows the first time step to count pre-simulation infections (once)
        self.has_results = False


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
