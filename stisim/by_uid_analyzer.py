import starsim as ss


class ByUIDAnalyzer(ss.Analyzer):
    """
    An analyzer that generates and reports custom time-index-series data by (custom) agent uid.

    Results are stored on a self.results dict per uid (as keys), each entry storing a ti (time index) series and
    corresponding value series.

    For example, if agent uid 5 & 10 are tracked and reports cd4 values for a three-timestep sim, self.results will
    store the following:
    {
        5:  {'ti': [0, 1, 2], 'val': [680.5, 672.2, 658.3]},
        10: {'ti': [0, 1, 2], 'val': [550.1, 532.2, 501.0]},
     }

    This result COULD be the result of specifying selector and subpop constructor arguments like so:

    # identifies the reportable data as agent CD4 values
    selector = lambda sim: sim.diseases.hiv.cd4

    # identifies HIV+ agents not on ART
    subpop = lambda sim: ( sim.diseases.hiv.infected & (~sim.diseases.hiv.on_art) ).uids

    """

    value_key = "val"
    ti_key = "ti"

    def __init__(self, selector, subpop=None, *args, **kwargs):
        """
        Args:
            selector: (callable) a function that returns desired data to be tracked by uid, given a sim object
                For example, the CD4 counts agents:
                    lambda sim: sim.diseases.hiv.cd4
            subpop: (optional) (callable, iterable) an iterable of uids OR a function that returns a an iterable of uids
                given a sim object.
                For example, HIV infected agents that are not on ART
                    lambda sim: ( sim.diseases.hiv.infected & (~sim.diseases.hiv.on_art) ).uids
                Default: all sim uids
        """
        super().__init__(*args, **kwargs)
        if subpop is None:
            subpop = lambda sim: sim.people.uids
        self.subpop = subpop
        self.selector = selector

    def step(self):
        pass

    def init_results(self):
        self.results = {}

    @classmethod
    def _init_result(cls):
        return {cls.value_key: [], cls.ti_key: []}

    @classmethod
    def _update_result(cls, result, value, ti):
        result[cls.value_key].append(value)
        result[cls.ti_key].append(ti)

    def update_results(self):
        uids = self.subpop(sim=self.sim)
        selected_data = self.selector(sim=self.sim)

        # record selected data by uid (by agent)
        for uid in uids:
            if uid not in self.results:
                self.results[uid] = self._init_result()
            self._update_result(result=self.results[uid], value=selected_data[uid], ti=self.ti)
