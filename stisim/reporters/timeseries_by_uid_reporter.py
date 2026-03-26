from stisim.reporters.itimeseries_reporter import ITimeSeriesReporter


class TimeSeriesByUIDReporter(ITimeSeriesReporter):
    """
    An analyzer that generates and reports custom time-index-series data by (custom) agent uid.

    Results are stored on a self.results dict per uid (as keys), each entry storing a ti (time index) series and
    corresponding value series.

    For example, if agent uid 5 & 10 are tracked and report cd4 values for a three-timestep sim, self.results['ts'] will
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
    def set_result(self, sim, data, uids):
        # record selected data by uid (by agent)
        for uid in uids:
            if uid not in self.ts_results:
                self.ts_results[uid] = self._init_result()
            self._set_result(result=self.ts_results[uid], value=data[uid], ti=self.ti)
