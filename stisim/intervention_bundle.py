import starsim as ss

class InterventionBundle(ss.Intervention):
    """
    This class simply abstracts a group of Intervention objects into a single "Intervention". This enables a form
    of intervention dependency, where contained interventions can reference each other. Because a user defines an
    InterventionBundle, they know for sure the order involved and can more easily setup conditional effects during
    execution of the interventions.
    """
    def __init__(self, interventions: list, states: list = None, pars=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.interventions = interventions
        for intv in self.interventions:
            intv.parent = self

        self.update_pars(pars, **kwargs)

        states = [] if states is None else states
        states = [ss.BoolState(name=name) for name in states]
        self.define_states(*states)

    def init_pre(self, sim):
        super().init_pre(sim=sim)
        for intervention in self.interventions:
            intervention.init_pre(sim=sim)

    def init_post(self):
        super().init_post()
        for intervention in self.interventions:
            intervention.init_post()

    def start_step(self):
        super().start_step()
        for intervention in self.interventions:
            intervention.start_step()

    def step(self):
        for intervention in self.interventions:
            intervention.step()

    def finish_step(self):
        super().finish_step()
        for intervention in self.interventions:
            intervention.finish_step()

    def update_results(self):
        super().update_results()
        for intervention in self.interventions:
            intervention.update_results()

    def finalize(self):
        super().finalize()
        for intervention in self.interventions:
            intervention.finalize()

    def check_method_calls(self):
        missing = super().check_method_calls()
        for intervention in self.interventions:
            missing += intervention.check_method_calls()
        return missing