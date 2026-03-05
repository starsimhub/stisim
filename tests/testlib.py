from typing import Union

import starsim as ss
import stisim as sti

DEFAULT_STR = 'default'


def build_testing_sim(n_agents: int = 500,
                      duration: int = 40,  # years
                      diseases:         Union[str, list, None]                 = DEFAULT_STR,
                      maternal_network: Union[str, ss.MaternalNet, None]       = DEFAULT_STR,
                      prior_network:    Union[str, sti.PriorPartners, None]    = DEFAULT_STR,
                      sexual_network:   Union[str, sti.StructuredSexual, None] = DEFAULT_STR,
                      pregnancy:        Union[str, ss.Pregnancy, None]         = DEFAULT_STR,
                      death:            Union[str, ss.Deaths, None]            = DEFAULT_STR,
                      interventions: list = None,
                      analyzers: list = None):
    """DEFAULT_STR values mean use default, None means off, otherwise use what is passed in"""

    # # setup disease(s) like HIV
    if diseases == DEFAULT_STR:
        diseases = [sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=0.05)]

    # setup networks
    networks = []
    if sexual_network == DEFAULT_STR:
        sexual_network = sti.StructuredSexual(recall_prior=True)
    if sexual_network is not None:
        networks.append(sexual_network)

    if prior_network == DEFAULT_STR:
        prior_network = sti.PriorPartners()
    if prior_network is not None:
        networks.append(prior_network)

    if maternal_network == DEFAULT_STR:
        maternal_network = ss.MaternalNet()
    if maternal_network is not None:
        networks.append(maternal_network)

    # setup demographics
    demographics = []
    if pregnancy == DEFAULT_STR:
        pregnancy = ss.Pregnancy(fertility_rate=10)
    if pregnancy is not None:
        demographics.append(pregnancy)

    if death == DEFAULT_STR:
        death = ss.Deaths(death_rate=10)
    if death is not None:
        demographics.append(death)

    # interventions and analyzers
    if interventions is None:
        interventions = []

    if analyzers is None:
        analyzers = []

    # build time
    sim = sti.Sim(
        start=1990,
        dur=duration,
        n_agents=n_agents,
        diseases=diseases,
        networks=networks,
        demographics=demographics,
        interventions=interventions,
        analyzers=analyzers
    )
    return sim
