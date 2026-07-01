import numpy as np
import starsim as ss

from stisim.intervention_bundle import InterventionBundle
from stisim.logistics import Product, ProductCategory, DeliveryMode, Supply, Supplies

from hivsim_examples.logistics.clinic.enter_clinic import EnterClinic
from hivsim_examples.logistics.clinic.queue_for_hiv_test import QueueForHivTest
from hivsim_examples.logistics.clinic.take_hiv_test import TakeHivTest
from hivsim_examples.logistics.clinic.art_uptake_choice import ArtUptakeChoice
from hivsim_examples.logistics.clinic.initiate_art import InitiateArt
from hivsim_examples.logistics.clinic.set_art_return import SetArtReturn
from hivsim_examples.logistics.clinic.prep_uptake_choice import PrepUptakeChoice
from hivsim_examples.logistics.clinic.initiate_prep import InitiatePrep
from hivsim_examples.logistics.clinic.set_prep_return import SetPrepReturn
from hivsim_examples.logistics.clinic.condoms_uptake_choice import CondomsUptakeChoice
from hivsim_examples.logistics.clinic.initiate_condoms import InitiateCondoms
from hivsim_examples.logistics.clinic.set_condoms_return import SetCondomsReturn


# ---------------------------------------------------------------------------
# Products
# ---------------------------------------------------------------------------

hiv_test_product = Product(
    name='hiv_test',
    category=ProductCategory.DIAGNOSTIC,
    delivery_mode=DeliveryMode.RAPID_TEST,
    cost=0,  # User note: set the actual cost
    eff_by_ti=[1.0],  # User note: set the actual eff_by_ti values
)

dolutegravir_product = Product(
    name='dolutegravir',
    category=ProductCategory.ART,
    delivery_mode=DeliveryMode.PILL,
    cost=0,  # User note: set the actual cost
    eff_by_ti=[1.0],  # User note: set the actual eff_by_ti values
)

oral_product = Product(
    name='oral',
    category=ProductCategory.PREP,
    delivery_mode=DeliveryMode.PILL,
    cost=0,  # User note: set the actual cost
    eff_by_ti=[1.0],  # User note: set the actual eff_by_ti values
)

condoms_product = Product(
    name='condoms',
    category=ProductCategory.BARRIER,
    delivery_mode=DeliveryMode.GIVE,
    cost=0,  # User note: set the actual cost
    eff_by_ti=[1.0],  # User note: set the actual eff_by_ti values
)


# ---------------------------------------------------------------------------
# Supplies
# ---------------------------------------------------------------------------

supplies = Supplies(supplies=[
    Supply(quantity=np.inf, product=hiv_test_product),
    Supply(quantity=np.inf, product=dolutegravir_product),
    Supply(quantity=np.inf, product=oral_product),
    Supply(quantity=np.inf, product=condoms_product),
])


# ---------------------------------------------------------------------------
# Event 1: enter_clinic
# ---------------------------------------------------------------------------

# User note: adjust as needed
def enter_clinic_eligibility(sim):
    is_female = sim.people.female
    in_age_range = (sim.people.age >= 15) & (sim.people.age < 50)
    not_hiv_diagnosed = ~sim.diseases.hiv.diagnosed
    return is_female & in_age_range & not_hiv_diagnosed

enter_clinic = EnterClinic(
    name='enter_clinic',
    eligibilities=[enter_clinic_eligibility],
    coverages=[0.25],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 2: queue_for_hiv_test
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: queue_for_hiv_test
# link_value: True
# link_event_name_downcased: enter_clinic
# bundle_name: clinic

def queue_for_hiv_test_eligibilities(sim):
    selected = sim.interventions.clinic.enter_clinic_selected
    results = sim.interventions.clinic.enter_clinic_results
    link_value = True
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

queue_for_hiv_test = QueueForHivTest(
    name='queue_for_hiv_test',
    eligibilities=[queue_for_hiv_test_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 3: take_hiv_test
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: take_hiv_test
# link_value: True
# link_event_name_downcased: queue_for_hiv_test
# bundle_name: clinic

def take_hiv_test_eligibilities(sim):
    selected = sim.interventions.clinic.queue_for_hiv_test_selected
    results = sim.interventions.clinic.queue_for_hiv_test_results
    link_value = True
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

take_hiv_test = TakeHivTest(
    name='take_hiv_test',
    eligibilities=[take_hiv_test_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 4: art_uptake_choice
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: art_uptake_choice
# link_value: True
# link_event_name_downcased: take_hiv_test
# bundle_name: clinic

def art_uptake_choice_eligibilities(sim):
    selected = sim.interventions.clinic.take_hiv_test_selected
    results = sim.interventions.clinic.take_hiv_test_results
    link_value = True
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

art_uptake_choice = ArtUptakeChoice(
    name='art_uptake_choice',
    eligibilities=[art_uptake_choice_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 5: initiate_art
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: 'dolutegravir'
# product_quantity: 1
# event_name_downcased: initiate_art
# link_value: True
# link_event_name_downcased: art_uptake_choice
# bundle_name: clinic

def initiate_art_eligibilities(sim):
    selected = sim.interventions.clinic.art_uptake_choice_selected
    results = sim.interventions.clinic.art_uptake_choice_results
    link_value = True
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

initiate_art = InitiateArt(
    name='initiate_art',
    eligibilities=[initiate_art_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 6: set_art_return
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: set_art_return
# link_value: 'all'
# link_event_name_downcased: initiate_art
# bundle_name: clinic

def set_art_return_eligibilities(sim):
    selected = sim.interventions.clinic.initiate_art_selected
    results = sim.interventions.clinic.initiate_art_results
    link_value = 'all'
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

set_art_return = SetArtReturn(
    name='set_art_return',
    eligibilities=[set_art_return_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 7: prep_uptake_choice
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: prep_uptake_choice
# link_value: False
# link_event_name_downcased: take_hiv_test
# bundle_name: clinic

def prep_uptake_choice_eligibilities(sim):
    selected = sim.interventions.clinic.take_hiv_test_selected
    results = sim.interventions.clinic.take_hiv_test_results
    link_value = False
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

prep_uptake_choice = PrepUptakeChoice(
    name='prep_uptake_choice',
    eligibilities=[prep_uptake_choice_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 8: initiate_prep
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: 'oral'
# product_quantity: 1
# event_name_downcased: initiate_prep
# link_value: True
# link_event_name_downcased: prep_uptake_choice
# bundle_name: clinic

def initiate_prep_eligibilities(sim):
    selected = sim.interventions.clinic.prep_uptake_choice_selected
    results = sim.interventions.clinic.prep_uptake_choice_results
    link_value = True
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

initiate_prep = InitiatePrep(
    name='initiate_prep',
    eligibilities=[initiate_prep_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 9: set_prep_return
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: set_prep_return
# link_value: 'all'
# link_event_name_downcased: initiate_prep
# bundle_name: clinic

def set_prep_return_eligibilities(sim):
    selected = sim.interventions.clinic.initiate_prep_selected
    results = sim.interventions.clinic.initiate_prep_results
    link_value = 'all'
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

set_prep_return = SetPrepReturn(
    name='set_prep_return',
    eligibilities=[set_prep_return_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 10: condoms_uptake_choice
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: condoms_uptake_choice
# link_value: False
# link_event_name_downcased: prep_uptake_choice
# bundle_name: clinic

def condoms_uptake_choice_eligibilities(sim):
    selected = sim.interventions.clinic.prep_uptake_choice_selected
    results = sim.interventions.clinic.prep_uptake_choice_results
    link_value = False
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

condoms_uptake_choice = CondomsUptakeChoice(
    name='condoms_uptake_choice',
    eligibilities=[condoms_uptake_choice_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 11: initiate_condoms
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: 'condoms'
# product_quantity: 25
# event_name_downcased: initiate_condoms
# link_value: True
# link_event_name_downcased: condoms_uptake_choice
# bundle_name: clinic

def initiate_condoms_eligibilities(sim):
    selected = sim.interventions.clinic.condoms_uptake_choice_selected
    results = sim.interventions.clinic.condoms_uptake_choice_results
    link_value = True
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

initiate_condoms = InitiateCondoms(
    name='initiate_condoms',
    eligibilities=[initiate_condoms_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# Event 12: set_condoms_return
# ---------------------------------------------------------------------------

# template variables and their settings
# product_name: None
# product_quantity: 0
# event_name_downcased: set_condoms_return
# link_value: 'all'
# link_event_name_downcased: initiate_condoms
# bundle_name: clinic

def set_condoms_return_eligibilities(sim):
    selected = sim.interventions.clinic.initiate_condoms_selected
    results = sim.interventions.clinic.initiate_condoms_results
    link_value = 'all'
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities

set_condoms_return = SetCondomsReturn(
    name='set_condoms_return',
    eligibilities=[set_condoms_return_eligibilities],
    supplies=supplies,
)


# ---------------------------------------------------------------------------
# InterventionBundle states
# ---------------------------------------------------------------------------

# User note: default type is ss.BoolState; update types as needed for each state
states = [
    ss.BoolState('enter_clinic_selected'),
    ss.BoolState('enter_clinic_results'),
    ss.BoolState('queue_for_hiv_test_selected'),
    ss.BoolState('queue_for_hiv_test_results'),
    ss.BoolState('take_hiv_test_selected'),
    ss.BoolState('take_hiv_test_results'),
    ss.BoolState('art_uptake_choice_selected'),
    ss.BoolState('art_uptake_choice_results'),
    ss.BoolState('initiate_art_selected'),
    ss.BoolState('initiate_art_results'),
    ss.BoolState('set_art_return_selected'),
    ss.BoolState('set_art_return_results'),
    ss.BoolState('prep_uptake_choice_selected'),
    ss.BoolState('prep_uptake_choice_results'),
    ss.BoolState('initiate_prep_selected'),
    ss.BoolState('initiate_prep_results'),
    ss.BoolState('set_prep_return_selected'),
    ss.BoolState('set_prep_return_results'),
    ss.BoolState('condoms_uptake_choice_selected'),
    ss.BoolState('condoms_uptake_choice_results'),
    ss.BoolState('initiate_condoms_selected'),
    ss.BoolState('initiate_condoms_results'),
    ss.BoolState('set_condoms_return_selected'),
    ss.BoolState('set_condoms_return_results'),
]


# ---------------------------------------------------------------------------
# InterventionBundle instance
# ---------------------------------------------------------------------------

clinic = InterventionBundle(
    interventions=[
        enter_clinic,
        queue_for_hiv_test,
        take_hiv_test,
        art_uptake_choice,
        initiate_art,
        set_art_return,
        prep_uptake_choice,
        initiate_prep,
        set_prep_return,
        condoms_uptake_choice,
        initiate_condoms,
        set_condoms_return,
    ],
    states=states,
    name='clinic',
)
