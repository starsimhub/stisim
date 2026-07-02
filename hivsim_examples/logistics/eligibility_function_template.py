# template variables and their settings
# product_name: {product_name}
# product_quantity: {product_quantity}
# event_name_downcased: {event_name_downcased}
# link_value: {link_value}
# link_event_name_downcased: {link_event_name_downcased}
# bundle_name: {bundle_name}

def {event_name_downcased}_eligibilities(sim):
    bundle = getattr(sim.interventions, {bundle_name})
    selected = getattr(bundle, {link_event_name_downcased} + '_selected')
    results = getattr(bundle, {link_event_name_downcased} + '_results')
    link_value = {link_value}

    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities