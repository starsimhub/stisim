# template variables and their settings
# product_name: {product_name}
# product_quantity: {product_quantity}
# event_name_downcased: {event_name_downcased}
# link_value: {link_value}
# link_event_name_downcased: {link_event_name_downcased}
# bundle_name: {bundle_name}

def {event_name_downcased}_eligibilities(sim):
    selected = sim.interventions.{bundle_name}.{link_event_name_downcased}_selected  # TODO: update to use getattr to allow consistent use of repr() in .format() call
    results = sim.interventions.{bundle_name}.{link_event_name_downcased}_results  # TODO: update to use getattr to allow consistent use of repr() in .format() call
    link_value = {link_value}
    eligibilities = []
    for s, r in zip(selected, results):
        if link_value == "all":
            eligibilities.append(s)
        else:
            eligibilities.append(s[r == link_value])
    return eligibilities