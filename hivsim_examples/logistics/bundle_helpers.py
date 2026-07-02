_VALUE_VARS = frozenset({'product_name', 'product_quantity', 'link_value', 'bundle_name', 'link_event_name_downcased'})
_IDENTIFIER_VARS = frozenset({'event_name_downcased'})


def fill_template(template_text: str, **kwargs) -> str:
    """Fill a bundle code template.

    Value variables (product_name, product_quantity, link_value, bundle_name,
    link_event_name_downcased) are passed through repr() so they render as Python
    literals.  The one identifier variable (event_name_downcased) is inserted as-is
    because it appears as a function/variable name in the eligibility template.
    """
    processed = {}
    for key, value in kwargs.items():
        if key in _VALUE_VARS:
            processed[key] = repr(value)
        elif key in _IDENTIFIER_VARS:
            processed[key] = value
        else:
            raise ValueError(f"Unknown template variable: {key!r}")
    return template_text.format(**processed)
