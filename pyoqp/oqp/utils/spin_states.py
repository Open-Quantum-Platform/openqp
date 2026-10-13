"""Spin manifolds of a response calculation, read from the config alone.

A run can answer two spin ladders at once, so "which state" is only meaningful
together with "which ladder".  These readers are the single place that maps the
input keys onto that pair, shared by the drivers that solve the manifolds and the
input checker that judges the request.

Pure config access: no engine import, so the parser, the checker and the library
can all use it.
"""


def _states(value) -> list[int]:
    """A gradient selector as written in an input file, or an already parsed list."""
    if value is None or value == '':
        return []
    if isinstance(value, str):
        value = value.replace(',', ' ').split()
    return [int(state) for state in value]


def soc_requested(config) -> bool:
    """True when this run computes spin-orbit coupling as a step of the flow.

    Read tolerantly: scripting callers and partial configurations never pass
    through schema normalization, so the value can still be the input string.
    """
    value = (config.get('properties') or {}).get('soc', False)
    if isinstance(value, str):
        return value.strip().lower() in ('true', 't', '1', 'yes', 'on')
    return bool(value)


def configured_spin(config) -> str:
    """Spin label of the manifold this run answers on its own.

    Any multiplicity other than 1 or 3 reads as the singlet ladder; only a
    spin-orbit run is restricted, and the input checker says so there.  This
    reader must not raise, because every gradient run passes through it.
    """
    multiplicity = int((config.get('tdhf') or {}).get('multiplicity', 1))
    return 't' if multiplicity == 3 else 's'


def other_spin(spin_label: str) -> str:
    """The manifold a run answers in addition to the one it configured."""
    return 's' if spin_label == 't' else 't'


def grad_requests(config) -> dict:
    """Gradient request per spin manifold, falling back to [properties] grad.

    grad_s and grad_t name a ladder outright.  The historical selector describes
    the ladder this run already answers, so it fills only that one when nothing
    named it; that is what keeps every existing deck asking for the same states.
    """
    properties = config.get('properties') or {}
    requests = {label: _states(properties.get(f'grad_{label}', []))
                for label in ('s', 't')}
    configured = configured_spin(config)
    if not requests[configured]:
        requests[configured] = _states(properties.get('grad', []))
    return requests
