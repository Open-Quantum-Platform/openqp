from typing import Dict, Tuple
from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
from oqp.utils.input_parser import SECTION_OPTION_ALIASES

def build_kw_map(schema: dict) -> Dict[str, Tuple[str, str]]:
    """
    Build a flat mapping { 'section.option': (section, option) } from OQP_CONFIG_SCHEMA.
    """
    flat: Dict[str, Tuple[str, str]] = {}
    for section, options in schema.items():
        for opt in options.keys():
            flat[f"{section}.{opt}"] = (section, opt)
    return flat

KW_MAP_PURE = build_kw_map(OQP_CONFIG_SCHEMA)


def canonical_param_key(user_key: str) -> str:
    """Return the schema key for a public key or a retained legacy alias."""
    key = str(user_key)
    if "." not in key:
        return key
    section, option = key.split(".", 1)
    option = SECTION_OPTION_ALIASES.get(section, {}).get(option, option)
    return f"{section}.{option}"

def resolve_param_key(user_key: str) -> Tuple[str, str]:
    """
    """
    user_key = canonical_param_key(user_key)
    if user_key in KW_MAP_PURE:
        return KW_MAP_PURE[user_key]
    raise KeyError(f"Unknown parameter '{user_key}'. "
                   f"Valid keys are: {', '.join(sorted(KW_MAP_PURE.keys()))}")
