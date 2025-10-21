from typing import Dict, Tuple
from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA

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

def resolve_param_key(user_key: str) -> Tuple[str, str]:
    """
    """
    if user_key in KW_MAP_PURE:
        return KW_MAP_PURE[user_key]
    raise KeyError(f"Unknown parameter '{user_key}'. "
                   f"Valid keys are: {', '.join(sorted(KW_MAP_PURE.keys()))}")
