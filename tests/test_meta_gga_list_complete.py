"""The input checker's meta-GGA name set must cover every meta-GGA OpenQP knows.

hess.type=auto must not pick the analytic Hessian (no tau channel) for a
meta-GGA.  The runtime asks LibXC for the functional family
(oqp_functional_needs_tau); the static pre-check in input_checker uses a name
set.  This test keeps that set complete: every functional name in
source/dftlib/libxc.F90 whose definition adds an XC_MGGA_* or XC_HYB_MGGA_*
LibXC component must be listed.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_every_libxc_meta_gga_name_is_in_the_checker_set():
    src = (ROOT / "source/dftlib/libxc.F90").read_text()
    checker = (ROOT / "pyoqp/oqp/utils/input_checker.py").read_text()
    names = set(re.search(r'_UMRSF_META_GGA_FUNCTIONALS = frozenset\("""(.*?)"""',
                          checker, re.S).group(1).split())
    missing = []
    for part in re.split(r"\n\s*case\s*\(", src)[1:]:
        label = re.match(r"([^)]*)\)", part)
        labels = [x.lower() for x in re.findall(r'"([^"]+)"', label.group(1))] if label else []
        body = part.split("\n    case", 1)[0]
        if labels and re.search(r"XC_(HYB_)?MGGA", body):
            missing += [name for name in labels if name not in names]
    assert not missing, f"meta-GGA functionals missing from _UMRSF_META_GGA_FUNCTIONALS: {sorted(set(missing))}"
