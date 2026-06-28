"""Single source of truth for what the example test suite regression-checks.

Historically ``Molecule.check_ref`` compared *every* key in a reference JSON
except a hand-maintained ``skip_keys`` blocklist. That silently failed to test
quantities that live in sidecar files (IR/Raman intensities) or that were only
stored under internal ``OQP::`` arrays (non-adiabatic couplings), and nothing
flagged a newly added example whose reference was missing a value it should
carry.

This module replaces that with an explicit *allowlist* registry. Each entry
declares a physics quantity, where it lives, the runtypes/methods/properties it
applies to, and whether it is required (must be present and non-empty). The same
registry drives three things so they cannot drift apart:

  * ``keys_to_compare``  -> what ``check_ref`` compares
  * ``lean_keep``        -> what ``save_data(lean=True)`` keeps
  * ``missing_required`` -> the validation gate run when examples are added

Adding a new regression value is a one-line entry here.
"""

import glob
import json
import os
import re
from dataclasses import dataclass, field
from typing import FrozenSet, Optional

# ---------------------------------------------------------------------------
# Runtype groupings
# ---------------------------------------------------------------------------

# Runtypes that produce a meaningful nuclear gradient to regression-test.
#
# 'irc' is deliberately EXCLUDED: an IRC stops after a fixed number of steps at
# a NON-stationary point along the reaction path, whose geometry (and therefore
# gradient) is history-dependent and amplifies ULP-level BLAS/compiler
# differences into ~1e-4 cross-platform drift. A geometry *optimization*
# converges to a stationary point and stays reproducible, so those keep coord.
# For IRC we regression-test the reproducible quantity -- the energy.
GRAD_RUNTYPES = frozenset({
    'grad', 'optimize', 'meci', 'mecp', 'mep', 'tci', 'ts', 'neb',
})

# Runtypes whose final geometry is NOT a reliable cross-platform regression
# target (see GRAD_RUNTYPES note). 'coord' is skipped for these.
GEOMETRY_UNSTABLE_RUNTYPES = frozenset({'irc'})

# Keys that identify the structure / are bookkeeping. ``atoms`` is an exact
# identity check; the rest are metadata.
IDENTITY_KEYS = frozenset({'atoms', 'coord'})
METADATA_KEYS = frozenset({
    'json', 'symmetry_metadata', 'mass', 'hessian_metadata',
    'vibrational_intensity_metadata',
})


@dataclass(frozen=True)
class RegKey:
    """One regression quantity."""
    key: str
    # Runtypes this key applies to. '*' means every runtype.
    runtypes: object = '*'
    required: bool = False           # gate: must be present AND non-empty
    source: str = 'primary'          # 'primary' or 'sidecar'
    sidecar_field: Optional[str] = None  # field name in the sidecar (defaults to key)
    phase_invariant: bool = False    # compare magnitudes (sign/phase ambiguous)
    skip_sub: tuple = ()             # dict sub-keys to ignore (e.g. EKT orbitals)
    rtol: float = 0.0                # relative tolerance (0 = exact abs compare)
    needs_excited: bool = False      # only when an excited-state method is active
    needs_prop: Optional[str] = None  # only when this scf_prop was requested
    exclude_runtypes: frozenset = frozenset()  # runtypes where this is NOT checked

    def field(self):
        """The field name to read from its source (sidecar field or the key)."""
        return self.sidecar_field or self.key

    def applies(self, runtype, excited, props):
        if self.runtypes != '*' and runtype not in self.runtypes:
            return False
        if runtype in self.exclude_runtypes:
            return False
        if self.needs_excited and not excited:
            return False
        if self.needs_prop is not None and self.needs_prop not in props:
            return False
        return True


# ---------------------------------------------------------------------------
# THE REGISTRY -- add a row to regression-test a new quantity.
#
# HOW TO ADD A NEW REGRESSION VALUE (the one place you touch):
#
#   1. Add a single ``RegKey(...)`` row below describing the quantity:
#        - key       : the JSON key it is written under
#        - runtypes  : which runtypes carry it ('*' for all, or a frozenset)
#        - required  : True if every applicable example MUST have it non-empty
#        - source    : 'sidecar' if it lives in <input>.hess.json, else 'primary'
#        - needs_excited / needs_prop : gate on excited-state method or scf_prop
#        - phase_invariant / skip_sub : for sign-ambiguous or structured values
#   2. If it is a newly computed quantity, emit it from Molecule.get_results()
#      under the same key.
#   3. Regenerate the affected references once (OQP_LEAN_JSON=1); lean_keep()
#      keeps the new key automatically.
#
#   From then on it is compared by check_ref, kept by lean dumps, and ENFORCED
#   by the validation gate: any example missing a required key fails loudly
#   (``openqp --validate_examples`` / the CI step), so a new example or a new
#   value cannot silently go untested.
# ---------------------------------------------------------------------------
REGISTRY = (
    # Identity / geometry. atoms is an exact match; coord is the (input or
    # optimized) geometry, skipped for IRC where the path point is unstable.
    RegKey('atoms', runtypes='*', required=True),
    RegKey('coord', runtypes='*', required=True,
           exclude_runtypes=GEOMETRY_UNSTABLE_RUNTYPES),
    RegKey('energy', runtypes='*', required=True),
    RegKey('grad', runtypes=GRAD_RUNTYPES, required=True),
    # Hessian + derived vibrational quantities live in the <input>.hess.json
    # sidecar. 'hess' is the raw matrix (mapped from sidecar 'hessian').
    RegKey('hess', runtypes=frozenset({'hess'}), required=True, source='sidecar',
           sidecar_field='hessian'),
    RegKey('freqs', runtypes=frozenset({'hess'}), required=True, source='sidecar'),
    # IR intensities and Raman activities are second-order response properties
    # (dipole / polarizability derivatives over the analytic Hessian). They are
    # far more sensitive to the SCF convergence *path* than the energy or the
    # Hessian itself: an SCF that converges to the same energy via a different
    # route (e.g. the coarse->fine XC grid ramp, integral screening, a different
    # guess/BLAS order) shifts the converged density at ~1e-6 and these
    # quantities amplify that to ~1e-4 (IR) / ~1e-3 (Raman) -- above the global
    # round(diff,4) ~5e-5 gate, while the energy stays bit-identical. Compare
    # them with a small relative tolerance (a genuine regression is orders of
    # magnitude larger). Same rationale as the SOC rtol below.
    RegKey('infrared_intensities', runtypes=frozenset({'hess'}),
           required=True, source='sidecar', rtol=1e-3),
    RegKey('raman_activities', runtypes=frozenset({'hess'}),
           required=True, source='sidecar', rtol=1e-4),
    # Excitation energies: meaningful only for excited-state methods; a ground
    # state run stores the placeholder [0].
    RegKey('td_energies', runtypes='*', required=True, needs_excited=True),
    # SOC matrix elements are large (10^4-10^5 cm^-1), so the default absolute
    # round(diff,4) gate (~5e-5) would demand ~10 significant figures -- tighter
    # than ULP-level BLAS/compiler/integral-screening noise (~1e-8 relative),
    # which makes the test fail on harmless numerical differences across builds
    # (Debug vs Release, OpenBLAS vs Accelerate, response integral cutoff). Use a
    # relative tolerance instead; a real SOC regression is >>1e-6 relative, and
    # the response itself stays exact-checked via td_*_energies below.
    RegKey('soc', runtypes=frozenset({'soc'}), required=True, rtol=1e-6),
    # A SOC run computes both spin-resolved excitation ladders; the public
    # td_energies mirrors only one, so test the singlet and triplet explicitly.
    RegKey('td_singlet_energies', runtypes=frozenset({'soc'}), required=True),
    RegKey('td_triplet_energies', runtypes=frozenset({'soc'}), required=True),
    # Non-adiabatic / derivative couplings are sign/phase ambiguous between
    # runs, so compare magnitudes.
    RegKey('nac', runtypes=frozenset({'nac', 'nacme'}), required=True,
           phase_invariant=True),
    RegKey('mrsf_ekt', runtypes=frozenset({'ekt'}), required=True,
           skip_sub=('orbitals_mo', 'dyson_orbitals_mo')),
    # GIAO NMR shielding is a second-order response property and, like IR/Raman
    # above, is sensitive to the SCF convergence path well beyond the energy;
    # compare with a small relative tolerance (a real regression is far larger).
    RegKey('nmr_shielding', runtypes='*', required=True, needs_prop='nmr', rtol=1e-4),
    # SCF property results, each gated on its requested scf_prop value.
    RegKey('dipole', runtypes='*', required=True, needs_prop='el_mom'),
    RegKey('mulliken_charges', runtypes='*', required=True, needs_prop='mulliken'),
    RegKey('lowdin_charges', runtypes='*', required=True, needs_prop='lowdin'),
    RegKey('resp_charges', runtypes='*', required=True, needs_prop='resp'),
)

_BY_KEY = {e.key: e for e in REGISTRY}


def _ctx(runtype, excited, props):
    props = frozenset(props or ())
    return runtype, bool(excited), props


def keys_to_compare(runtype, excited=False, props=None):
    """Registry keys that ``check_ref`` should compare in this context."""
    rt, ex, pr = _ctx(runtype, excited, props)
    return [e for e in REGISTRY if e.applies(rt, ex, pr)]


def missing_required(present_nonempty, runtype, excited=False, props=None):
    """Required registry keys absent or empty in this reference -> gate failures.

    ``present_nonempty`` is the set of keys that exist and are non-empty in the
    reference (primary + sidecar merged).
    """
    rt, ex, pr = _ctx(runtype, excited, props)
    return [e.key for e in REGISTRY
            if e.required and e.applies(rt, ex, pr)
            and e.key not in present_nonempty]


def lean_keep(key):
    """True if a JSON key survives a lean (test-reference) dump.

    Keeps identity, metadata, and every registered physics key; drops internal
    ``OQP::`` arrays and anything else not declared a regression target.
    """
    return key in IDENTITY_KEYS or key in METADATA_KEYS or key in _BY_KEY


def is_phase_invariant(key):
    e = _BY_KEY.get(key)
    return bool(e and e.phase_invariant)


def skip_sub_keys(key):
    e = _BY_KEY.get(key)
    return e.skip_sub if e else ()


def sidecar_keys(runtype, excited=False, props=None):
    """Registry keys sourced from the .hess.json sidecar in this context."""
    return [e for e in keys_to_compare(runtype, excited, props)
            if e.source == 'sidecar']


# ---------------------------------------------------------------------------
# Validation gate -- run when examples are added (openqp --validate_examples).
# ---------------------------------------------------------------------------

def _safe_load(path):
    try:
        with open(path, 'r') as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return {}


def _context_from_input(inp_path):
    """Lightweight (runtype, excited, props) read straight from a .inp file."""
    try:
        text = open(inp_path, 'r').read()
    except OSError:
        return 'energy', False, []
    m = re.search(r'^\s*runtype\s*=\s*(\w+)', text, re.I | re.M)
    runtype = m.group(1).lower() if m else 'energy'
    excited = re.search(r'^\s*\[tdhf\]', text, re.I | re.M) is not None
    props = []
    pm = re.search(r'^\s*scf_prop\s*=\s*([^\n#]+)', text, re.I | re.M)
    if pm:
        props = [p.strip().lower() for p in re.split(r'[,\s]+', pm.group(1)) if p.strip()]
    return runtype, excited, props


def _present_nonempty(ref_path, inp_path):
    """Registry keys that are present and non-empty in a reference (primary +
    sidecar), keyed by registry key name."""
    primary = _safe_load(ref_path)
    sidecar = _safe_load(inp_path[:-4] + '.hess.json') if inp_path.endswith('.inp') else {}
    present = set()
    for e in REGISTRY:
        src = sidecar if e.source == 'sidecar' else primary
        v = src.get(e.field())
        if v not in (None, [], {}, ''):
            present.add(e.key)
    return present


def validate_examples(examples_dir):
    """Validate every example reference under ``examples_dir``.

    Returns a list of ``(relative_inp_path, [missing_required_keys])`` for
    references that lack a required regression value for their runtype/method/
    properties. An empty list means every example carries everything the
    registry says it must -- this is the gate that fails when a new example (or
    a newly registered value) is added without its reference value.
    """
    failures = []
    for inp in sorted(glob.glob(os.path.join(examples_dir, '**', '*.inp'),
                                recursive=True)):
        ref = inp[:-4] + '.json'
        if not os.path.exists(ref):
            continue  # not a regression example (no committed reference)
        runtype, excited, props = _context_from_input(inp)
        miss = missing_required(_present_nonempty(ref, inp), runtype, excited, props)
        if miss:
            failures.append((os.path.relpath(inp, examples_dir), miss))
    return failures


# ---------------------------------------------------------------------------
# Feature-coverage gate -- every opt-in capability must be exercised by an
# example, so a new feature cannot ship without a test.
#
# Mechanism: a feature almost always adds a boolean keyword to
# OQP_CONFIG_SCHEMA (default False = opt-in). This gate requires every such
# flag to be set true by at least one committed example .inp, UNLESS it is:
#   * EXEMPT_FLAGS      -- not a distinct testable capability (IO/restart/
#                          scratch-cleanup/debug/test-harness/placeholder).
#   * KNOWN_UNCOVERED   -- a real capability that is grandfathered as a tracked
#                          gap (TODO: add an example). New flags may NOT be
#                          added here without review.
# A new opt-in flag that is neither exercised nor listed fails the gate.
# Classifications were derived by auditing each flag's code usage.
# ---------------------------------------------------------------------------

# section.option -> reason. Not a distinct testable capability.
EXEMPT_FLAGS = {
    'guess.continue_geom': 'JSON-restart geometry-selection IO convenience',
    'symmetry.strict': 'validation guard (mismatch -> fatal error), not a capability',
    'scf.trh_ls': 'TRAH line-search numerical sub-knob',
    'properties.td_prop': 'unfinished placeholder; input checker warns to disable',
    'properties.export': 'writes already-computed data to text (pure IO toggle)',
    'properties.back_door': 'programmatic data-injection hook for external drivers',
    'hess.read': 'reload cached .hess.json instead of recomputing (IO restart)',
    'hess.restart': 'resume numerical Hessian from cached scratch (IO restart)',
    'hess.clean': 'rmtree numerical-Hessian scratch after run (IO cleanup)',
    'nac.restart': 'reload cached .dcme scratch in numerical NAC (IO restart)',
    'nac.clean': 'rmtree numerical-NAC scratch after run (IO cleanup)',
    'tests.exception': 'test-harness toggle (raise vs exit on non-convergence)',
}

# section.option -> reason. Real capability, tracked gap (needs an example).
KNOWN_UNCOVERED = {
    'input.d4': 'native DFT-D4 dispersion (PR #235); add a d4=true example',
    'symmetry.use_response_symmetry': 'irrep-blocked Davidson TDDFT response solver',
    'scf.trh_stab': 'TRAH SCF stability-following (escape unstable solutions)',
    'optimize.init_scf': 'fresh initial-guess SCF each optimization step',
    'nac.bp': 'branching-plane analysis at conical intersections',
}

_SCHEMA_FILE = os.path.join(os.path.dirname(__file__), os.pardir,
                            'molecule', 'oqpdata.py')


def _optin_bool_flags(schema_file=None):
    """{`section.option`: default} for every bool schema keyword default False."""
    import ast
    src = schema_file or _SCHEMA_FILE
    tree = ast.parse(open(src).read())
    node = next((n.value for n in ast.walk(tree)
                 if isinstance(n, ast.Assign)
                 and any(getattr(t, 'id', '') == 'OQP_CONFIG_SCHEMA'
                         for t in n.targets)), None)
    out = {}
    if not isinstance(node, ast.Dict):
        return out
    for sk, sv in zip(node.keys, node.values):
        section = ast.literal_eval(sk)
        if not isinstance(sv, ast.Dict):
            continue
        for ok, ov in zip(sv.keys, sv.values):
            option = ast.literal_eval(ok)
            spec = {}
            for k, v in zip(ov.keys, ov.values):
                key = ast.literal_eval(k)
                if isinstance(v, ast.Name):
                    spec[key] = v.id
                elif isinstance(v, ast.Constant):
                    spec[key] = v.value
            if spec.get('type') == 'bool' and str(spec.get('default')).lower() in ('false',):
                out['%s.%s' % (section, option)] = spec.get('default')
    return out


def _flag_exercised_true(option, examples_dir):
    """True if any example .inp sets ``option`` to a truthy value."""
    pat = re.compile(r'^\s*' + re.escape(option) + r'\s*=\s*'
                     r'(true|t|1|yes|\.true\.)\s*(#.*)?$', re.I | re.M)
    for inp in glob.glob(os.path.join(examples_dir, '**', '*.inp'), recursive=True):
        try:
            if pat.search(open(inp).read()):
                return True
        except OSError:
            continue
    return False


def feature_coverage(examples_dir, schema_file=None):
    """Check that every opt-in feature flag is exercised by an example.

    Returns ``(failures, grandfathered)``:
      * ``failures``     -- [(flag, 'no example exercises it')] for flags that
                            are neither exempt nor known-uncovered -> gate fails.
      * ``grandfathered``-- [(flag, reason)] tracked gaps that still need a test.
    """
    flags = _optin_bool_flags(schema_file)
    failures, grandfathered = [], []
    for flag in sorted(flags):
        _, option = flag.split('.', 1)
        if flag in EXEMPT_FLAGS:
            continue
        if _flag_exercised_true(option, examples_dir):
            continue
        if flag in KNOWN_UNCOVERED:
            grandfathered.append((flag, KNOWN_UNCOVERED[flag]))
        else:
            failures.append((flag, 'opt-in feature flag not exercised by any '
                                   'example (add an example, or classify it in '
                                   'EXEMPT_FLAGS / KNOWN_UNCOVERED with a reason)'))
    return failures, grandfathered
