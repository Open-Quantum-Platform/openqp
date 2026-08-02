"""AFQMC option marshalling.

The AFQMC compute path lives entirely in Fortran (``source/afqmc/``); Python's
only job is to turn the ``[afqmc]`` input section into the packed option vector
the entry point reads, and to stage the optional trial files under the fixed
names the vendored readers expect.

The option layout below is a mirror of the one declared in
``source/afqmc/afqmc_oqp.F90``. The two must be changed together; the length is
checked on the Fortran side.
"""

import os
import shutil

import numpy as np

#: ``[afqmc] trial`` values and the codes the Fortran side branches on.
TRIAL_KINDS = {
    'mean_field': 0,   # single closed/open-shell SCF determinant
    'sf': 1,           # SF-CIS determinant expansion, from trial_file
    'mrsf': 2,         # MRSF-CIS determinant expansion, from trial_file
    'mrsf_state': 3,   # an MRSF-TDDFT root computed by OpenQP itself, in memory
}

#: Trial kinds that consume an externally supplied determinant file.
TRIAL_FILES = {1: 'SFDAT', 2: 'MRSFDAT'}

#: Order of the packed real64 option vector. Mirrors afqmc_oqp.F90.
AFQMC_OPTION_ORDER = (
    'trial',                     # 1  trial kind code
    'chol_tol',                  # 2
    'nwalkers',                  # 3
    'nsteps',                    # 4
    'timestep',                  # 5
    'seed',                      # 6
    'stabilize_every',           # 7
    'population_control_every',  # 8
    'estimate_every',            # 9
    'accumulate_after',          # 10
    'force_bias_bound',          # 11
    'oo_orbitals',               # 12
    'nlow',                      # 13
    'low_max',                   # 14
    'low_cap',                   # 15
    'low_start',                 # 16
    'state',                     # 17 1-based MRSF root for trial=mrsf_state
    'trial_threshold',           # 18 drop trial determinants below this |c|
)


def trial_kind(config):
    """Resolve ``[afqmc] trial`` to its integer code."""
    name = str(config.get('trial', 'mean_field')).strip().lower()
    try:
        return TRIAL_KINDS[name]
    except KeyError:
        known = ', '.join(sorted(TRIAL_KINDS))
        raise ValueError(f'Unknown [afqmc] trial {name!r}. Use one of: {known}') from None


def pack_afqmc_options(config):
    """Pack the ``[afqmc]`` config section into the Fortran option vector."""
    values = dict(config)
    values['trial'] = trial_kind(config)
    values['oo_orbitals'] = 1.0 if _as_bool(config.get('oo_orbitals', False)) else 0.0
    packed = np.zeros(len(AFQMC_OPTION_ORDER), dtype=np.float64)
    for index, key in enumerate(AFQMC_OPTION_ORDER):
        packed[index] = float(values[key])
    return packed


def stage_trial_files(config):
    """Copy user-supplied trial/orbital/lower-state files to their fixed names.

    The vendored readers take 8-character file names from the parent program's
    input convention, so the driver stages whatever the user pointed at into the
    working directory under the expected name rather than plumbing paths through
    the Fortran interface.
    """
    kind = trial_kind(config)
    staged = {}

    target = TRIAL_FILES.get(kind)
    if target:
        source = str(config.get('trial_file', '')).strip()
        if not source:
            raise ValueError(f'[afqmc] trial = {config["trial"]} requires trial_file')
        staged[target] = _stage(source, target)

    if _as_bool(config.get('oo_orbitals', False)):
        source = str(config.get('oo_file', '')).strip()
        if not source:
            raise ValueError('[afqmc] oo_orbitals requires oo_file')
        staged['OOORBDAT'] = _stage(source, 'OOORBDAT')

    if int(config.get('nlow', 0)) > 0:
        source = str(config.get('low_file', '')).strip()
        if not source:
            raise ValueError('[afqmc] nlow requires low_file')
        staged['LOWSTATE'] = _stage(source, 'LOWSTATE')

    return staged


def _stage(source, target):
    if not os.path.isfile(source):
        raise FileNotFoundError(f'[afqmc] cannot read {source}')
    if os.path.abspath(source) != os.path.abspath(target):
        shutil.copyfile(source, target)
    return os.path.abspath(target)


def _as_bool(value):
    if isinstance(value, str):
        return value.strip().lower() in ('1', 'true', 'yes', 'on')
    return bool(value)
