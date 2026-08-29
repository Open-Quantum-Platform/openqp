"""Production analytic MRSF-TDDFT NAC -- the nac-lagrangian assembly (v3).

For every ordered pair (I != J), all terms from certified components (see
tools/nac_lagrangian/MRSF_NAC_DERIVATION.md, Secs. 4 and 7.24-7.41):

  d_IJ = antisym[ T1 + zeta(X):B^R + X:V + gamma:Sk ]
    T1    = [amp2e + esum](ytil_IJ, X_J)      slot-injected engines
    X     = MT_frozen + MT_response + gamma
    Z:B^R = one native ROHF/ROKS adjoint solve and analytic Fortran contraction
    X:V   = resident symmetric-U / overlap-derivative contraction
    gamma = exact-tlf state-metric derivative built in resident Fortran
    ytil  = X_I / (om_J - om_I)                 orthonormal-eigenvector identity

MT_frozen comes from the resident Fortran ``mrsf_nac_wpair`` closed-form
bilinear adjoint.  The single ``mrsf_nac_lagrangian`` call owns the actual
state count, streamed exact metric columns, eigenvector response, ordered-pair
loop, source injection, ROHF/ROKS Z solve, every nuclear-coordinate
contraction, antisymmetrization, and gap scaling.  Python only validates the
public calculation scope and reshapes the final records.  The computational
adjoint zeta is minus the Lee-gradient multiplier convention.
"""
import os
from functools import wraps

import numpy as np


_MISSING = object()
_NAC_MUTATED_ENV = ('NAC_DUMP_ROHF_RESPONSE',)


def _with_temporary_nac_state(function):
    """Restore process and molecule state changed by an analytic NAC call."""
    @wraps(function)
    def wrapped(mol, *args, **kwargs):
        environment = {
            name: os.environ[name] if name in os.environ else _MISSING
            for name in _NAC_MUTATED_ENV
        }
        try:
            td_bvec_mo = np.array(mol.data['OQP::td_bvec_mo'], copy=True)
        except Exception:
            td_bvec_mo = _MISSING
        try:
            control = mol.data._data.control
            int2e_cutoff = control.int2e_cutoff
        except Exception:
            control = None
            int2e_cutoff = _MISSING

        try:
            return function(mol, *args, **kwargs)
        finally:
            # Keep the cleanup nested so later state is restored even if an
            # earlier assignment unexpectedly fails.
            try:
                if td_bvec_mo is not _MISSING:
                    mol.data['OQP::td_bvec_mo'] = td_bvec_mo
            finally:
                try:
                    if control is not None and int2e_cutoff is not _MISSING:
                        control.int2e_cutoff = int2e_cutoff
                finally:
                    for name, old_value in environment.items():
                        if old_value is _MISSING:
                            os.environ.pop(name, None)
                        else:
                            os.environ[name] = old_value

    return wrapped


def _resident_pair_cartesian(raw, nstate, natom):
    """Expose a Fortran (3*natom,nstate,nstate) tensor as [I,J,atom,xyz]."""
    flat = np.asarray(raw, dtype=float).reshape(-1)
    expected = 3 * natom * nstate * nstate
    if flat.size != expected:
        raise RuntimeError(
            'Resident MRSF NAC pair tensor has an inconsistent size: '
            f'{flat.size} != {expected}'
        )
    # Fortran stores the first state index faster than the second.  A C-order
    # view therefore exposes [J,I,atom,xyz]; the transpose restores [I,J,...].
    return flat.reshape(nstate, nstate, natom, 3).transpose(1, 0, 2, 3).copy()


@_with_temporary_nac_state
def analytic_nac(mol):
    import oqp

    debug_path = os.environ.get('NAC_ANALYTIC_DEBUG')

    if mol.config['tdhf']['multiplicity'] != 1:
        raise NotImplementedError(
            'analytic MRSF NAC v3 currently implements the singlet fold only'
        )

    # The ROHF/ROKS stationarity residual is amplified by the near-degenerate
    # NAC response.  A 1e-6 SCF threshold produces millihartree/bohr-scale
    # coupling errors even though the Z equation itself is converged; the
    # same-process forward-response gate closes at 1e-5 only once the SCF is
    # converged to at least 1e-8.  Refuse a silently inaccurate result.
    if float(mol.config['scf']['conv']) > 1.0e-8:
        raise RuntimeError(
            'Analytic MRSF NAC requires [scf] conv <= 1e-8; '
            '1e-10 is recommended near a crossing.'
        )
    if float(mol.config['tdhf']['conv']) > 1.0e-8:
        raise RuntimeError(
            'Analytic MRSF NAC requires [tdhf] conv <= 1e-8; '
            '1e-10 is recommended near a crossing.'
        )

    # Own the record before the resident driver reserves/removes TagArray
    # outputs; a NumPy view into TagArray storage would otherwise be invalid.
    td_energies = np.array(
        mol.data['OQP::td_energies'], dtype=float, copy=True
    ).reshape(-1)
    nstate = td_energies.size
    if nstate < 2:
        raise RuntimeError('Analytic MRSF NAC requires at least two states')
    resident_nstate = int(mol.data._data.tddft.nstate)
    if resident_nstate != nstate:
        raise RuntimeError(
            'Resident MRSF state count and energy record disagree: '
            f'{resident_nstate} != {nstate}'
        )
    natom = int(mol.data['natom'])

    # Production has no nuclear-coordinate forward CPHF branch.  The separate
    # tools/nac_lagrangian/rohf_response_gate.py diagnostic owns that 3N solve.
    oqp.mrsf_nac_lagrangian(mol)

    if debug_path:
        # Debug export is deliberately observational: it copies resident
        # records after the single Fortran call and performs no pair algebra.
        debug = {
            'td_energies': td_energies.copy(),
            'td_bvec_mo': np.array(
                mol.data['OQP::td_bvec_mo'], copy=True
            ),
            'dp_ordered': np.array(
                mol.data['OQP::nac_dp_ordered'], copy=True
            ),
            'dcv': np.array(mol.data['OQP::nac_dcv'], copy=True),
            'nacv': np.array(mol.data['OQP::nac_nacv'], copy=True),
        }
        for key, tag in (
                ('predictor_dcv', 'OQP::nac_predictor_dcv'),
                ('predictor_nacv', 'OQP::nac_predictor_nacv')):
            try:
                debug[key] = np.array(mol.data[tag], copy=True)
            except (KeyError, RuntimeError, TypeError, ValueError):
                pass
        np.savez(debug_path, **debug)

    dcv = _resident_pair_cartesian(
        mol.data['OQP::nac_dcv'], nstate, natom
    )
    nacv = _resident_pair_cartesian(
        mol.data['OQP::nac_nacv'], nstate, natom
    )
    return nacv, dcv
