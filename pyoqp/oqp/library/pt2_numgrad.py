"""Central-difference numerical gradients for the native PT2 family.

The native PT2 family (method = caspt2 / ms-caspt2 / xms-caspt2 / mrmp2 /
mcqdpt2 / xmcqdpt2, all dispatched through
``oqp.library.caspt2_dyall.native_caspt2_energy``) is energy-only: there are
no analytic PT2 nuclear gradients.  This module supplies finite-difference
gradients by central differences over Cartesian displacements, re-running the
full SCF + (CASSCF/CASCI reference) + PT2 pipeline at each displaced geometry
through the exact same code path a single-point energy uses
(``SinglePoint.energy``).  ``Gradient.gradient`` in
``oqp.library.single_point`` dispatches here for the PT2 methods, which makes
runtype=grad and the gradient-driven runtypes (optimize, meci with the
penalty/ubp searches, mecp, ts, mep, neb, irc) work unchanged for PT2.

Conventions
-----------
* ``mol.energies`` for PT2 holds the PT2 total energies of the computed roots
  in ascending order (single-state methods: one entry).  Gradient state
  indices (``[properties] grad``, ``[optimize] istate/jstate``) index this
  list directly, exactly like the hf/tdhf convention.  Which roots are
  computed (and therefore what the list contains) follows the existing
  ``[pt2] root`` / ``target_roots`` / ``nroot`` convention.
* Coordinates are Bohr; energies Hartree; gradients Hartree/Bohr, shaped
  ``(nstate, natom, 3)`` like the hf/tdhf gradient arrays.

Configuration ([pt2] section, read with dict.get so missing schema rows are
harmless; see PT2_GRAD_NOTES.md for the schema/checker rows to add)
--------------------------------------------------------------------------
* ``grad_step``  (float, default 1.0e-3 Bohr): central-difference half-step.
  The energy pipeline is converged to ~1e-9..1e-10 Eh and bitwise
  deterministic for a fixed guess, so the FD noise floor is ~eps/h ~ 1e-6
  Eh/Bohr at h=1e-3 while the O(h^2) truncation error is ~1e-6 |E'''| --
  h ~ (eps)^(1/3) ~ 1e-3 balances the two.  It also matches the default
  ``[hess] dx`` and the vibrational-intensity displacement already used in
  this code base.
* ``grad_guess`` ('cold' | 'warm', default 'cold'): orbital start-up at each
  displaced geometry.  'cold' re-runs the configured ``[guess] type``
  (huckel/hcore/...) so every displaced energy is bit-for-bit the energy a
  fresh single point at that geometry would produce.  'warm' temporarily sets
  the guess type to 'json', which for the RHF references PT2 requires keeps
  the in-memory MOs from the previous geometry as the SCF start (no file I/O
  involved).  Warm starts inherit the semicanonical orbitals PT2 commits to
  OQP::VEC_MO_A; the SCF re-converges to the same solution through the stored
  density, but 'cold' is the default because it is exactly reproducible
  against independent single points.
* ``grad_gap_warn`` (float, default 1.0e-5 Hartree): root-swap warning floor,
  see below.

Root tracking (honest limitation)
---------------------------------
Multistate roots are energy-ordered.  Near a crossing, displaced geometries
can reorder the adiabatic states, which silently differentiates "state k by
energy order" instead of a diabatically-followed state and contaminates the
FD gradient with the neighbouring surface.  No CI-vector overlap data is
exposed by the PT2 kernels, so this module cannot re-identify states across
displacements; instead it WARNS (in the log) whenever the smallest adjacent
gap at any displaced geometry falls below max(grad_gap_warn,
2 * the largest per-state energy shift caused by the displacement) -- the
regime where an ordering swap is possible.  Swaps with roots outside the
computed set (e.g. single-state PT2 crossing an uncomputed root) are
undetectable.  Penalty-function MECI searches keep a finite gap and are the
recommended way to approach degeneracies with these gradients.
"""

import numpy as np

from oqp.utils.file_utils import dump_log

# Method labels dispatched to numerical PT2 gradients (normalized, lower-case;
# mirrors the PT2 branch of SinglePoint.energy in single_point.py).
PT2_NUMGRAD_METHODS = frozenset({
    'caspt2', 'ms-caspt2', 'mscaspt2', 'xms-caspt2', 'xmscaspt2',
    'mrmp2', 'mcqdpt2', 'xmcqdpt2',
})

DEFAULT_GRAD_STEP_BOHR = 1.0e-3
DEFAULT_GRAD_GUESS = 'cold'
DEFAULT_GAP_WARN = 1.0e-5


def _pt2_option(mol, key, default):
    """Read an optional [pt2] key via dict.get (schema rows may not exist)."""
    return mol.config.get('pt2', {}).get(key, default)


def pt2_numerical_gradient(mol, grad_list, sp=None):
    """Central-difference PT2 gradients for the states in ``grad_list``.

    Parameters
    ----------
    mol : oqp.molecule.Molecule
        Molecule whose current geometry defines the expansion point.  The
        central energies are taken from ``mol.energies`` when present (the
        callers -- compute_grad and the optimizer one_step functions -- run
        ``SinglePoint.energy`` immediately before the gradient); otherwise
        they are computed here.
    grad_list : iterable of int
        State indices (into ``mol.energies``) whose gradients are requested.
        All computed states share the same displaced energies, so the full
        ``(nstate, natom, 3)`` array is returned with every state filled.
    sp : SinglePoint, optional
        Reuse an existing SinglePoint driver (tests/advanced callers).

    Returns
    -------
    numpy.ndarray of shape (nstate, natom, 3), Hartree/Bohr.
    """
    # Lazy import: single_point imports this module inside Gradient.gradient,
    # so a module-level import here would be circular.
    from oqp.library.single_point import SinglePoint

    step = float(_pt2_option(mol, 'grad_step', DEFAULT_GRAD_STEP_BOHR))
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError(f'[pt2] grad_step must be a positive number, got {step}')
    guess_mode = str(_pt2_option(mol, 'grad_guess', DEFAULT_GRAD_GUESS)).strip().lower()
    if guess_mode not in ('cold', 'warm'):
        raise ValueError(f"[pt2] grad_guess must be 'cold' or 'warm', got '{guess_mode}'")
    gap_warn = float(_pt2_option(mol, 'grad_gap_warn', DEFAULT_GAP_WARN))

    if sp is None:
        sp = SinglePoint(mol)

    natom = int(mol.data['natom'])
    ncoord = 3 * natom
    x0 = np.asarray(mol.get_system(), dtype=float).reshape(-1).copy()

    # Central-point energies: reuse the ones the caller just computed.
    central = mol.energies
    if not central:
        central = sp.energy(do_init_scf=False)
    central = np.asarray([float(e) for e in central], dtype=float)
    nstate = central.size

    requested = [int(s) for s in np.atleast_1d(np.asarray(grad_list, dtype=int))]
    for s in requested:
        if s < 0 or s >= nstate:
            raise ValueError(
                f'PT2 gradient state {s} is out of range: this calculation '
                f'computed {nstate} PT2 state(s) (indices 0..{nstate - 1}, '
                'ascending in energy). PT2 state indices index mol.energies '
                'directly -- use [properties] grad=0 / [optimize] istate=0 '
                'for single-state PT2, or request more roots via '
                '[pt2] target_roots / nroot.'
            )

    dump_log(mol, title=(
        'PyOQP: PT2 numerical gradient (central differences, '
        'step=%.3e Bohr, %d displaced energies, %s restarts)'
        % (step, 2 * ncoord, guess_mode)))

    e_plus = np.zeros((ncoord, nstate))
    e_minus = np.zeros((ncoord, nstate))
    swap_flags = []  # (coord index, sign, min adjacent gap, max state shift)

    guess_type_saved = mol.config['guess']['type']
    try:
        if guess_mode == 'warm':
            # For the closed-shell RHF references PT2 enforces, guess type
            # 'json' keeps the in-memory MOs/density untouched as the SCF
            # start (see oqp.library.guess); no file is read.
            mol.config['guess']['type'] = 'json'

        for i in range(ncoord):
            for sign, store in ((1.0, e_plus), (-1.0, e_minus)):
                x = x0.copy()
                x[i] += sign * step
                mol.update_system(x)
                energies = sp.energy(do_init_scf=False)
                if len(energies) != nstate:
                    raise RuntimeError(
                        'PT2 numerical gradient: displaced geometry '
                        f'(coord {i}, {"+" if sign > 0 else "-"}{step:.1e} Bohr) '
                        f'returned {len(energies)} PT2 states, expected {nstate}. '
                        'The reference/active space changed across the '
                        'displacement; the FD gradient is not well defined here.'
                    )
                arr = np.asarray([float(e) for e in energies], dtype=float)
                store[i, :] = arr
                if nstate > 1:
                    min_gap = float(np.min(np.diff(arr)))
                    max_shift = float(np.max(np.abs(arr - central)))
                    if min_gap < max(gap_warn, 2.0 * max_shift):
                        swap_flags.append((i, sign, min_gap, max_shift))
    finally:
        # Restore the user guess and the central geometry, and re-converge the
        # pipeline at the expansion point so mol leaves this function holding
        # central-point data (orbitals, energies, PT2 tags) -- not the data of
        # the last displaced geometry.  Guarded so a failure inside the loop
        # is not masked by a failing restore.
        mol.config['guess']['type'] = guess_type_saved
        try:
            mol.update_system(x0)
            closure = np.asarray(
                [float(e) for e in sp.energy(do_init_scf=False)], dtype=float)
            drift = float(np.max(np.abs(closure - central))) if closure.size == nstate \
                else float('inf')
            if drift > 1.0e-10:
                dump_log(mol, title=(
                    'PyOQP: WARNING PT2 numerical gradient closure check: '
                    're-computed central energies drifted by %.3e Hartree '
                    '(state pipeline is not reproducing itself; the gradient '
                    'may be contaminated)' % drift))
        except Exception:  # pragma: no cover - keep the original error visible
            dump_log(mol, title=(
                'PyOQP: WARNING PT2 numerical gradient could not restore the '
                'central geometry state; molecule data corresponds to a '
                'displaced geometry'))

    if swap_flags:
        worst = min(swap_flags, key=lambda t: t[2])
        dump_log(mol, title=(
            'PyOQP: WARNING PT2 numerical gradient: possible root swapping -- '
            '%d of %d displaced points have adjacent PT2 states closer than '
            'twice the displacement-induced energy shift (worst: gap=%.3e '
            'Hartree at coord %d%s). Energy-ordered states may have reordered '
            'across displacements; treat these gradients with caution near '
            'crossings (see pt2_numgrad module docs).'
            % (len(swap_flags), 2 * ncoord, worst[2], worst[0],
               '+' if worst[1] > 0 else '-')))

    # g[state, coord] = (E+ - E-) / (2h), reshaped to (nstate, natom, 3)
    grads = ((e_plus - e_minus) / (2.0 * step)).T.reshape((nstate, natom, 3)).copy()

    for s in requested:
        dump_log(mol, title='PyOQP: PT2 Numerical Gradient of Root %s' % s)

    return grads
