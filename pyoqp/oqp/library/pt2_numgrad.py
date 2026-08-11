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

Configuration ([pt2] section; read with dict.get so an older config without
these rows still runs on the documented defaults)
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
* ``grad_ranks_per_group`` (int, default 0 = automatic): under MPI, how many
  ranks cooperate on ONE displaced energy.  See "MPI" below.  0 gives every
  displacement its own rank (``ngroup = min(world_size, 2*ncoord)``), which
  is the right choice whenever there are at least as many displacements as
  ranks; raise it when the world is much larger than ``2*ncoord`` or when one
  energy does not fit in one rank's memory.  Ignored without mpi4py, on one
  rank, and when the molecule was built with ``usempi=False``.

MPI
---
The displaced energies are independent, but they are NOT independent *ranks*
of work: OpenQP already parallelises the inside of a single energy evaluation
(the native layer receives a communicator through ``MPIManager.set_mpi_comm``
and decomposes the SCF, the integrals and the DFT grid collectively over it).
Handing rank r a different geometry therefore leaves every rank waiting in a
collective the others never enter; an earlier attempt at exactly that
returned a silently ZERO gradient on 4 ranks.

So the fan-out is one level up: ``MPIManager.task_groups`` splits COMM_WORLD
into per-displacement groups and installs the sub-communicator as the
manager's communicator (and as ``infos%mpiinfo%comm`` inside liboqp) for the
duration of the loop.  A group runs whole displacements by itself and MPI
inside one evaluation still works, over the group.  Reassembly zeroes the
non-root ranks of each group before the world reduction, because every rank
of a group holds the same energies and a plain sum would count each of them
``ranks_per_group`` times.  The log is written by world rank 0 only, so under
a split it reports the displacements of group 0 and not the others.

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
from oqp.utils.mpi_utils import MPIManager

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
    ranks_per_group = int(_pt2_option(mol, 'grad_ranks_per_group', 0))

    if sp is None:
        sp = SinglePoint(mol)

    natom = int(mol.data['natom'])
    ncoord = 3 * natom
    x0 = np.asarray(mol.get_system(), dtype=float).reshape(-1).copy()

    # Central-point energies: reuse the ones the caller just computed.
    central = mol.energies
    if not central:
        # Same init_scf policy as the displaced points below: the central and
        # displaced energies must come from the SAME pipeline or the finite
        # difference is not a derivative of anything.
        central = sp.energy(do_init_scf=(guess_mode != 'warm'))
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

    # The 2*ncoord displaced energies are independent, but they cannot be
    # dealt out to individual ranks of COMM_WORLD: OpenQP runs MPI *inside* a
    # single energy evaluation, so a rank given its own geometry sits in a
    # collective the rest of the world never enters (an earlier attempt at
    # that returned a silently ZERO gradient on 4 ranks).  MPIManager.
    # task_groups splits COMM_WORLD into per-displacement groups and installs
    # the sub-communicator both here and in liboqp, so the SCF decomposes
    # WITHIN a group.  See the MPI section of the module docstring.
    tasks = [(i, sign) for i in range(ncoord) for sign in (1.0, -1.0)]

    mpi = MPIManager()
    fanout = bool(getattr(mol, 'usempi', True)) and bool(mpi.use_mpi)
    # Warm restarts chain each displacement onto the previous one's converged
    # state, which is only well defined in a single sequence.  Split across task
    # groups, each group walks a DIFFERENT subset -- with two groups one takes
    # the + displacements and the other the -, where a serial run alternates
    # them -- so on a system with more than one SCF basin the gradient depends
    # on the rank layout.  The closure is worse: ranks from different groups
    # would enter the collective SCF holding different densities.  Cold mode
    # restores a common snapshot before every point and is unaffected, so only
    # warm gives up the fan-out.
    if fanout and guess_mode == 'warm':
        fanout = False
        dump_log(mol, title=(
            'PyOQP: PT2 numerical gradient: [pt2] grad_guess=warm disables the '
            'MPI displacement fan-out, because chaining restarts across task '
            'groups makes the result depend on the rank layout. Use '
            'grad_guess=cold to fan out across groups.'))

    e_plus = np.zeros((ncoord, nstate))
    e_minus = np.zeros((ncoord, nstate))
    swap_flags = []  # (coord index, sign, min adjacent gap, max state shift)

    guess_type_saved = mol.config['guess']['type']
    # Snapshot the starting orbitals.  `cold` promises every displacement starts
    # from the same guess, but a file-based guess ([guess] type=json) is loaded
    # once at runner setup and guess() is then a no-op for RHF -- so each cold
    # displacement actually started from the PREVIOUS displaced geometry's
    # PT2-mutated state.  That makes +h and -h able to land in different SCF
    # basins and the gradient order-dependent, which is exactly what cold is
    # supposed to rule out.  Restoring this before every displaced evaluation
    # makes the policy true for any guess type (it is a no-op for hcore/huckel,
    # which rebuild from scratch anyway).
    _restart_tags = ('OQP::VEC_MO_A', 'OQP::VEC_MO_B', 'OQP::DM_A', 'OQP::DM_B',
                     'OQP::E_MO_A', 'OQP::E_MO_B')
    _restart_start = {}
    for _tag in _restart_tags:
        try:
            _restart_start[_tag] = np.array(mol.data[_tag], dtype=float, copy=True)
        except Exception:
            pass
    # Whether the displaced-energy loop below is already unwinding.  The
    # restore in `finally` must not mask that failure -- but when nothing is
    # in flight, a failed restore is itself fatal and has to propagate, or
    # this function would hand back a gradient while `mol` still holds a
    # displaced geometry.
    loop_failed = False
    try:
        if guess_mode == 'warm':
            # For the closed-shell RHF references PT2 enforces, guess type
            # 'json' keeps the in-memory MOs/density untouched as the SCF
            # start (see oqp.library.guess); no file is read.
            mol.config['guess']['type'] = 'json'

        with mpi.task_groups(len(tasks), ranks_per_group=ranks_per_group,
                             enabled=fanout,
                             on_switch=lambda: mpi.set_mpi_comm(mol.data)) as groups:
            dump_log(mol, title=(
                'PyOQP: PT2 numerical gradient (central differences, '
                'step=%.3e Bohr, %d displaced energies, %s restarts, '
                '%d MPI rank group(s))'
                % (step, 2 * ncoord, guess_mode, groups.ngroup)))

            for task in groups.indices:
                i, sign = tasks[task]
                store = e_plus if sign > 0 else e_minus
                x = x0.copy()
                x[i] += sign * step
                mol.update_system(x)
                if guess_mode != 'warm':
                    # Restoring the orbitals alone was not enough: the RHF json
                    # guess is a no-op, so the SCF actually starts from the
                    # stored DENSITY, which still carried the previous
                    # displacement.  Restore every restart tag the SCF reads.
                    for _tag, _val in _restart_start.items():
                        mol.data[_tag][...] = _val
                # The central point runs SinglePoint.energy() with the
                # configured [scf] init_scf preconvergence stage; forcing it off
                # here differentiated a different pipeline from the one the
                # central energy came from, and removed the convergence aid
                # exactly where a displaced geometry is most likely to need it.
                # Skip it only under an explicit warm-start policy.
                energies = sp.energy(do_init_scf=(guess_mode != "warm"))
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

        # Record what the split actually did.  Tests and callers need to be
        # able to assert that the work really was spread (a run that quietly
        # evaluated every displacement on every rank produces the right
        # numbers and no speed-up at all).
        mol.pt2_numgrad_layout = {
            'ntask': len(tasks),
            'ngroup': groups.ngroup,
            'group': groups.group,
            'ntask_local': len(list(groups.indices)),
            'is_group_root': groups.is_group_root,
        }

        # Outside the split: COMM_WORLD is back, and reduce_sum drops the
        # duplicate copies the non-root ranks of each group are holding.
        e_plus = groups.reduce_sum(e_plus)
        e_minus = groups.reduce_sum(e_minus)
        swap_flags = groups.gather_list(swap_flags)
    except BaseException:
        loop_failed = True
        raise
    finally:
        # Restore the user guess and the central geometry, and re-converge the
        # pipeline at the expansion point so mol leaves this function holding
        # central-point data (orbitals, energies, PT2 tags) -- not the data of
        # the last displaced geometry.  Guarded so a failure inside the loop
        # is not masked by a failing restore.
        mol.config['guess']['type'] = guess_type_saved
        try:
            mol.update_system(x0)
            if guess_mode != 'warm':
                for _tag, _val in _restart_start.items():
                    mol.data[_tag][...] = _val
            # Restore through the configured pipeline too.  Forcing init_scf
            # off here left `mol` holding orbitals and energies from a
            # different pipeline than every other evaluation -- and if the
            # auxiliary stage selects or stabilizes a different SCF solution,
            # the drift check below only warns while that state persists.
            closure = np.asarray(
                [float(e) for e in sp.energy(do_init_scf=(guess_mode != 'warm'))],
                dtype=float)
            drift = float(np.max(np.abs(closure - central))) if closure.size == nstate \
                else float('inf')
            if drift > 1.0e-10:
                dump_log(mol, title=(
                    'PyOQP: WARNING PT2 numerical gradient closure check: '
                    're-computed central energies drifted by %.3e Hartree '
                    '(state pipeline is not reproducing itself; the gradient '
                    'may be contaminated)' % drift))
        except Exception:
            dump_log(mol, title=(
                'PyOQP: WARNING PT2 numerical gradient could not restore the '
                'central geometry state; molecule data corresponds to a '
                'displaced geometry'))
            # Only swallow this while another exception is already unwinding,
            # so the original cause stays visible.  Otherwise the displaced
            # energies were all fine and the restore is the only thing that
            # failed -- returning a gradient here would hand a
            # gradient-driven optimizer a molecule whose central electronic
            # state was never restored.
            if not loop_failed:
                raise

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
