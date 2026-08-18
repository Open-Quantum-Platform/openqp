"""Central-difference nuclear gradients for multireference wavefunctions.

OpenQP uses the analytic kernel for state-specific CASSCF.  This module
differentiates the converged total energies returned by ``SinglePoint.energy``
over Cartesian displacements for SA-CASSCF and the CASPT2/NEVPT2/QDPT2
families.  It can also be called explicitly for a state-specific CASSCF
finite-difference reference, but the production gradient dispatch does not use
that as a substitute for the analytic derivative.  The numerical driver is
used by ``grad``, ``optimize``, ``ts``, ``mep``, and ``irc``.

SA-CASSCF publishes the roots selected by ``[state_average]``; their indices in
``[properties] grad`` and ``[optimize] istate`` address that published list.
PT2 indices similarly address the energies selected by ``[pt2]
target_roots``/``nroot``.

The controls ``grad_step``, ``grad_guess``, ``grad_gap_warn``, and
``grad_ranks_per_group`` are read from ``[casscf]`` for CASSCF/SA-CASSCF and
from ``[pt2]`` for PT2.  ``grad_guess=cold`` restores the same central-point
SCF data before each displaced calculation; ``warm`` follows the preceding
solution serially.  Coordinates are Bohr, energies Hartree, and gradients
Hartree/Bohr with shape ``(nstate, natom, 3)``.

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

Root tracking limitation
------------------------
The displaced calculations do not yet match CI vectors between geometries.
When two published roots become closer than the displacement-induced energy
change, the log warns that their energy order may have exchanged.  Crossings
with an uncomputed root cannot be detected.
"""

import numpy as np

from oqp.utils.file_utils import dump_log
from oqp.utils.mpi_utils import MPIManager

# Method labels dispatched to central-difference nuclear gradients.
CASSCF_NUMGRAD_METHODS = frozenset({
    'casscf', 'sa-casscf', 'sacasscf',
})
PT2_NUMGRAD_METHODS = frozenset({
    'caspt2', 'ms-caspt2', 'mscaspt2', 'xms-caspt2', 'xmscaspt2',
    'mrmp2', 'mcqdpt2', 'xmcqdpt2',
})
WF_NUMGRAD_METHODS = CASSCF_NUMGRAD_METHODS | PT2_NUMGRAD_METHODS

DEFAULT_GRAD_STEP_BOHR = 1.0e-3
DEFAULT_GRAD_GUESS = 'cold'
DEFAULT_GAP_WARN = 1.0e-5


def _method_label(mol):
    return str(mol.config['input']['method']).strip().lower()


def _option_section(mol):
    return 'casscf' if _method_label(mol) in CASSCF_NUMGRAD_METHODS else 'pt2'


def _gradient_option(mol, key, default):
    """Read a method-family numerical-gradient option."""
    return mol.config.get(_option_section(mol), {}).get(key, default)


def _family_label(mol):
    return 'CASSCF' if _method_label(mol) in CASSCF_NUMGRAD_METHODS else 'PT2'


def wavefunction_numerical_gradient(mol, grad_list, sp=None):
    """Central-difference nuclear gradients for the states in ``grad_list``.

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

    method = _method_label(mol)
    family = _family_label(mol)
    section = _option_section(mol)
    if method not in WF_NUMGRAD_METHODS:
        raise ValueError(f'Numerical wavefunction gradient is unavailable for {method}')

    step = float(_gradient_option(mol, 'grad_step', DEFAULT_GRAD_STEP_BOHR))
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError(f'[{section}] grad_step must be a positive number, got {step}')
    guess_mode = str(_gradient_option(mol, 'grad_guess', DEFAULT_GRAD_GUESS)).strip().lower()
    if guess_mode not in ('cold', 'warm'):
        raise ValueError(
            f"[{section}] grad_guess must be 'cold' or 'warm', got '{guess_mode}'")
    gap_warn = float(_gradient_option(mol, 'grad_gap_warn', DEFAULT_GAP_WARN))
    ranks_per_group = int(_gradient_option(mol, 'grad_ranks_per_group', 0))

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
                f'{family} gradient state {s} is out of range: this '
                f'calculation published {nstate} state(s), with indices '
                f'0..{nstate - 1}. State indices address mol.energies directly.'
            )

    # The 2*ncoord displaced energies are independent, but they cannot be
    # dealt out to individual ranks of COMM_WORLD: OpenQP runs MPI *inside* a
    # single energy evaluation, so a rank given its own geometry sits in a
    # collective the rest of the world never enters (an earlier attempt at
    # that returned a silently ZERO gradient on 4 ranks).
    # MPIManager.task_groups splits COMM_WORLD into per-displacement groups
    # and installs the sub-communicator both here and in liboqp, so the SCF
    # decomposes WITHIN a group.  See the MPI section of the module docstring.
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
            f'PyOQP: {family} numerical gradient: [{section}] '
            'grad_guess=warm disables the '
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
                f'PyOQP: {family} numerical gradient (central differences, '
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
                        f'{family} numerical gradient: displaced geometry '
                        f'(coord {i}, {"+" if sign > 0 else "-"}{step:.1e} Bohr) '
                        f'returned {len(energies)} states, expected {nstate}. '
                        'The reference/active space changed across the '
                        'displacement; the FD gradient is not well defined here.'
                    )
                arr = np.asarray([float(e) for e in energies], dtype=float)
                store[i, :] = arr
                if nstate > 1:
                    min_gap = float(np.min(np.diff(np.sort(arr))))
                    max_shift = float(np.max(np.abs(arr - central)))
                    if min_gap < max(gap_warn, 2.0 * max_shift):
                        swap_flags.append((i, sign, min_gap, max_shift))

        # Record what the split actually did.  Tests and callers need to be
        # able to assert that the work really was spread (a run that quietly
        # evaluated every displacement on every rank produces the right
        # numbers and no speed-up at all).
        layout = {
            'ntask': len(tasks),
            'ngroup': groups.ngroup,
            'group': groups.group,
            'ntask_local': len(list(groups.indices)),
            'is_group_root': groups.is_group_root,
        }
        mol.wf_numgrad_layout = layout
        if method in PT2_NUMGRAD_METHODS:
            # Backward-compatible diagnostic attribute used by existing tests
            # and external scripts written for the original PT2-only driver.
            mol.pt2_numgrad_layout = layout

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
                    f'PyOQP: WARNING {family} numerical gradient closure check: '
                    're-computed central energies drifted by %.3e Hartree '
                    '(state pipeline is not reproducing itself; the gradient '
                    'may be contaminated)' % drift))
        except Exception:
            dump_log(mol, title=(
                f'PyOQP: WARNING {family} numerical gradient could not restore the '
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
            f'PyOQP: WARNING {family} numerical gradient: possible root swapping -- '
            '%d of %d displaced points have adjacent states closer than '
            'twice the displacement-induced energy shift (worst: gap=%.3e '
            'Hartree at coord %d%s). Energy-ordered states may have reordered '
            'across displacements; treat these gradients with caution near '
            'crossings (see wf_numgrad module docs).'
            % (len(swap_flags), 2 * ncoord, worst[2], worst[0],
               '+' if worst[1] > 0 else '-')))

    # g[state, coord] = (E+ - E-) / (2h), reshaped to (nstate, natom, 3)
    grads = ((e_plus - e_minus) / (2.0 * step)).T.reshape((nstate, natom, 3)).copy()

    for s in requested:
        dump_log(mol, title=f'PyOQP: {family} Numerical Gradient of Root {s}')

    return grads
