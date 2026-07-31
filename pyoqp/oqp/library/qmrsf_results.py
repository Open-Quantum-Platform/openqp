"""QMRSF-icPT2 results post-processing.

The live Fortran routine ``tdhf_qmrsf_icpt2`` writes a fixed-name validation
dump ``qmrsf_icpt2_full_live.dat`` into the run's working directory.  This module
turns that raw dump into clean, consumable results: a parser, a results builder,
a JSON writer, and a formatted log table.

Dump file format (free-format ``es24.16`` text, written by
``source/modules/tdhf_qmrsf_icpt2.F90``)::

    line 1            : "<norb_w> <nPd>"        two ints (nPd = number of dressed roots, 36)
    next norb_w lines : h_win rows               (norb_w values each)
    next norb_w**3    : eri_win blocks (p,q,r,:) (norb_w values each)
    1 line            : ecore                    (nuc. repulsion + frozen-core const, Hartree)
    1 line            : eps_win                  (norb_w ROHF MO energies)
    1 line            : eP                        (nPd bare CAS-root ELECTRONIC energies)
    1 line            : edr_en                    (nPd icPT2 Epstein-Nesbet ELECTRONIC energies)
    1 line            : edr_dy                    (nPd icPT2 Dyall ELECTRONIC energies)

All ``e*`` arrays are electronic energies; add ``ecore`` for totals.
"""

import json

# CODATA Hartree -> eV (matches the value used elsewhere in the QMRSF stack).
HARTREE_TO_EV = 27.211386245988


def _read_floats(line):
    """Parse a whitespace-separated line of Fortran ``es24.16`` floats."""
    return [float(tok) for tok in line.split()]


def parse_qmrsf_icpt2_dump(path):
    """Read the scalar results from a ``qmrsf_icpt2_full_live.dat`` dump.

    Only the scalar spectra are returned (norb_w, nPd, ecore, eps_win, eP,
    edr_en, edr_dy).  The large ``h_win`` / ``eri_win`` integral blocks are
    skipped, but the exact number of their lines is consumed so the file
    pointer reaches the trailing scalar records.

    Parameters
    ----------
    path : str
        Filesystem path to the dump file.

    Returns
    -------
    dict
        ``{"norb_w", "nPd", "ecore", "eps_win", "eP", "edr_en", "edr_dy"}``.

    Raises
    ------
    ValueError
        If the file is truncated or the header is malformed.
    """
    with open(path, 'r') as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError('QMRSF-icPT2 dump is empty: %s' % path)

    header = lines[0].split()
    if len(header) < 2:
        raise ValueError('QMRSF-icPT2 dump header malformed (expected "<norb_w> <nPd>"): %r' % lines[0])
    norb_w = int(header[0])
    nPd = int(header[1])
    if norb_w <= 0 or nPd <= 0:
        raise ValueError('QMRSF-icPT2 dump header has non-positive sizes: norb_w=%d nPd=%d' % (norb_w, nPd))

    # Skip the integral blocks: norb_w rows of h_win + norb_w**3 rows of eri_win.
    skip = norb_w + norb_w ** 3
    idx = 1 + skip  # first trailing scalar record (ecore)

    # 5 trailing scalar records: ecore, eps_win, eP, edr_en, edr_dy.
    needed = idx + 5
    if len(lines) < needed:
        raise ValueError(
            'QMRSF-icPT2 dump is truncated: have %d lines, need at least %d '
            '(norb_w=%d, nPd=%d)' % (len(lines), needed, norb_w, nPd)
        )

    ecore = _read_floats(lines[idx])[0]
    eps_win = _read_floats(lines[idx + 1])
    eP = _read_floats(lines[idx + 2])
    edr_en = _read_floats(lines[idx + 3])
    edr_dy = _read_floats(lines[idx + 4])

    for name, arr, n in (
        ('eps_win', eps_win, norb_w),
        ('eP', eP, nPd),
        ('edr_en', edr_en, nPd),
        ('edr_dy', edr_dy, nPd),
    ):
        if len(arr) != n:
            raise ValueError(
                'QMRSF-icPT2 dump record %s has %d values, expected %d' % (name, len(arr), n)
            )

    # Optional <S^2> records (s2_cas, s2_en, s2_dy) -- present for engine-path runs;
    # absent for the CASCI-fallback path. Backward compatible.
    def _opt(j):
        return _read_floats(lines[idx + j]) if len(lines) > idx + j else None
    s2_cas, s2_en, s2_dy = _opt(5), _opt(6), _opt(7)

    return {
        'norb_w': norb_w,
        'nPd': nPd,
        'ecore': ecore,
        'eps_win': eps_win,
        'eP': eP,
        'edr_en': edr_en,
        'edr_dy': edr_dy,
        's2_cas': s2_cas,
        's2_en': s2_en,
        's2_dy': s2_dy,
    }


def build_qmrsf_results(dump, ref_energy):
    """Assemble a clean results dict from a parsed dump and the reference energy.

    Total energies are ``E_electronic + ecore``.  Excitation energies are
    measured relative to each method's own ground state (state 0) and converted
    Hartree -> eV.

    Parameters
    ----------
    dump : dict
        Output of :func:`parse_qmrsf_icpt2_dump`.
    ref_energy : float
        The converged quintet ROHF reference energy (Hartree).

    Returns
    -------
    dict
        Results dict with per-state totals (Hartree) and excitation energies (eV).
    """
    ecore = dump['ecore']
    nPd = dump['nPd']
    eP = dump['eP']
    edr_en = dump['edr_en']
    edr_dy = dump['edr_dy']

    # Each method's ground-state total (state 0) anchors its excitation energies.
    cas0 = eP[0] + ecore
    en0 = edr_en[0] + ecore
    dy0 = edr_dy[0] + ecore

    s2c, s2e, s2d = dump.get('s2_cas'), dump.get('s2_en'), dump.get('s2_dy')

    def _mult(arr, i):
        if not arr:
            return None
        return int(round((1.0 + 4.0 * max(arr[i], 0.0)) ** 0.5))   # 2S+1 from <S^2>=S(S+1)

    states = []
    for i in range(nPd):
        e_cas = eP[i] + ecore
        e_en = edr_en[i] + ecore
        e_dy = edr_dy[i] + ecore
        states.append({
            'index': i,
            'E_CAS': e_cas,
            'E_icPT2_EN': e_en,
            'E_icPT2_Dyall': e_dy,
            'exc_CAS_eV': (e_cas - cas0) * HARTREE_TO_EV,
            'exc_EN_eV': (e_en - en0) * HARTREE_TO_EV,
            'exc_Dyall_eV': (e_dy - dy0) * HARTREE_TO_EV,
            'mult_cas': _mult(s2c, i),
            'mult_en': _mult(s2e, i),
            'mult_dy': _mult(s2d, i),
        })

    return {
        'method': 'QMRSF-icPT2',
        'reference_energy': float(ref_energy),
        'ecore': ecore,
        'n_window_orbitals': dump['norb_w'],
        'n_states': nPd,
        'states': states,
    }


def write_qmrsf_json(results, json_path):
    """Write the results dict to ``json_path`` as pretty-printed JSON."""
    with open(json_path, 'w') as fh:
        json.dump(results, fh, indent=2)
    return json_path


def parse_qmrsf_dk_dump(path):
    """Read the scalar results from a ``qmrsf_dk_full_live.dat`` dump.

    Dump file format (free-format ``es24.16`` text, written by
    ``source/modules/tdhf_qmrsf_dk.F90``)::

        line 1            : "<nact> <ndet> <nopen> <nclosed>"  four ints
        next nact lines   : h_act rows                          (nact values each)
        next nact**3      : eri_act blocks (p,q,r,:)            (nact values each)
        1 line            : ecore                               (Hartree)
        1 line            : omega_d                             (nclosed bare 0OS doubles)
        1 line            : cas_ref                             (ndet CAS electronic energies)
        1 line            : dressed                             (ndet DK electronic energies)
        1 line            : adiab                               (nopen adiabatic electronic energies)
        1 line            : g1_cas g1_exact gap_adiab gap_dressed   (4 gate metrics)
        1 line            : nmiss                               (int; #0OS doubles adiabatic misses)

    All ``e*`` arrays are electronic energies; add ``ecore`` for totals.

    Returns
    -------
    dict
        ``{"nact","ndet","nopen","nclosed","ecore","omega_d","cas_ref",
           "dressed","adiab","gate1_dk_cas","gate1_dk_exact","gap_adiab",
           "gap_dressed","nmiss"}``.
    """
    with open(path, 'r') as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError('QMRSF-DK dump is empty: %s' % path)

    if lines[0].split()[0] == 'QMRSF_DK_NATIVE_V1':
        return _parse_qmrsf_dk_native(lines, path)

    header = lines[0].split()
    if len(header) < 4:
        raise ValueError(
            'QMRSF-DK dump header malformed (expected "<nact> <ndet> <nopen> <nclosed>"): %r'
            % lines[0])
    nact, ndet, nopen, nclosed = (int(header[0]), int(header[1]),
                                  int(header[2]), int(header[3]))
    if min(nact, ndet, nopen, nclosed) <= 0:
        raise ValueError('QMRSF-DK dump header has non-positive sizes: %r' % lines[0])

    # Skip the integral blocks: nact rows of h_act + nact**3 rows of eri_act.
    idx = 1 + nact + nact ** 3

    needed = idx + 6  # ecore, omega_d, cas_ref, dressed, adiab, gates, nmiss = 7 records
    if len(lines) < needed + 1:
        raise ValueError(
            'QMRSF-DK dump is truncated: have %d lines, need at least %d '
            '(nact=%d, ndet=%d)' % (len(lines), needed + 1, nact, ndet))

    ecore = _read_floats(lines[idx])[0]
    omega_d = _read_floats(lines[idx + 1])
    cas_ref = _read_floats(lines[idx + 2])
    dressed = _read_floats(lines[idx + 3])
    adiab = _read_floats(lines[idx + 4])
    gates = _read_floats(lines[idx + 5])
    nmiss = int(lines[idx + 6].split()[0])

    for name, arr, n in (
        ('omega_d', omega_d, nclosed),
        ('cas_ref', cas_ref, ndet),
        ('dressed', dressed, ndet),
        ('adiab', adiab, nopen),
    ):
        if len(arr) != n:
            raise ValueError(
                'QMRSF-DK dump record %s has %d values, expected %d' % (name, len(arr), n))
    if len(gates) < 4:
        raise ValueError('QMRSF-DK dump gate record has %d values, expected 4' % len(gates))

    # --- optional DFT-dressed extension (records D0..D3 appended after nmiss) ---
    #   D0: "<is_dft 0|1> <kscale>" ; D1: cas_dft (ndet) ; D2: dk_dft (ndet) ;
    #   D3: adiab_dft (nopen).  Older dumps lack these -> is_dft=False.
    is_dft, kscale = False, 1.0
    cas_dft, dk_dft, adiab_dft = None, None, None
    d0 = idx + 7
    if len(lines) >= d0 + 4:
        head = lines[d0].split()
        try:
            is_dft = (int(head[0]) == 1)
            kscale = float(head[1])
            cas_dft = _read_floats(lines[d0 + 1])
            dk_dft = _read_floats(lines[d0 + 2])
            adiab_dft = _read_floats(lines[d0 + 3])
            if len(cas_dft) != ndet or len(dk_dft) != ndet or len(adiab_dft) != nopen:
                is_dft, cas_dft, dk_dft, adiab_dft = False, None, None, None
        except (ValueError, IndexError):
            is_dft, cas_dft, dk_dft, adiab_dft = False, None, None, None

    # Optional trailing <S^2> record (bare-CAS spin label), appended after the
    # DFT records (or absent on older dumps). Backward compatible.
    s2_dk = _read_floats(lines[d0 + 4]) if len(lines) > d0 + 4 else None
    if s2_dk is not None and len(s2_dk) != ndet:
        s2_dk = None

    # Optional D5/D6: CSF spin-adapted DRESSED spectrum + exact multiplicity
    # (multiplicity-block projection of the dressed H -- spin-pure by construction).
    sa_eval = _read_floats(lines[d0 + 5]) if len(lines) > d0 + 5 else None
    sa_s2 = _read_floats(lines[d0 + 6]) if len(lines) > d0 + 6 else None
    if sa_eval is not None and len(sa_eval) != ndet:
        sa_eval = None
    if sa_s2 is not None and len(sa_s2) != ndet:
        sa_s2 = None

    return {
        'nact': nact,
        'ndet': ndet,
        'nopen': nopen,
        'nclosed': nclosed,
        'ecore': ecore,
        'omega_d': omega_d,
        'cas_ref': cas_ref,
        'dressed': dressed,
        'adiab': adiab,
        'gate1_dk_cas': gates[0],
        'gate1_dk_exact': gates[1],
        'gap_adiab': gates[2],
        'gap_dressed': gates[3],
        'nmiss': nmiss,
        'is_dft': is_dft,
        'kscale': kscale,
        'cas_dft': cas_dft,
        'dk_dft': dk_dft,
        'adiab_dft': adiab_dft,
        's2_dk': s2_dk,
        'sa_eval': sa_eval,
        'sa_s2': sa_s2,
    }


def _parse_qmrsf_dk_native(lines, path):
    """Parse the compact dump written by the OpenQP-native paper builder."""
    header = lines[0].split()
    if len(header) != 5:
        raise ValueError('Malformed native QMRSF-DK header in %s: %r' % (path, lines[0]))

    try:
        nact, nsinglet, ntriplet, nquintet = map(int, header[1:])
    except ValueError as exc:
        raise ValueError('Non-integer native QMRSF-DK dimensions in %s' % path) from exc

    if (nact, nsinglet, ntriplet, nquintet) != (4, 20, 15, 1):
        raise ValueError(
            'Unsupported native QMRSF-DK dimensions in %s: %s'
            % (path, (nact, nsinglet, ntriplet, nquintet)))
    if len(lines) < 7:
        raise ValueError('Native QMRSF-DK dump is truncated: %s' % path)

    scales = _read_floats(lines[1])
    singlet = _read_floats(lines[3])
    triplet = _read_floats(lines[4])
    quintet = _read_floats(lines[5])
    diagnostics = _read_floats(lines[6])
    if len(scales) != 2 or len(singlet) != nsinglet or len(triplet) != ntriplet:
        raise ValueError('Native QMRSF-DK dump has inconsistent record lengths: %s' % path)
    if len(quintet) != nquintet or len(diagnostics) != 2:
        raise ValueError('Native QMRSF-DK dump has inconsistent diagnostics: %s' % path)

    a_quintet = _read_floats(lines[2])[0]
    if abs(quintet[0] - a_quintet) > 1.0e-10:
        raise ValueError(
            'Native QMRSF-DK quintet anchor is inconsistent in %s (%.16g vs %.16g)'
            % (path, a_quintet, quintet[0]))

    return {
        'format': 'native_v1',
        'nact': nact,
        'nsinglet': nsinglet,
        'ntriplet': ntriplet,
        'nquintet': nquintet,
        'c_h': scales[0],
        'c_ref': scales[1],
        'a_quintet': a_quintet,
        'singlet': singlet,
        'triplet': triplet,
        'quintet': quintet,
        'orthonormal_error': diagnostics[0],
        'discarded_cross_spin_block_max': diagnostics[1],
    }


def build_qmrsf_dk_results(dump, ref_energy):
    """Assemble a clean DK results dict from a parsed dump and the reference energy.

    Total energies are ``E_electronic + ecore``.  Excitation energies are
    measured relative to the DK ground state (state 0) and converted Hartree->eV.
    The validation gates are carried through verbatim.

    Parameters
    ----------
    dump : dict
        Output of :func:`parse_qmrsf_dk_dump`.
    ref_energy : float
        The converged quintet ROHF reference energy (Hartree).

    Returns
    -------
    dict
        Results dict with per-state DK/CAS totals (Hartree), excitation energies
        (eV), the 0OS classification, and the validation-gate summary.
    """
    if dump.get('format') == 'native_v1':
        return _build_qmrsf_dk_native_results(dump, ref_energy)

    ecore = dump['ecore']
    ndet = dump['ndet']
    cas_ref = dump['cas_ref']
    dressed = dump['dressed']

    dk0 = dressed[0] + ecore
    cas0 = cas_ref[0] + ecore

    is_dft = bool(dump.get('is_dft'))
    cas_dft = dump.get('cas_dft')
    dk_dft = dump.get('dk_dft')
    dk_dft0 = (dk_dft[0] + ecore) if is_dft and dk_dft else None
    cas_dft0 = (cas_dft[0] + ecore) if is_dft and cas_dft else None

    s2 = dump.get('s2_dk')

    # CSF spin-adapted DRESSED spectrum (spin-pure by construction): multiplicity
    # is EXACT (block label), not the nominal bare assignment used for the
    # full-determinant diagonalization.  This is the production singlet spectrum.
    sa_eval = dump.get('sa_eval')
    sa_s2 = dump.get('sa_s2')
    sa0 = (sa_eval[0] + ecore) if (is_dft and sa_eval) else None

    def _mult(arr, i):
        if not arr:
            return None
        return int(round((1.0 + 4.0 * max(arr[i], 0.0)) ** 0.5))   # 2S+1 from <S^2>=S(S+1)

    states = []
    for i in range(ndet):
        e_dk = dressed[i] + ecore
        e_cas = cas_ref[i] + ecore
        st = {
            'index': i,
            'E_DK': e_dk,
            'E_CAS': e_cas,
            'exc_DK_eV': (e_dk - dk0) * HARTREE_TO_EV,
            'exc_CAS_eV': (e_cas - cas0) * HARTREE_TO_EV,
            's2': (s2[i] if s2 else None),
            'mult': _mult(s2, i),       # bare-CAS spin label (nominal for DFT-dressed)
        }
        if is_dft and dk_dft and cas_dft:
            e_dk_dft = dk_dft[i] + ecore
            e_cas_dft = cas_dft[i] + ecore
            st['E_DK_DFT'] = e_dk_dft
            st['E_CAS_DFT'] = e_cas_dft
            st['exc_DK_DFT_eV'] = (e_dk_dft - dk_dft0) * HARTREE_TO_EV
            st['exc_CAS_DFT_eV'] = (e_cas_dft - cas_dft0) * HARTREE_TO_EV
        if is_dft and sa_eval and sa_s2:
            e_sa = sa_eval[i] + ecore
            st['E_DK_DFT_SA'] = e_sa                       # spin-adapted dressed total
            st['exc_DK_DFT_SA_eV'] = (e_sa - sa0) * HARTREE_TO_EV
            st['s2_sa'] = sa_s2[i]                         # exact <S^2> (block label)
            st['mult_sa'] = _mult(sa_s2, i)               # exact multiplicity
        states.append(st)

    gate1 = (dump['gate1_dk_cas'] < 1e-9 and dump['gate1_dk_exact'] < 1e-9)
    gate2 = (dump['nmiss'] == dump['nclosed']
             and dump['gap_adiab'] > 1e-3 and dump['gap_dressed'] < 1e-9)

    return {
        'method': 'QMRSF-DK',
        'reference_energy': float(ref_energy),
        'ecore': ecore,
        'is_dft_dressed': is_dft,
        'kscale': float(dump.get('kscale', 1.0)),
        'n_open_singles': dump['nopen'],
        'n_0os_doubles': dump['nclosed'],
        'n_states': ndet,
        'omega_d_0os': dump['omega_d'],
        'states': states,
        'gates': {
            'gate1_DK_vs_CAS_max_abs': dump['gate1_dk_cas'],
            'gate1_DK_vs_augmented_exact_max_abs': dump['gate1_dk_exact'],
            'gate1_pass': bool(gate1),
            'gate2_n_0os_doubles_missed_by_adiabatic': dump['nmiss'],
            'gate2_worst_adiabatic_gap': dump['gap_adiab'],
            'gate2_worst_dressed_gap': dump['gap_dressed'],
            'gate2_pass': bool(gate2),
        },
    }


def _build_qmrsf_dk_native_results(dump, ref_energy):
    """Apply the paper energy law to the native spin-adapted eigenvalues."""
    anchor = dump['a_quintet']

    def _states(eigenvalues, multiplicity):
        total = [float(ref_energy) + (value - anchor) for value in eigenvalues]
        base = total[0]
        return [
            {
                'index': i,
                'multiplicity': multiplicity,
                'mult': multiplicity,
                'matrix_eigenvalue': value,
                'E_DK': energy,
                'exc_DK_eV': (energy - base) * HARTREE_TO_EV,
            }
            for i, (value, energy) in enumerate(zip(eigenvalues, total))
        ]

    singlets = _states(dump['singlet'], 1)
    triplets = _states(dump['triplet'], 3)
    quintets = _states(dump['quintet'], 5)
    orth_pass = dump['orthonormal_error'] < 1.0e-12

    return {
        'method': 'QMRSF-DK',
        'implementation': 'OpenQP-native paper matrix',
        'reference_energy': float(ref_energy),
        'is_dft_dressed': True,
        'c_H': dump['c_h'],
        'c_ref': dump['c_ref'],
        'a_quintet': anchor,
        'n_states': dump['nsinglet'],
        'n_singlets': dump['nsinglet'],
        'n_triplets': dump['ntriplet'],
        'n_quintets': dump['nquintet'],
        # PyOQP's public excited-state result remains the target singlet block.
        'states': singlets,
        'triplet_states': triplets,
        'quintet_states': quintets,
        'diagnostics': {
            'csf_orthonormal_error': dump['orthonormal_error'],
            'discarded_cross_spin_block_max': dump['discarded_cross_spin_block_max'],
            'csf_orthonormal_pass': orth_pass,
            # Exchange scaling and the determinant-resolved v_xc lift are
            # projected into separate spin-adapted blocks by construction;
            # their discarded cross-block norm is informative, not a failure.
            'pass': bool(orth_pass),
        },
    }


def format_qmrsf_dk_log_table(results, max_states=10):
    """Render an aligned text table of the lowest DK states for the log."""
    if results.get('implementation') == 'OpenQP-native paper matrix':
        return _format_qmrsf_dk_native_log_table(results, max_states)

    states = results.get('states', [])
    n_show = min(max_states, len(states))
    g = results.get('gates', {})

    header_lines = [
        'QMRSF-DK results',
        'reference (quintet ROHF) = %18.10f Hartree' % results['reference_energy'],
        'E_core (nuc + frozen)    = %18.10f Hartree' % results['ecore'],
        'open-shell singles = %d   0OS doubles = %d   states = %d   (showing lowest %d)'
        % (results['n_open_singles'], results['n_0os_doubles'],
           results['n_states'], n_show),
        '',
    ]

    dft_on = results.get('is_dft_dressed') and \
        all('E_DK_DFT' in st for st in states[:n_show])
    if dft_on:
        header_lines.insert(
            4, 'DFT-dressed: KS reference, active exchange scaled by kscale = %.4f'
            % results.get('kscale', 1.0))
        col = '%5s %5s %18s %14s %18s %14s'
        rule = '-' * 80
        table = [
            col % ('state', '2S+1', 'E_DK-DFT(Eh)', 'exc-DFT(eV)', 'E_DK-bare(Eh)', 'exc-bare(eV)'),
            rule,
        ]
        rowfmt = '%5d %5s %18.10f %14.4f %18.10f %14.4f'
        for st in states[:n_show]:
            mult = st.get('mult')
            table.append(rowfmt % (st['index'], ('-' if mult is None else str(mult)),
                                   st['E_DK_DFT'], st['exc_DK_DFT_eV'],
                                   st['E_DK'], st['exc_DK_eV']))
    else:
        col = '%5s %5s %18s %18s %14s'
        rule = '-' * 66
        table = [
            col % ('state', '2S+1', 'E_DK(Eh)', 'E_CAS(Eh)', 'exc-DK(eV)'),
            rule,
        ]
        rowfmt = '%5d %5s %18.10f %18.10f %14.4f'
        for st in states[:n_show]:
            mult = st.get('mult')
            table.append(rowfmt % (st['index'], ('-' if mult is None else str(mult)),
                                   st['E_DK'], st['E_CAS'], st['exc_DK_eV']))

    gate_lines = [
        '',
        'GATE 1 (DK == CAS): %s  (max|E_DK-E_CAS| = %.2e)'
        % ('PASS' if g.get('gate1_pass') else 'FAIL', g.get('gate1_DK_vs_CAS_max_abs', float('nan'))),
        'GATE 2 (adiabatic misses %d 0OS doubles): %s  (worst adiab gap = %.2e)'
        % (g.get('gate2_n_0os_doubles_missed_by_adiabatic', 0),
           'PASS' if g.get('gate2_pass') else 'FAIL',
           g.get('gate2_worst_adiabatic_gap', float('nan'))),
    ]

    return '\n'.join(header_lines + table + gate_lines)


def _format_qmrsf_dk_native_log_table(results, max_states=10):
    states = results.get('states', [])
    n_show = min(max_states, len(states))
    diag = results.get('diagnostics', {})
    header = [
        'QMRSF-DK results (OpenQP-native paper matrix)',
        'reference quintet energy = %18.10f Hartree' % results['reference_energy'],
        'active exchange c_H = %.6f   reference exchange c_ref = %.6f'
        % (results['c_H'], results['c_ref']),
        'spin-adapted roots: %d singlets, %d triplets, %d quintet'
        % (results['n_singlets'], results['n_triplets'], results['n_quintets']),
        '',
    ]
    col = '%5s %5s %18s %14s'
    table = [col % ('state', '2S+1', 'E_DK(Eh)', 'exc-DK(eV)'), '-' * 48]
    row = '%5d %5d %18.10f %14.4f'
    for state in states[:n_show]:
        table.append(row % (state['index'], state['mult'], state['E_DK'],
                            state['exc_DK_eV']))
    footer = [
        '',
        'CSF diagnostics: %s  (orthonormal %.2e, discarded cross-spin %.2e)'
        % ('PASS' if diag.get('pass') else 'FAIL',
           diag.get('csf_orthonormal_error', float('nan')),
           diag.get('discarded_cross_spin_block_max', float('nan'))),
    ]
    return '\n'.join(header + table + footer)


def format_qmrsf_log_table(results, max_states=10):
    """Render an aligned text table of the lowest states for the log.

    Columns: state index, the three total energies (Hartree), and the
    Epstein-Nesbet icPT2 excitation energy (eV, relative to state 0).
    """
    states = results.get('states', [])
    n_show = min(max_states, len(states))

    header_lines = [
        'QMRSF-icPT2 results',
        'reference (quintet ROHF) = %18.10f Hartree' % results['reference_energy'],
        'E_core (nuc + frozen)    = %18.10f Hartree' % results['ecore'],
        'window orbitals = %d   dressed roots = %d   (showing lowest %d)'
        % (results['n_window_orbitals'], results['n_states'], n_show),
        '',
    ]

    col = '%5s %18s %18s %18s %14s'
    rule = '-' * 76
    table = [
        col % ('state', 'E_CAS(Eh)', 'E_icPT2-EN(Eh)', 'E_icPT2-Dyall(Eh)', 'exc-EN(eV)'),
        rule,
    ]
    rowfmt = '%5d %18.10f %18.10f %18.10f %14.4f'
    for st in states[:n_show]:
        table.append(rowfmt % (
            st['index'], st['E_CAS'], st['E_icPT2_EN'], st['E_icPT2_Dyall'], st['exc_EN_eV'],
        ))

    return '\n'.join(header_lines + table)
