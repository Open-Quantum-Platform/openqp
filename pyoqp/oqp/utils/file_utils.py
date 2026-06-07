"""Helper utilities for file manipulation"""

import os
import time
import datetime
import numpy as np
from oqp.molden.moldenwriter import write_frequency
from oqp.periodic_table import SYMBOL_MAP, ELEMENTS_NAME
from oqp.utils.constants import ANGSTROM_TO_BOHR
from oqp.utils.mpi_utils import mpi_dump


def try_basis(basis, path=None, fallback='6-31g'):
    """try various basis file locations and return the matching one"""

    if path:
        basis_path = path
#    elif os.environ["OPENQP_ROOT"]:
#        basis_path = os.environ["OPENQP_ROOT"] + "/share/basis_sets"
    else:
        try:
            os.environ["OPENQP_ROOT"]
        except KeyError:
            os.environ["OPENQP_ROOT"] = os.path.abspath(os.path.dirname(__file__), os.pardir)
        basis_path = os.environ["OPENQP_ROOT"] + "/share/basis_sets"

    if not basis:
        basis = fallback

    tryfile = basis
    if os.path.isfile(tryfile):
        return tryfile

    tryfile = f'{basis_path}/{basis}'
    if os.path.isfile(tryfile):
        return tryfile

    tryfile = f'{basis_path}/{basis}.basis'
    if os.path.isfile(tryfile):
        return tryfile

    raise FileNotFoundError(f"Basis `{basis}` is not available")


def try_data_file(name):
    """Resolve a data file shipped under share/basis_sets (installed) or the
    source basis_sets/ tree (development)."""

    try:
        root = os.environ["OPENQP_ROOT"]
    except KeyError:
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

    for candidate in (
        name,
        os.path.join(root, "share", "basis_sets", name),
        os.path.join(root, "basis_sets", name),
    ):
        if os.path.isfile(candidate):
            return candidate

    raise FileNotFoundError(f"Data file `{name}` is not available")


def what_is_time():
    # This function return current time

    return datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')


def how_long(start, end):
    # This function calculate time between start and end

    walltime = end - start
    walltime = '%5d days %5d hours %5d minutes %5d seconds' % (
        int(walltime / 86400),
        int((walltime % 86400) / 3600),
        int(((walltime % 86400) % 3600) / 60),
        int(((walltime % 86400) % 3600) % 60))
    return walltime


def _to_yes_no(value):
    return 'yes' if bool(value) else 'no'
@mpi_dump


def dump_log(mol, title=None, section=None, info=None, must_print=False):
    # function to write information to main log
    logfile = mol.log
    method = mol.config['input']['method']
    basis = mol.config['input']['basis']
    natom = mol.data._data.mol_prop.natom
    charge = mol.data._data.mol_prop.charge
    functional = mol.config['input']['functional']
    if not functional:
        functional = 'hf'

    scf_type = mol.data._data.control.scftype
    scf_maxit = mol.data._data.control.maxit
    scf_mult = mol.data._data.mol_prop.mult
    scf_forced_attempt = mol.config['scf']['forced_attempt']
    scf_conv = mol.config['scf']['conv']
    scf_incre = mol.config['scf']['incremental']
    diis_type = mol.config['scf']['diis_type']
    vdiis_cdiis_switch = mol.config['scf']['vdiis_cdiis_switch']
    vdiis_vshift_switch = mol.config['scf']['vdiis_vshift_switch']
    vshift_cdiis_switch = mol.config['scf']['vshift_cdiis_switch']
    vshift = mol.config['scf']['vshift']

    td_type = mol.config['tdhf']['type']
    td_maxit = mol.config['tdhf']['maxit']
    td_maxit_zv = mol.config['tdhf']['maxit_zv']
    td_mult = mol.data._data.tddft.mult
    td_conv = mol.config['tdhf']['conv']
    td_nstate = mol.config['tdhf']['nstate']
    td_zvconv = mol.config['tdhf']['zvconv']
    td_nvdav = mol.config['tdhf']['nvdav']

    scftypes = {1: "rhf", 2: "uhf", 3: "rohf"}

    mode = 'a'
    loginfo = """
   ==============================================
   %s
   ==============================================
""" % title

    if section == 'start':
        mode = 'w'
        loginfo = """
   PyOQP started at %s

""" % what_is_time()

    if section == 'end':
        start = mol.start_time
        end = time.time()
        loginfo = '\n   %s\n   PyOQP terminated at %s in %s\n' % (title, what_is_time(), how_long(start, end))

    if section == 'guess':
        loginfo += """
   PyOQP guess type:                   %s
   PyOQP guess file:                   %s
   PyOQP guess alpha:                  %s
   PyOQP guess beta:                   %s
   PyOQP guess swapmo:                 %s

""" % (info['guess_type'], info['guess_file'], info['guess_alpha'], info['guess_beta'], info['guess_swapmo'])

    if section == 'basis_overlap':
        coord = mol.config['input']['system2']
        file2 = mol.config['guess']['file2']
        align_type = mol.config['nac']['align']

        if len(file2) > 0:
            file_type = file2
            coord_type = file2
        else:
            file_type = 'compute'
            if len(coord) == 1:
                coord_type = coord[0]
            else:
                coord_type = 'input'

        loginfo += """
   PyOQP previous coordinates          %s
   PyOQP previous data file            %s
   PyOQP align type                    %s

""" % (coord_type, file_type, align_type)

    if section == 'nacme':
        nac_type = mol.config['nac']['type']
        nac_dt = mol.config['nac']['dt']
        nac_dim = mol.config["tdhf"]["nstate"]

        loginfo += """
   PyOQP nac type                      %s
   PyOQP nacme delta                   %s
   PyOQP nac matrix dimension          %s

   Note:
        derivative coupling:      d_ij = (s_ij - s_ji) / delt
        non-adiabatic coupling:   h_ij = d_ij * e_ji
        
        d_ij is time derivative if delt in atomic unit of time
        d_ij is distance derivative if delt in atomic unit of distance (bohr)
        d_ij is unitless if delt is 1

""" % (nac_type, nac_dt, nac_dim)

    if section == 'input':
        loginfo += """
   PyOQP natom:                        %s
   PyOQP charge:                       %s
    
""" % (natom, charge)


    if section == 'symmetry':
        metadata = mol.symmetry_metadata if isinstance(getattr(mol, 'symmetry_metadata', {}), dict) else {}
        symmetry_status = metadata.get('status', 'disabled')
        loginfo += """
   PyOQP symmetry status:                      %s
   PyOQP symmetry requested point group:       %s
   PyOQP symmetry requested subgroup:          %s
   PyOQP symmetry detected point group:        %s
   PyOQP symmetry detected subgroup:           %s
   PyOQP symmetry label MO:                   %s
   PyOQP symmetry label states:               %s
   PyOQP symmetry label modes:                %s
   PyOQP symmetry use integral symmetry:      %s
   PyOQP symmetry use response symmetry:      %s
   PyOQP symmetry strict:                     %s
   PyOQP symmetry tolerance:                  %s
""" % (
            symmetry_status,
            metadata.get('requested_point_group', metadata.get('point_group', 'auto')),
            metadata.get('requested_subgroup', metadata.get('subgroup', 'auto')),
            metadata.get('detected_point_group', metadata.get('point_group', 'c1')),
            metadata.get('detected_subgroup', metadata.get('subgroup', 'c1')),
            _to_yes_no(metadata.get('label_mo', True)),
            _to_yes_no(metadata.get('label_states', True)),
            _to_yes_no(metadata.get('label_modes', True)),
            _to_yes_no(metadata.get('use_integral_symmetry', False)),
            _to_yes_no(metadata.get('use_response_symmetry', False)),
            _to_yes_no(metadata.get('strict', False)),
            metadata.get('tolerance', 1.0e-5),
        )

    if section in ['scf']:
        loginfo += """
   PyOQP method:                       %s
   PyOQP hf/functional:                %s
   PyOQP basis:                        %s
   PyOQP scf type:                     %s
   PyOQP scf maxit:                    %s
   PyOQP scf forced attempt:           %s
   PyOQP scf multiplicity:             %s
   PyOQP scf convergence:              %s
   PyOQP scf incremental:              %s
   PyOQP diis type:                    %s
   PyOQP vdiis/cdiis switch:           %s
   PyOQP vdiis/vshift switch:          %s
   PyOQP cdiis/vshift_switch           %s
   PyOQP vshift:                       %s
   
""" % (
            method, functional, basis, scftypes[scf_type], scf_maxit, scf_forced_attempt, scf_mult, scf_conv, scf_incre,
            diis_type,
            vdiis_cdiis_switch,
            vdiis_vshift_switch,
            vshift_cdiis_switch,
            vshift
        )

    if section == 'tdhf':
        loginfo += """
   PyOQP method:                       %s
   PyOQP functional:                   %s
   PyOQP td type:                      %s
   PyOQP td maxit:                     %s
   PyOQP td maxit z-vector:            %s
   PyOQP td multiplicity:              %s
   PyOQP td convergence:               %s
   PyOQP td number of states:          %s
   PyOQP td z-vector of convergence:   %s
   PyOQP td dimension of Davidson:     %s
    
""" % (method, functional, td_type, td_maxit, td_maxit_zv, td_mult, td_conv, td_nstate, td_zvconv, td_nvdav)

    if section == 'dftd':
        loginfo += """
   PyOQP dftd correction:                   %14s
   PyOQP dftd method:                       %14s
   PyOQP dftd functional:                   %14s

""" % (info['d4'], info['type'], functional)

    if section == 'energy':
        loginfo += '   PyOQP electronic energies\n'
        for n, energy in enumerate(info['el']):
            loginfo += f'   PyOQP state {n:<6} {energy:<16.8f}\n'

        d4 = float(info['d4'])
        loginfo += f'\n   PyOQP dftd correction {d4:<16.8f}\n\n'
        loginfo += '   PyOQP dispersion corrected energies\n'
        for n, energy in enumerate(mol.energies):
            loginfo += f'   PyOQP state {n:<6} {energy:<16.8f}\n'

    if section == 'grad':
        atoms = mol.get_atoms()
        loginfo += '   PyOQP electronic gradients\n'
        for n in info['grad_list']:
            grad = write_grad(atoms, info['el'][n])
            loginfo += f'   PyOQP state {n:<6}\n{grad}\n'

        d4 = info['d4']
        loginfo += f'\n   PyOQP dftd correction\n'
        loginfo += f'{write_grad(atoms, d4)}\n\n'
        loginfo += '   PyOQP dispersion corrected gradients\n'
        for n in info['grad_list']:
            grad = write_grad(atoms, mol.grads[n])
            loginfo += f'   PyOQP state {n:<6}\n{grad}\n'

    if section == 'opt':
        loginfo += """
   PyOQP follow state:                 %14s
   PyOQP energy shift:                 %14.6f %14.6f %s
   PyOQP rmsd step:                    %14.6f %14.6f %s
   PyOQP max step:                     %14.6f %14.6f %s
   PyOQP rmsd grad:                    %14.6f %14.6f %s
   PyOQP max grad:                     %14.6f %14.6f %s

""" % (
            info['istate'], info['de'], info['energy_shift'], np.abs(info['de']) <= info['energy_shift'],
            info['rmsd_step'], info['target_rmsd_step'], info['rmsd_step'] <= info['target_rmsd_step'],
            info['max_step'], info['target_max_step'], info['max_step'] <= info['target_max_step'],
            info['rmsd_grad'], info['target_rmsd_grad'], info['rmsd_grad'] <= info['target_rmsd_grad'],
            info['max_grad'], info['target_max_grad'], info['max_grad'] <= info['target_max_grad'],
        )

    if section == 'cons_sphere':
        loginfo += """
   PyOQP follow state:                 %14s
   PyOQP target step:                  %14.6f
   PyOQP energy shift:                 %14.6f %14.6f %s
   PyOQP rmsd step:                    %14.6f %14.6f %s
   PyOQP max step:                     %14.6f %14.6f %s
   PyOQP constraint step size:         %14.6f %14.6f %s
   PyOQP constraint rmsd grad:         %14.6f %14.6f %s
   PyOQP constraint max grad:          %14.6f %14.6f %s   
""" % (
            info['istate'], info['step_size'],
            info['de'], info['energy_shift'], np.abs(info['de']) <= info['energy_shift'],
            info['rmsd_step'], info['target_rmsd_step'], info['rmsd_step'] <= info['target_rmsd_step'],
            info['max_step'], info['target_max_step'], info['max_step'] <= info['target_max_step'],
            info['radius'], info['step_tol'], info['radius'] > info['step_tol'],
            info['rmsd_grad'], info['target_rmsd_grad'], info['rmsd_grad'] <= info['target_rmsd_grad'],
            info['max_grad'], info['target_max_grad'], info['max_grad'] <= info['target_max_grad'],
        )

    if section == 'penalty':
        loginfo += """
   PyOQP follow state:                 %14s %14s
   PyOQP meci search algorithm:        %14s
   PyOQP penalty sigma:                %14.6f
   PyOQP penalty alpha:                %14.6f
   PyOQP penalty increase:             %14.6f
   PyOQP energy shift:                 %14.6f %14.6f %s
   PyOQP energy gap:                   %14.6f %14.6f %s
   PyOQP rmsd step:                    %14.6f %14.6f %s
   PyOQP max step:                     %14.6f %14.6f %s
   PyOQP rmsd grad:                    %14.6f %14.6f %s
   PyOQP max grad:                     %14.6f %14.6f %s
       
""" % (
            info['istate'], info['jstate'],
            info['meci_search'],
            info['sigma'], info['alpha'], info['incre'],
            info['de'], info['energy_shift'], np.abs(info['de']) <= info['energy_shift'],
            info['gap'], info['energy_gap'], info['gap'] <= info['energy_gap'],
            info['rmsd_step'], info['target_rmsd_step'], info['rmsd_step'] <= info['target_rmsd_step'],
            info['max_step'], info['target_max_step'], info['max_step'] <= info['target_max_step'],
            info['rmsd_grad'], info['target_rmsd_grad'], info['rmsd_grad'] <= info['target_rmsd_grad'],
            info['max_grad'], info['target_max_grad'], info['max_grad'] <= info['target_max_grad'],
        )

    if section == 'ubp':
        loginfo += """
   PyOQP follow state:                 %14s %14s
   PyOQP meci search algorithm:        %14s
   PyOQP cgv norm:                     %14.6f
   PyOQP cgv orth:                     %14.6f
   PyOQP energy shift:                 %14.6f %14.6f %s
   PyOQP energy gap:                   %14.6f %14.6f %s
   PyOQP rmsd step:                    %14.6f %14.6f %s
   PyOQP max step:                     %14.6f %14.6f %s
   PyOQP rmsd grad:                    %14.6f %14.6f %s
   PyOQP max grad:                     %14.6f %14.6f %s

""" % (
            info['istate'], info['jstate'],
            info['meci_search'],
            info['norm'], info['orth'],
            info['de'], info['energy_shift'], np.abs(info['de']) <= info['energy_shift'],
            info['gap'], info['energy_gap'], info['gap'] <= info['energy_gap'],
            info['rmsd_step'], info['target_rmsd_step'], info['rmsd_step'] <= info['target_rmsd_step'],
            info['max_step'], info['target_max_step'], info['max_step'] <= info['target_max_step'],
            info['rmsd_grad'], info['target_rmsd_grad'], info['rmsd_grad'] <= info['target_rmsd_grad'],
            info['max_grad'], info['target_max_grad'], info['max_grad'] <= info['target_max_grad'],
        )

    if section == 'quad':
        loginfo += """
   PyOQP follow state:                 %14s %14s
   PyOQP mecp search algorithm:        %14s
   PyOQP energy shift:                 %14.6f %14.6f %s
   PyOQP energy gap:                   %14.6f %14.6f %s
   PyOQP rmsd step:                    %14.6f %14.6f %s
   PyOQP max step:                     %14.6f %14.6f %s
   PyOQP rmsd grad:                    %14.6f %14.6f %s
   PyOQP max grad:                     %14.6f %14.6f %s

""" % (
            info['istate'], info['jstate'] + info['nstate'],
            info['mecp_search'],
            info['de'], info['energy_shift'], np.abs(info['de']) <= info['energy_shift'],
            info['gap'], info['energy_gap'], np.abs(info['gap']) <= info['energy_gap'],
            info['rmsd_step'], info['target_rmsd_step'], info['rmsd_step'] <= info['target_rmsd_step'],
            info['max_step'], info['target_max_step'], info['max_step'] <= info['target_max_step'],
            info['rmsd_grad'], info['target_rmsd_grad'], info['rmsd_grad'] <= info['target_rmsd_grad'],
            info['max_grad'], info['target_max_grad'], info['max_grad'] <= info['target_max_grad'],
        )

    if section == 'mep':
        loginfo += """
   PyOQP MEP follow state:             %14s
   PyOPQ MEP opt steps:                %14s
   PyOQP MEP opt status:               %14s
   PyOQP MEP radius:                   %14.6f
   PyOQP MEP energy:                   %14.6f
   PyOQP MEP energy shift:             %14.6f
   
""" % (info['istate'], info['itr'], info['status'], info['radius'], info['energy'], info['de'])

    if section == 'num_nacv':
        ndim, dx, restart, jobs, nproc, threads = info
        loginfo = """
   PyOQP nac type                  %14s
   PyOQP number of displacements       %14s
   PyOQP size of displacements         %14s
   PyOQP calculation restart           %14s
   PyOQP number of nacme               %14s
   PyOQP number of processes           %14s
   PyOQP number of threads             %14s
   
""" % ('numerical', ndim, dx, restart, jobs, nproc, threads)

    if section == 'nacv_worker':
        order, idx, flag, timing = info
        start, end, rank, threads, host = timing
        loginfo = f'   PyOQP step: {order:<8} displacement: {idx:<8} {flag:<10} in {end - start:<16.0f} sec' \
                  f' from rank {rank:<3} with {threads:<3} threads on node {host}\n'

    if section == 'nacv':
        atoms = mol.get_atoms()
        states = mol.config['nac']['states']
        energies = mol.energies
        for ij in states:
            i, j = np.sort(ij)
            gap = energies[j] - energies[i]
            nac = write_grad(atoms, info[i - 1, j - 1])
            loginfo += f'   PyOQP NAC vector between state {i:<6} {j:<6} in Hartree/Bohr gap: {gap:16.8f}\n{nac}\n'

    if section == 'dcv':
        atoms = mol.get_atoms()
        states = mol.config['nac']['states']
        energies = mol.energies
        for ij in states:
            i, j = np.sort(ij)
            gap = energies[j] - energies[i]
            dc = write_grad(atoms, info[i - 1, j - 1])
            loginfo += f'   PyOQP DC vector between state {i:<6} {j:<6} in 1/Bohr gap: {gap:16.8f}\n{dc}\n'

    if section == 'nacm':
        for i in info:
            loginfo += ' '.join('%16.8f' % x for x in i) + '\n'

    if section == 'bp':
        atoms = mol.get_atoms()
        g1, g2, h, x, y, sx, sy, pitch, tilt, peak, bifu = info
        x = write_grad(atoms, x)
        y = write_grad(atoms, y)
        loginfo += f'   X vector\n{x}\n    Y vector\n{y}\n'
        loginfo += f'   average energy:    {sx:16.8f} * x + {sy:16.8f} * y\n'
        loginfo += f'   energy difference: {pitch:16.8f} * [(x^2 + y^2) + {tilt:16.8f} * (x^2 - y^2)] ** 0.5\n'
        loginfo += f'   P: {peak:16.8f} <1 peaked; >1 sloped\n'
        loginfo += f'   P: {bifu:16.8f} <1 bifurcating; >1 single-path\n'

    if section == 'read_hess':
        hess_file = mol.log.replace('.log', 'hess.json')
        loginfo = """
   PyOQP read hessian file            %14s
""" % hess_file

    if section == 'num_hess':
        state, ndim, dx, restart, jobs, nproc, threads = info
        loginfo = """
   PyOQP hessian type                  %14s
   PyOQP hessian follow state          %14s
   PyOQP number of displacements       %14s
   PyOQP size of displacements         %14s
   PyOQP calculation restart           %14s
   PyOQP number of grad                %14s
   PyOQP number of processes           %14s
   PyOQP number of threads             %14s

""" % ('numerical', state, ndim, dx, restart, jobs, nproc, threads)

    if section == 'hess_worker':
        order, idx, flag, timing = info
        start, end, rank, threads, host = timing
        loginfo = f'   PyOQP step: {order:<8} displacement: {idx:<8} {flag:<10} in {end - start:<16.0f} sec' \
                  f' from rank {rank:<3} with {threads:<3} threads on node {host}\n'

    if section == 'freq':
        ir = np.asarray(getattr(mol, 'infrared_intensities', []), dtype=float)
        raman = np.asarray(getattr(mol, 'raman_activities', []), dtype=float)
        if ir.size == len(info) and raman.size == len(info):
            loginfo += '   Mode       Frequency(cm-1)      IR(km/mol)        Raman(activity)\n'
            for n, f in enumerate(info):
                loginfo += f'   {n + 1:4d} {f:20.2f} {ir[n]:16.6f} {raman[n]:20.6f}\n'
        else:
            for n, f in enumerate(info):
                loginfo += f'   PyOQP freq {n + 1}:  {f:12.2f}\n'

    if section == 'freq_modes':
        atoms, freqs, modes = info
        atoms = np.asarray(atoms, dtype=int)
        freqs = np.asarray(freqs, dtype=float)
        modes = np.asarray(modes, dtype=float)
        natom = len(atoms)
        loginfo += """
   Normal mode eigenvectors (Cartesian, mass-unweighted)
   Frequencies -- values are in cm^-1; X/Y/Z columns are normal-mode components.
"""
        for start in range(0, len(freqs), 1):
            stop = min(start + 1, len(freqs))
            block = range(start, stop)
            loginfo += '\n                 ' + ''.join(f'{mode_index + 1:>12d}' for mode_index in block) + '\n'
            loginfo += '   Frequencies --' + ''.join(f'{freqs[mode_index]:12.4f}' for mode_index in block) + '\n'
            loginfo += '     Atom AN      ' + ''.join(f'{axis:>12s}' for _mode_index in block for axis in ('X', 'Y', 'Z')) + '\n'
            for atom_index, atomic_number in enumerate(atoms):
                symbol = ELEMENTS_NAME[atomic_number] if 0 <= atomic_number < len(ELEMENTS_NAME) else str(atomic_number)
                row = f'   {atom_index + 1:6d} {atomic_number:2d} {symbol:>2s}'
                for mode_index in block:
                    vec = modes[mode_index].reshape((natom, 3))[atom_index]
                    row += ''.join(f'{component:12.8f}' for component in vec)
                loginfo += row + '\n'

    if section == 'thermo':
        temp = info['temp']
        mass = info['mass']
        rc = info['rc']
        rt = info['rt']
        el = info['el']
        zpe = info['zpe']
        u_trans = info['u_trans']
        u_rot = info['u_rot']
        u_vib = info['u_vib']
        pv = info['pv']
        st_el = info['st_el']
        st_trans = info['st_trans']
        st_rot = info['st_rot']
        st_vib = info['st_vib']

        u_el = u_trans + u_rot + u_vib + zpe
        u = u_el + el
        h_el = u_el + pv
        h = h_el + el
        st = st_el + st_trans + st_rot + st_vib
        g_el = h_el + st
        g = g_el + el

        loginfo += """
   temperature K:                    %16.2f
   pressure atm:                     %16.2f
   total mass amu:                   %16.2f
   rotational constant cm-1:   %12.4f %12.4f %12.4f 
   rotational temperature K:   %12.4f %12.4f %12.4f 

   ====================================================
   summary of internal energy (U)
   ====================================================
   U = E(el) + E(trans) + E(rot) + E(vib) + E(ZPE)

   E(el) electronic energy:          %16.8f
   E(trans) translational energy:    %16.8f
   E(rot) rotational energy:         %16.8f
   E(vib) vibrational energy:        %16.8f
   E(ZPE) zero-point energy:         %16.8f
   ----------------------------------------------------
   total correction to internal:     %16.8f   
   total internal energy:            %16.8f
   
   ====================================================
   summary of enthalpy (H)
   ====================================================
   H = U + pV
   
   E(el) electronic energy:          %16.8f
   U - E(el) correction:             %16.8f
   pV enthalpy correction:           %16.8f
   ----------------------------------------------------
   total correction to enthalpy:     %16.8f
   total enthalpy:                   %16.8f
   
   ====================================================
   summary of entropy (TS)
   ====================================================
   TS = TS(el) + TS(trans) + TS(rot) + TS(vib)
   
   TS(el) electronic entropy:        %16.8f
   TS(trans) translational entropy:  %16.8f
   TS(rot) rotational entropy:       %16.8f
   TS(vib) vibrational entropy:      %16.8f
   ----------------------------------------------------
   total entropy:                    %16.8f   
   
   ====================================================
   summary of Gibbs free energy (G)
   ====================================================
   G = H - TS
   
   E(el) electronic energy:          %16.8f
   H - E(el) correction:             %16.8f
   TS total entropy:                 %16.8f
   ----------------------------------------------------
   total correction to Gibbs:        %16.8f
   total Gibbs free energy:          %16.8f

""" % (
            temp, 1.0, mass,
            rc[0], rc[1], rc[2],
            rt[0], rt[1], rt[2],
            el, u_trans, u_rot, u_vib, zpe,
            u_el, u,
            el, u_el, pv,
            h_el, h,
            st_el, st_trans, st_rot, st_vib,
            st,
            el, h_el, st,
            g_el, g
        )

    with open(logfile, mode) as out:
        out.write(loginfo)

@mpi_dump
def dump_data(mol, data, title=None, fpath='.'):
    # function to write data in specific logs
    if title == 'ENERGY':
        energies, title = data
        filename = 'energies'

        if title:
            filename = f'{title}.{filename}'

        np.savetxt(f'{fpath}/{filename}', np.array(energies).reshape((-1, 1)), fmt='%24.16f')

    if title == 'GRADIENT':
        grads, title, grad_list = data
        filename = 'grad'

        if title:
            filename = f'{title}.{filename}'

        for i in grad_list:
            np.savetxt(f'{fpath}/{filename}_{i}', grads[i], fmt='%24.16f')

    if title == 'NACME':
        nacme, title = data
        filename = 'nacme'

        if title:
            filename = f'{title}.{filename}'

        np.savetxt(f'{fpath}/{filename}', nacme, fmt='%24.16f')

    if title == 'DCME':
        dcme, title = data
        filename = 'dcme'

        if title:
            filename = f'{title}.{filename}'

        np.savetxt(f'{fpath}/{filename}', dcme, fmt='%24.16f')

    if title == 'NACV':
        mol, nacv, title = data
        filename = 'nac'
        states = mol.config['nac']['states']

        if title:
            filename = f'{title}.{filename}'

        for ij in states:
            i, j = np.sort(ij)
            filename_ext = f'{filename}_{i}_{j}'
            np.savetxt(f'{fpath}/{filename_ext}', nacv[i - 1, j - 1], fmt='%24.16f')

    if title == 'BP':
        mol, g1, g2, h, x, y, i, j = data
        molden = write_frequency(mol, np.array([101, 201, 301, 401, 300, 400]), np.array([g1, g2, g2 - g1, h, x, y]))

        with open(f'{fpath}/{mol.project_name}.{i}_{j}.gh.molden', 'w') as out:
            out.write(molden)

    if title == 'OPTIMIZATION':
        itr, atoms, coordinates, energy, de, rmsd_step, max_step, rmsd_grad, max_grad = data

        if itr == 1:
            mode = 'w'
            status = """%5s %16s %16s %14s %14s %14s %14s
----------------------------------------------------------------------------------------------------------------------
""" % (
                'Step', 'Energy', 'Shift', 'RMSD Step', 'Max Step', 'RMSD Grad', 'Max Grad'
            )
        else:
            mode = 'a'
            status = ''

        xyz = write_xyz(atoms, coordinates, (itr, energy))
        status += '%5s %16.8f %16.8f %14.6f %14.6f %14.6f %14.6f\n' % (
            itr, energy, de, rmsd_step, max_step, rmsd_grad, max_grad,
        )

        with open(f'{fpath}/opt.xyz', 'w') as out:
            out.write(xyz)

        with open(f'{fpath}/opt_geom.xyz', mode) as out:
            out.write(xyz)

        with open(f'{fpath}/opt_status.txt', mode) as out:
            out.write(status)

    if title == 'MECI':
        itr, atoms, coordinates, energy, de, gap, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df, max_df = data
        if itr == 1:
            mode = 'w'
            status = """%5s %16s %16s %16s %14s %14s %14s %14s %14s %14s
--------------------------------------------------------------------------------------------------------------------------------------------------
""" % (
                'Step', 'Energy', 'Shift', 'Gap', 'RMSD Step', 'Max Step', 'RMSD dP', 'Max dP', 'RMSD dE', 'Max dE'
            )
        else:
            mode = 'a'
            status = ''

        xyz = write_xyz(atoms, coordinates, (itr, energy))
        status += '%5s %16.8f %16.8f %16.8f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f\n' % (
            itr, energy, de, gap, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df, max_df
        )

        with open(f'{fpath}/opt.xyz', 'w') as out:
            out.write(xyz)

        with open(f'{fpath}/opt_geom.xyz', mode) as out:
            out.write(xyz)

        with open(f'{fpath}/opt_status.txt', mode) as out:
            out.write(status)

    if title == 'MECP':
        itr, atoms, coordinates, energy, de, gap, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df2, max_df2 = data
        if itr == 1:
            mode = 'w'
            status = """%5s %16s %16s %16s %14s %14s %14s %14s %14s %14s
--------------------------------------------------------------------------------------------------------------------------------------------------
""" % (
                'Step', 'Energy', 'Shift', 'Gap', 'RMSD Step', 'Max Step', 'RMSD G', 'Max G', 'RMSD dG', 'Max dG'
            )
        else:
            mode = 'a'
            status = ''

        xyz = write_xyz(atoms, coordinates, (itr, energy))
        status += '%5s %16.8f %16.8f %16.8f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f\n' % (
            itr, energy, de, gap, rmsd_step, max_step, rmsd_grad, max_grad, rmsd_df2, max_df2
        )

        with open(f'{fpath}/opt.xyz', 'w') as out:
            out.write(xyz)

        with open(f'{fpath}/opt_geom.xyz', mode) as out:
            out.write(xyz)

        with open(f'{fpath}/opt_status.txt', mode) as out:
            out.write(status)

    if title == 'CONS_SPHERE':
        itr, atoms, coordinates, energy, de, rmsd_step, max_step, rmsd_grad, max_grad, radius, dist = data
        if itr == 1:
            mode = 'w'
            status = """%5s %16s %16s %14s %14s %14s %14s %14s %14s
------------------------------------------------------------------------------------------------------------------------------------
""" % (
                'Step', 'Energy', 'Shift', 'RMSD Step', 'Max Step', 'RMSD Grad', 'Max Grad', 'Radius', 'Dist',
            )
        else:
            mode = 'a'
            status = ''

        xyz = write_xyz(atoms, coordinates, (itr, energy))
        status += '%5s %16.8f %16.8f %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f\n' % (
            itr, energy, de, rmsd_step, max_step, rmsd_grad, max_grad, radius, dist
        )

        with open(f'{fpath}/opt.xyz', 'w') as out:
            out.write(xyz)

        with open(f'{fpath}/opt_geom.xyz', mode) as out:
            out.write(xyz)

        with open(f'{fpath}/opt_status.txt', mode) as out:
            out.write(status)

    if title == 'MEP':
        itr, atoms, coordinates, energy, de = data
        if itr == 0:
            mode = 'w'
            status = """%5s %16s %16s
----------------------------------------
%5s %16.8f %16.8f
""" % ('Step', 'Energy', 'Shift', itr, energy, de)

        else:
            mode = 'a'
            status = '%5s %16.8f %16.8f\n' % (itr, energy, de)

        xyz = write_xyz(atoms, coordinates, (itr, energy))

        with open(f'{fpath}/mep.xyz', 'w') as out:
            out.write(xyz)

        with open(f'{fpath}/mep_geom.xyz', mode) as out:
            out.write(xyz)

        with open(f'{fpath}/mep_status.txt', mode) as out:
            out.write(status)

    if title == 'NUM_NACV':
        order, idx, norm, timing = data
        start, end, rank, threads, node = timing

        if order == 1:
            mode = 'w'
            status = """%8s %8s %8s %8s %8s %16s %16s
----------------------------------------------------------------------------------
%8s %8s %8s %8s %8s %16d %16.8f
""" % ('Step', 'Index', 'Rank', 'Threads', 'Node', 'Time', 'Norm', order, idx, rank, threads, node, end - start, norm)

        else:
            mode = 'a'
            status = '%8s %8s %8s %8s %8s %16d %16.8f\n' % (order, idx, rank, threads, node, end - start, norm)

        with open(f'{fpath}/nacv.status', mode) as out:
            out.write(status)

    if title == 'NUM_HESS':
        order, idx, norm, timing = data
        start, end, rank, threads, node = timing

        if order == 1:
            mode = 'w'
            status = """%8s %8s %8s %8s %8s %16s %16s
----------------------------------------------------------------------------------
%8s %8s %8s %8s %8s %16d %16.8f
""" % ('Step', 'Index', 'Rank', 'Threads', 'Node', 'Time', 'Norm', order, idx, rank, threads, node, end - start, norm)

        else:
            mode = 'a'
            status = '%8s %8s %8s %8s %8s %16d %16.8f\n' % (order, idx, rank, threads, node, end - start, norm)

        with open(f'{fpath}/hess.status', mode) as out:
            out.write(status)

    if title == 'FREQ':
        mol, freqs, modes = data
        molden = write_frequency(mol, freqs, modes)

        with open(f'{fpath}/{mol.project_name}.freq.molden', 'w') as out:
            out.write(molden)


def write_xyz(atoms, coord, info):
    # coord in Bohr
    coord = coord.reshape((-1, 3))
    atoms = np.asarray(atoms).reshape(-1)
    natom = len(coord)
    xyz = '%s\nGeom %s\n' % (natom, ' '.join([str(x) for x in info]))
    for n, line in enumerate(coord):
        a = np.asarray(atoms[n]).reshape(-1)[0]
        x, y, z = line[0: 3]
        xyz += '%-5s %24.16f %24.16f %24.16f\n' % (
            ELEMENTS_NAME[SYMBOL_MAP[int(a)]],
            x * ANGSTROM_TO_BOHR,
            y * ANGSTROM_TO_BOHR,
            z * ANGSTROM_TO_BOHR
        )

    return xyz


def write_grad(atoms, grad):
    # grad in Hartree/Bohr
    grad = grad.reshape((-1, 3))
    atoms = np.asarray(atoms).reshape(-1)
    xyz = ''
    for n, line in enumerate(grad):
        a = atoms[n]
        x, y, z = line[0: 3]
        xyz += '%5s %16.8f %16.8f %16.8f\n' % (ELEMENTS_NAME[SYMBOL_MAP[int(a)]], x, y, z)

    return xyz


def write_config(config):
    input_file = ''
    input_dict = {}
    for section in config.keys():
        if section == 'test':
            continue

        input_file += f'[{section}]\n'
        input_dict[section] = {}
        for key, value in config[section].items():
            if not value:
                continue

            if isinstance(value, list):
                if isinstance(value[0], list):
                    value = ','.join(['%s %s' % (x[0], x[1]) for x in value])
                elif isinstance(value[0], int):
                    value = ','.join(['%s' % x for x in value])
                elif isinstance(value[0], str):
                    value = value[0]
                else:
                    continue

            input_file += f'{key}={value}\n'
            input_dict[section][key] = str(value)

        input_file += '\n'

    return input_file, input_dict
