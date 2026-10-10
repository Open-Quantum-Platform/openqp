"""OQP run type functions"""
import os

import numpy as np

import oqp
import oqp.library
from oqp.utils.tb_backends import is_tb_method, tb_section_name
from oqp.library.single_point import (
    SinglePoint, Gradient, Hessian, LastStep,
    BasisOverlap, NACME, NAC
)

from oqp.library.libscipy import StateSpecificOpt, MECIOpt, MECPOpt, MEP, QMMMOpt
from oqp.library.libgeometric import (
    GeometricIRCOpt,
    GeometricMECIOpt,
    GeometricMECPOpt,
    GeometricNEBOpt,
    GeometricOpt,
    GeometricTSOpt,
)
from oqp.library.liboqp import (
    OQPOpt, OQPTSOpt, OQPMECIOpt, OQPMECPOpt, OQPTCIOpt,
    OQPNEBOpt, OQPIRCOpt, OQPMEPOpt,
)
#from oqp.library.libopenmm import QMMM_MD


def compute_namd(mol):
    # Tully fewest-switches surface-hopping nonadiabatic molecular dynamics.
    # Gas-phase or QM/MM (electrostatic ESPF embedding + OpenMM MM region).
    qmmm = mol.config['input'].get('qmmm_flag')
    soc = mol.config['md'].get('soc')
    soc_basis = str(mol.config['md'].get('soc_basis', 'adiabatic')).lower()
    if soc and soc_basis not in ('adiabatic', 'mch'):
        raise ValueError("[md] soc_basis must be 'adiabatic' or 'mch'")
    if qmmm and soc:
        if soc_basis == 'mch':
            from oqp.library.namd import NAMD_SOC_MCH_QMMM
            NAMD_SOC_MCH_QMMM(mol).run()
        else:
            from oqp.library.namd import NAMD_SOC_QMMM
            NAMD_SOC_QMMM(mol).run()
    elif qmmm:
        from oqp.library.namd import NAMD_QMMM
        NAMD_QMMM(mol).run()
    elif soc:
        if soc_basis == 'mch':
            from oqp.library.namd import NAMD_SOC_MCH
            NAMD_SOC_MCH(mol).run()
        else:
            from oqp.library.namd import NAMD_SOC
            NAMD_SOC(mol).run()
    else:
        from oqp.library.namd import NAMD
        NAMD(mol).run()


def compute_energy(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute properties
    compute_scf_prop(mol)

    # re-save mol data so property results (e.g. OQP::nmr_shielding) reach the
    # JSON; the save inside SinglePoint.energy() runs before properties exist
    if mol.config['guess']['save_mol']:
        mol.save_data()


def compute_scf_prop(mol):
    # compute HF/DFT properties
    properties = mol.config["properties"]["scf_prop"]
    for prop in properties:
        if prop == 'el_mom':
            oqp.electric_moments(mol)
        elif prop == 'mulliken':
            oqp.mulliken(mol)
        elif prop == 'lowdin':
            oqp.lowdin(mol)
        elif prop == 'resp':
            oqp.resp_charges(mol)
        elif prop == 'nmr':
            scf_type = mol.config.get("scf", {}).get("type", "rhf")
            if isinstance(scf_type, str):
                scf_type = scf_type.lower()

            nmr_gauge = mol.config.get("properties", {}).get("nmr_gauge", "cgo")
            if isinstance(nmr_gauge, str):
                nmr_gauge = nmr_gauge.lower()
            if nmr_gauge == "cgo":
                if scf_type in ("uhf", "rohf"):
                    raise NotImplementedError(
                        "CGO NMR shielding supports closed-shell RHF references only. "
                        "Use properties.nmr_gauge=giao for open-shell (UHF/ROHF) NMR."
                    )
                oqp.nmr_shielding(mol)
            elif nmr_gauge == "giao":
                oqp.nmr_giao_shielding(mol)
            else:
                raise ValueError(
                    f"Unknown NMR gauge formulation {nmr_gauge!r}; expected 'cgo' or 'giao'"
                )
        else:
            raise ValueError(f'Unknown property: {prop}')

#    if 'resp' not in properties and mol.config['input']['qmmm_flag']:
#        oqp.resp_charges(mol)



def compute_grad(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute gradient
    Gradient(mol).gradient()

    # compute properties
    compute_scf_prop(mol)
#    if mol.config['input']['qmmm_flag']:
#       oqp.resp_charges(mol)

    # compute dftd4
    LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])

def compute_md(mol):

    # prepare guess orbital
    prep_guess(mol)
    
    #Run MD
#    qmmm_md = QMMM_MD(mol)
    qmmm_md.run_md()

def compute_nacme(mol):
    # compute reference energy
    sp = SinglePoint(mol)
    ref_energy = sp.reference()

    # compute mo overlap and apply alignment
    BasisOverlap(mol).overlap()

    # compute excitation energy
    sp.excitation(ref_energy)

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute nacme
    # Numerical-NAC displacement workers must compare the same physical root
    # on the +dx and -dx sides even if the adiabatic solver exchanges their
    # energy order.  Ordinary NACME/MD steps retain energy-root order.
    spatial_worker = os.environ.get('OQP_NUM_NAC_WORKER', '') == '1'
    NACME(mol).nacme(
        reorder_x=spatial_worker,
        normalize_retained=spatial_worker,
    )


def compute_nac(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute nacme
    NAC(mol).nac()

def compute_soc(mol):
    if is_tb_method(mol.config['input']['method']):
        if tb_section_name(mol.config) == 'xtb':
            from oqp.library.openqp_xtb import xtb_soc as tb_soc
        else:
            from oqp.library.openqp_dftb import dftb_soc as tb_soc
        tb_soc(mol)
        LastStep(mol).compute(mol)
        return

    sp = SinglePoint(mol)
    ref_energy = sp.reference()          # SCF один раз

    td = mol.config['tdhf']
    saved_mult = int(td.get('multiplicity', 1))
    saved_nstate = int(td.get('nstate', 1))
    ns = int(td.get('nstate_s', 0)) or saved_nstate
    nt = int(td.get('nstate_t', 0)) or saved_nstate

    def select_manifold(mult, nstate):
        # Keep config and native state synchronized so Python/native log labels
        # describe the manifold that is actually being calculated.
        td['multiplicity'] = mult
        td['nstate'] = nstate
        mol.data.set_tdhf_multiplicity(mult)
        mol.data.set_tdhf_nstate(nstate)
        sp.nstate = nstate
        mol.data['OQP::state_sign'] = np.ones(nstate)

    try:
        select_manifold(1, ns)
        mol.singlet_energies = sp.excitation(ref_energy)
        mol.data['OQP::td_singlet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_s'] = mol.data['OQP::td_bvec_mo'].copy()

        select_manifold(3, nt)
        mol.triplet_energies = sp.excitation(ref_energy)
        mol.data['OQP::td_triplet_energies'] = mol.data['OQP::td_energies'].copy()
        mol.data['OQP::td_bvec_mo_t'] = mol.data['OQP::td_bvec_mo'].copy()
    finally:
        select_manifold(saved_mult, saved_nstate)

    oqp.soc_mrsf(mol)

    LastStep(mol).compute(mol) 

def compute_hess(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute hessian
    Hessian(mol).hessian()


def compute_thermo(mol):
    Hessian(mol).hessian()


def compute_geom(mol):
    # initialize optimizer
    optimizer = get_optimizer(mol)

    # optimize coordinates
    optimizer.optimize()

    # compute properties
    compute_scf_prop(mol)


def compute_properties(mol):
    # compute reference energy
    sp = SinglePoint(mol)
    ref_energy = sp.reference()

    # compute mo overlap and apply alignment
    BasisOverlap(mol).overlap()

    # compute excitation energy
    sp.excitation(ref_energy)

    # Transport the central-geometry response-vector gauge before gradients,
    # numerical-NAC workers, and the checkpoint JSON are produced.  PyRAI2MD
    # feeds that JSON into the next MD step; aligning only inside displaced
    # workers leaves the time-series phase undefined.
    nacme = NACME(mol)
    nacme.align_x()

    # compute gradient
    Gradient(mol).gradient()

    # compute dftd4
    LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])

    # compute nac or nacme
    nac_type = mol.config['properties']['nac']
    if nac_type == 'nac':
        NAC(mol).nac()
    elif nac_type == 'nacme':
        nacme.nacme(align=False)
    else:
        pass

def compute_data(mol):
    # compute reference energy
    SinglePoint(mol).energy()

    # compute gradient
    Gradient(mol).gradient()

    # compute dftd4
    LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])

    # compute nac or nacme
    nac_type = mol.config['properties']['nac']
    if nac_type == 'nac':
        NAC(mol).nac()
    else:
        pass



def get_optimizer(mol):
    runtype = mol.config['input']['runtype']
    # Normal parsed inputs receive the schema default, but scripting callers
    # and lightweight embedding APIs may provide a partial configuration.
    # Keep the native backend genuinely optional in user input: omitting both
    # [optimize] and optimize.lib must still select oqp.
    optimize_config = mol.config.setdefault('optimize', {})
    if 'lib' not in optimize_config:
        istate_omitted = 'istate' not in optimize_config
        from oqp.molecule.oqpdata import OQP_CONFIG_SCHEMA
        from oqp.utils.input_parser import schema_section_defaults
        defaults = schema_section_defaults(OQP_CONFIG_SCHEMA, 'optimize')
        for key, value in defaults.items():
            optimize_config.setdefault(key, value)
        # The schema default targets the first response state.  Lightweight
        # scripting configurations do not pass through schema normalization,
        # whose validator treats an omitted selector as the ground state for a
        # ground-state HF/DFT optimization.  Single-state CASPT2/MRMP2 likewise
        # publish only energy index zero.  Preserve those public behaviours
        # when materializing the remaining native defaults here.
        method = str(mol.config.get('input', {}).get('method', 'hf')).lower()
        if istate_omitted and method in {'hf', 'caspt2', 'mrmp2'}:
            optimize_config['istate'] = 0
    lib = str(optimize_config['lib']).strip().lower()

    # BaekA generalizes the adjacent-gap adaptive penalty from two to N states.
    # It is selected as a MECI search algorithm rather than as a separate
    # backend; the historical three-state ``tci`` spelling remains below.
    meci_search = str(
        optimize_config.get('meci_search', 'auto')
    ).strip().lower()
    if lib == 'oqp' and runtype == 'meci' and meci_search == 'auto':
        # Lazy import preserves compatibility with lightweight dispatcher test
        # doubles and keeps the automatic orchestration native-only.
        from oqp.library.liboqp import OQPAutoMECIOpt
        return OQPAutoMECIOpt(mol)
    if lib == 'oqp' and runtype == 'mecp' and \
            str(optimize_config.get('mecp_search', 'auto')).strip().lower() in ('auto', 'sqp'):
        # SQP supplies its own KKT step control instead of an objective for the
        # native engine, so it is dispatched like BaekA rather than through
        # OQPMECPOpt.
        from oqp.library.libscipy import SQPMECPOpt
        return SQPMECPOpt(mol)

    if lib == 'oqp' and runtype == 'meci' and meci_search == 'baeka':
        # Lazy import keeps lightweight dispatcher test doubles compatible;
        # normal OpenQP installations load the real native class here.
        from oqp.library.liboqp import OQPBaekAOpt
        return OQPBaekAOpt(mol)

    opt_lib = {
        'scipy': {
            'optimize': StateSpecificOpt,
            'meci': MECIOpt,
            'mecp': MECPOpt,
            'mep': MEP,
            'ts': None,
            'tci': None,
            'irc': None,
            'neb': None,
        },
        'geometric': {
            'optimize': GeometricOpt,
            'meci': GeometricMECIOpt,
            'mecp': GeometricMECPOpt,
            'mep': None,
            'ts': GeometricTSOpt,
            'tci': None,
            'irc': GeometricIRCOpt,
            'neb': GeometricNEBOpt,
        },
        'oqp': {
            'optimize': OQPOpt,
            'meci': OQPMECIOpt,
            'mecp': OQPMECPOpt,
            'mep': OQPMEPOpt,
            'ts': OQPTSOpt,
            'tci': OQPTCIOpt,
            'irc': OQPIRCOpt,
            'neb': OQPNEBOpt,
        },
    }

    if opt_lib[lib].get(runtype):
        return opt_lib[lib][runtype](mol)

    else:
        raise ValueError(f'optimization library {lib} does not support runtype {runtype}')
