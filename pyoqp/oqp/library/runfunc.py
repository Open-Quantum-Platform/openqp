"""OQP run type functions"""

import oqp
import oqp.library
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
from oqp.library.libdlfind import DLFindMin, DLFindTS, DLFindMECI, DLFindQMMM
from oqp.library.namd import NAMD
#from oqp.library.libopenmm import QMMM_MD


def compute_namd(mol):
    # Tully fewest-switches surface-hopping nonadiabatic molecular dynamics.
    # Gas-phase or QM/MM (electrostatic ESPF embedding + OpenMM MM region).
    if mol.config['input'].get('qmmm_flag'):
        from oqp.library.namd import NAMD_QMMM
        NAMD_QMMM(mol).run()
    else:
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
    NACME(mol).nacme()


def compute_nac(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute nacme
    NAC(mol).nac()

def compute_soc(mol):
    pass


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

    # compute gradient
    Gradient(mol).gradient()

    # compute dftd4
    LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])

    # compute nac or nacme
    nac_type = mol.config['properties']['nac']
    if nac_type == 'nac':
        NAC(mol).nac()
    elif nac_type == 'nacme':
        NACME(mol).nacme()
    else:
        pass

    # compute soc
    soc_type = mol.config['properties']['soc']
    if soc_type:
        pass
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

    # compute soc
    soc_type = mol.config['properties']['soc']
    if soc_type:
        pass
    else:
        pass


def get_optimizer(mol):
    runtype = mol.config['input']['runtype']
    lib = mol.config['optimize']['lib']

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

    # QM/MM dispatch (wired up later with the qmmm_flag input option):
    #    if runtype == 'optimize' and mol.config['input']['qmmm_flag']:
    #       if lib == 'dlfind': opt_lib[lib][runtype] = DLFindQMMM
    #       elif lib == 'scipy': opt_lib[lib][runtype] = QMMMOpt

    if opt_lib[lib].get(runtype):
        return opt_lib[lib][runtype](mol)

    else:
        raise ValueError(f'optimization library {lib} does not support runtype {runtype}')
