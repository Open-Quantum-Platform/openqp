"""OQP run type functions"""

import oqp
import oqp.library
from oqp.library.single_point import (
    SinglePoint, Gradient, Hessian, LastStep,
    BasisOverlap, NACME, NAC
)

from oqp.library.libscipy import StateSpecificOpt, MECIOpt, MECPOpt, MEP
from oqp.library.libdlfind import DLFindMin, DLFindTS, DLFindMECI
from oqp.library.libgeometric import (
    GeometricIRCOpt,
    GeometricMECIOpt,
    GeometricMECPOpt,
    GeometricOpt,
    GeometricTSOpt,
)


def compute_energy(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute properties
    compute_scf_prop(mol)


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
        else:
            raise ValueError(f'Unknown property: {prop}')


def compute_grad(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute gradient
    Gradient(mol).gradient()

    # compute properties
    compute_scf_prop(mol)

    # compute dftd4
    LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])

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
    sp = SinglePoint(mol)
    ref_energy = sp.reference()          # SCF один раз

    mol.data.set_tdhf_multiplicity(1)
    mol.singlet_energies = sp.excitation(ref_energy)

    mol.data['OQP::td_singlet_energies'] = mol.data['OQP::td_energies']
    mol.data['OQP::td_bvec_mo_s'] = mol.data['OQP::td_bvec_mo'].copy()

    mol.data.set_tdhf_multiplicity(3)
    mol.triplet_energies = sp.excitation(ref_energy)

    mol.data['OQP::td_triplet_energies'] = mol.data['OQP::td_energies']
    mol.data['OQP::td_bvec_mo_t'] = mol.data['OQP::td_bvec_mo'].copy()

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
    lib = mol.config['optimize']['lib']

    opt_lib = {
        'scipy': {
            'optimize': StateSpecificOpt,
            'meci': MECIOpt,
            'mecp': MECPOpt,
            'mep': MEP,
            'ts': None,
            'irc': None,
            'neb': None,
        },
        'dlfind': {
            'optimize': DLFindMin,
            'meci': DLFindMECI,
            'mecp': None,
            'mep': None,
            'ts': DLFindTS,
            'irc': None,
            'neb': None,
        },
        'geometric': {
            'optimize': GeometricOpt,
            'meci': GeometricMECIOpt,
            'mecp': GeometricMECPOpt,
            'mep': None,
            'ts': GeometricTSOpt,
            'irc': GeometricIRCOpt,
            'neb': None,
        },
    }

    if opt_lib[lib][runtype]:
        return opt_lib[lib][runtype](mol)

    else:
        raise ValueError(f'optimization library {lib} does not support runtype {runtype}')
