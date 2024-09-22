"""OQP run type functions"""

import oqp
import oqp.library
from oqp.library.single_point import (
    SinglePoint, Gradient, Hessian, LastStep,
    BasisOverlap, NACME, NAC
)

from oqp.library.libscipy import StateSpecificOpt, MECIOpt, MECPOpt, MEP
from oqp.library.libdlfind import DLFindMin, DLFindTS, DLFindMECI


def prep_guess(mol):
    oqp.library.ints_1e(mol)
    oqp.library.guess(mol)


def compute_energy(mol):
    # prepare guess orbital
    prep_guess(mol)

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
    # prepare guess orbital
    prep_guess(mol)

    # compute energy
    SinglePoint(mol).energy()

    # compute gradient
    Gradient(mol).gradient()

    # compute dftd4
    LastStep(mol).compute(mol, grad_list=mol.config['properties']['grad'])

def compute_nacme(mol):
    # prepare guess orbital
    prep_guess(mol)

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
    # prepare guess orbital
    prep_guess(mol)

    # compute energy
    SinglePoint(mol).energy()

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute nacme
    NAC(mol).nac()

def compute_soc(mol):
    pass


def compute_hess(mol):
    # prepare guess orbital
    prep_guess(mol)

    # compute energy
    SinglePoint(mol).energy()

    # compute dftd4
    LastStep(mol).compute(mol)

    # compute hessian
    Hessian(mol).hessian()


def compute_thermo(mol):
    Hessian(mol).hessian()


def compute_geom(mol):
    # prepare guess orbital
    prep_guess(mol)

    # initialize optimizer
    optimizer = get_optimizer(mol)

    # optimize coordinates
    optimizer.optimize()


def compute_properties(mol):
    # prepare guess orbital
    prep_guess(mol)

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
    # prepare guess orbital
    prep_guess(mol)

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
            'neb': None,
        },
        'dlfind': {
            'optimize': DLFindMin,
            'meci': DLFindMECI,
            'mecp': None,
            'mep': None,
            'ts': DLFindTS,
            'neb': None,
        },
    }

    if opt_lib[lib][runtype]:
        return opt_lib[lib][runtype](mol)

    else:
        raise ValueError(f'optimization library {lib} does not support runtype {runtype}')
