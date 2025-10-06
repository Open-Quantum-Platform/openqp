"""OQP run type functions"""

import oqp
import oqp.library
from oqp.library.single_point import (
    SinglePoint, Gradient, Hessian, LastStep,
    BasisOverlap, NACME, NAC
)

from oqp.library.libscipy import StateSpecificOpt, MECIOpt, MECPOpt, MEP
from oqp.library.libdlfind import DLFindMin, DLFindTS, DLFindMECI


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

def compute_dyson(mol):
    """
    Compute Dyson orbitals with minimal overhead:
    1. SCF (ROHF) calculation
    2. TDHF (MRSF-TDDFT) calculation  
    3. Z-vector calculation (only what's needed for Dyson)
    4. Dyson orbital calculation
    """
    # Step 1: Compute energy (SCF + TDHF)
    SinglePoint(mol).energy()
    
    # Step 2: Compute Z-vector (not full gradient)
    # This mirrors the gradient workflow but only computes Z-vector
    td = mol.config["tdhf"]["type"]
    target_state = mol.config.get('dyson', {}).get('target_state', 1)
    
    if target_state == 0:
        target_state = mol.config.get('tdhf', {}).get('target', 1)
    
    mol.data.set_tdhf_target(target_state)
    
    zvec_func = {
        'rpa': oqp.tdhf_z_vector,
        'tda': oqp.tdhf_z_vector,
        'sf': oqp.tdhf_sf_z_vector,
        'mrsf': oqp.tdhf_mrsf_z_vector,
    }
    
    # Compute Z-vector
    zvec_func[td](mol)
    
    # Step 3: Compute Dyson orbitals (will build relaxed matrices internally)
    if td in ['mrsf', 'sf']:
        oqp.dyson_mrsf(mol)
    else:
        oqp.dyson_orbital_calculation(mol)
    
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


def compute_dyson(mol):
    # compute energy
    SinglePoint(mol).energy()

    # compute gradient for relaxed density
    Gradient(mol).gradient()

    # compute dyson orbitals
    dyson_config = mol.config.get('dyson', {})
    compute_ip = dyson_config.get('ip', True)
    compute_ea = dyson_config.get('ea', False)
    target_state = dyson_config.get('target_state', mol.config.get('tdhf', {}).get('target', 1))
    
    oqp.compute_dyson_ekt(mol, target_state, compute_ip, compute_ea)


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
