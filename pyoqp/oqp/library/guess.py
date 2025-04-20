"""Set up initial guess density"""

import os
import copy
import oqp
from oqp.utils.file_utils import try_basis
from oqp.utils.file_utils import dump_log


def update_guess(mol):
    if mol.config['json']['scf_type'] == 'rhf':
        mol.data["OQP::VEC_MO_B"] = copy.deepcopy(mol.data["OQP::VEC_MO_A"])
        mol.data["OQP::E_MO_B"] = copy.deepcopy(mol.data["OQP::E_MO_A"])
        mol.data["OQP::DM_B"] = copy.deepcopy(mol.data["OQP::DM_A"])
    oqp.guess_json(mol)

def guess(mol):
    """Set up initial guess density"""

    guess_type = mol.config["guess"]["type"]
    guess_file = 'compute orbitals'
    swapmo= mol.config["guess"]["swapmo"]

    if guess_type == "huckel":
        hubas = try_basis("MINI_huckel", fallback=None)
        mol.data["OQP::hbasis_filename"] = hubas
        oqp.guess_huckel(mol)
        alpha = 'computed'
        beta = 'computed'

    elif guess_type == "hcore":
        oqp.guess_hcore(mol)
        alpha = 'computed'
        beta = 'computed'

    elif guess_type == 'json':
        guess_file = mol.config["guess"]["file"]
        if mol.config['scf']['type'] != 'rhf':
            update_guess(mol)
        alpha = 'reloaded'
        beta = 'reloaded'

    elif guess_type == 'auto':
        guess_file = mol.config["guess"]["file"]
        if os.path.exists(guess_file):
            alpha = 'reloaded'
            beta = 'reloaded'
        else:
            hubas = try_basis("MINI_huckel", fallback=None)
            mol.data["OQP::hbasis_filename"] = hubas
            oqp.guess_huckel(mol)
            alpha = 'computed'
            beta = 'computed'

#    # molden does not have sufficient numerical accuracy
#    elif guess_type == "molden":
#        # Check if the molden file from the input exists
#        if os.path.isfile(mol.config["guess"]["file"]):
#            guess_file = mol.config["guess"]["file"]
#        # Check if the default name with '.molden' exists
#        elif os.path.isfile(mol.input_file.replace('.inp', '.molden')):
#            guess_file = mol.input_file.replace('.inp', '.molden')
#        # If neither exists, raise an error
#        else:
#            raise FileNotFoundError(f'Molden file not found.')
#
#        reader = MoldenReader(guess_file)
#        mo_data = reader.read_mo()
#        mol.data["OQP::VEC_MO_A"] = mo_data.get("mo_vec_a", None)
#        mol.data["OQP::E_MO_A"] = mo_data.get("mo_e_a", None)
#        mol.data["OQP::DM_A"] = np.array(mo_data.get("dens_a", None))
#        alpha = 'read'
#        beta = 'read'
#        if mo_data.get("mo_vec_b", None) is not None:
#            # raise ValueError(f"Beta orbitals are missing in the molden file '{guess_file}'")
#            # copy alpha to beta
#            mol.data["OQP::VEC_MO_B"] = mo_data.get("mo_vec_b", None)
#            mol.data["OQP::E_MO_B"] = mo_data.get("mo_e_b", None)
#            mol.data["OQP::DM_B"] = np.array(mo_data.get("dens_b", None))

    else:
        raise ValueError(f'Unknown guess type={guess_type}')

    try:
        mol.data["OQP::VEC_MO_B"]

    except AttributeError:
        mol.data["OQP::VEC_MO_B"] = copy.deepcopy(mol.data["OQP::VEC_MO_A"])
        mol.data["OQP::E_MO_B"] = copy.deepcopy(mol.data["OQP::E_MO_A"])
        mol.data["OQP::DM_B"] = copy.deepcopy(mol.data["OQP::DM_A"])
        beta = 'copied'

    guess_info = {
        'guess_type': guess_type,
        'guess_file': guess_file,
        'guess_alpha': alpha,
        'guess_beta': beta,
        'guess_swapmo': swapmo,
    }

    dump_log(mol, title='   PyOQP: Orbital Guess', section='guess', info=guess_info)
