import oqp
import copy
def project_basis(mol):

    oqp.library.ints_1e(mol)
    basis = mol.data.get_basis()
    nbf = int(basis['nbf'])
    nbf2 = int(nbf * (nbf + 1)/2)
    mol.data["OQP::VEC_MO_A_tmp"] = [[0] * nbf for _ in range(nbf)]
    mol.data["OQP::VEC_MO_B_tmp"] = [[0] * nbf for _ in range(nbf)]
    mol.data["OQP::DM_A_tmp"] = [0] * nbf2
    mol.data["OQP::DM_B_tmp"] = [0] * nbf2

    oqp.proj_dm_newbas(mol)

    mol.data["OQP::VEC_MO_A"] = copy.deepcopy(mol.data["OQP::VEC_MO_A_tmp"])
    mol.data["OQP::VEC_MO_B"] = copy.deepcopy(mol.data["OQP::VEC_MO_B_tmp"])
    mol.data["OQP::DM_A"] = copy.deepcopy(mol.data["OQP::DM_A_tmp"])
    mol.data["OQP::DM_B"] = copy.deepcopy(mol.data["OQP::DM_B_tmp"])
    mol.data["OQP::E_MO_A"] = [0] * nbf
    mol.data["OQP::E_MO_B"] = [0] * nbf
    mol.data["OQP::FOCK_A"] = [0] * nbf2
    mol.data["OQP::FOCK_B"] = [0] * nbf2
    mol.data["E_MO_A"] = copy.deepcopy(mol.data["OQP::E_MO_A_tmp"])
    mol.data["E_MO_B"] = copy.deepcopy(mol.data["OQP::E_MO_B_tmp"])
