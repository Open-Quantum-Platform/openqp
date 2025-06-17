"""Wrapper for OQP data"""
import os.path
from pathlib import Path
import numpy as np
from oqp import ffi, lib
from oqp.periodic_table import MASSES, SYMBOL_MAP
from oqp.utils.constants import ANGSTROM_TO_BOHR


def sarray(strng):
    """Convert array of parameters to list of strngs"""
    return list(s.strip().lower() for s in strng.split(',')) if strng else ()


def farray(strng):
    """Convert array of parameters to list of floats"""
    return list(float(s) for s in strng.split(',')) if strng else ()


def iarray(strng):
    """Convert array of parameters to list of integers"""
    return list(int(s) for s in strng.split(',')) if strng else ()


def barray(strng):
    """Convert array of parameters to list of booleans"""
    return list(bool(s) for s in strng.split(',')) if strng else ()


def parray(strng):
    """Convert array of parameters to list of integer pairs"""
    """e.g. 1 2, 3 4 -> [[1, 2], [3, 4]]"""
    return list([int(s.split()[0]), int(s.split()[1])] for s in strng.split(',')) if strng else ()


def string(strng):
    """Handle string parameters"""
    return strng.lower()


def path(strng):
    """Convert string to Path"""
    return Path(strng)


OQP_CONFIG_SCHEMA = {
    'input': {
        'charge': {'type': int, 'default': '0'},
        'basis': {'type': string, 'default': ''},
        'library': {'type': string, 'default': ''},
        'functional': {'type': string, 'default': ''},
        'method': {'type': string, 'default': 'hf'},
        'runtype': {'type': string, 'default': 'energy'},
        'system': {'type': str, 'default': ''},
        'system2': {'type': str, 'default': ''},
        'd4': {'type': bool, 'default': 'False'},
    },
    'guess': {
        'type': {'type': string, 'default': 'huckel'},
        'file': {'type': str, 'default': ''},
        'file2': {'type': str, 'default': ''},
        'save_mol': {'type': bool, 'default': 'False'},
        'continue_geom': {'type': bool, 'default': 'False'},
        'swapmo': {'type' : string, 'default' : ''},
    },
    'scf': {
        'type': {'type': string, 'default': 'rhf'},
        'maxit': {'type': int, 'default': '30'},
        'forced_attempt': {'type': int, 'default': '1'},
        'maxdiis': {'type': int, 'default': '7'},
        'diis_reset_mod': {'type': int, 'default': '10'},
        'diis_reset_conv': {'type': float, 'default': '0.005'},
        'diis_method_threshold': {'type': float, 'default': '2.0'},
        'diis_type': {'type': string, 'default': 'cdiis'},
        'vdiis_cdiis_switch': {'type': float, 'default': '0.3'},
        'vdiis_vshift_switch': {'type': float, 'default': '0.003'},
        'vshift_cdiis_switch': {'type': float, 'default': '0.3'},
        'vshift': {'type': float, 'default': '0.0'},
        'mom': {'type': bool, 'default': 'False'},
        'mom_switch': {'type': float, 'default': '0.003'},
        'pfon': {'type': bool, 'default': 'False'},
        'pfon_start_temp': {'type': float, 'default': '2000.0'},
        'pfon_cooling_rate': {'type': float, 'default': '50.0'},
        'pfon_nsmear': {'type': float, 'default': '5.0'},
        'multiplicity': {'type': int, 'default': '1'},
        'conv': {'type': float, 'default': '1.0e-6'},
        'incremental': {'type': bool, 'default': 'True'},
        'init_scf': {'type': string, 'default': 'no'},
        'init_basis': {'type': string, 'default': 'none'},
        'init_library': {'type': string, 'default': ''},
        'init_it': {'type': int, 'default': '15'},
        'init_conv': {'type': float, 'default': '0.001'}, 
        'init_converger': {'type': int, 'default': '0'},
        'save_molden': {'type': bool, 'default': 'True'},
        'rstctmo': {'type': bool, 'default': 'False'},
        'soscf_type': {'type': int, 'default': '0'},
        'soscf_reset_mod': {'type': int, 'default': '0'},
        'soscf_lvl_shift': {'type': float, 'default': '0'},
        'alternative_scf': {'type': bool, 'default': 'False'},
        'verbose': {'type': int, 'default': '1'},
    },
    'dftgrid': {
        'hfscale': {'type': float, 'default': '-1.0'},
        'cam_flag': {'type': bool, 'default': 'False'},
        'cam_alpha': {'type': float, 'default': '-1.0'},
        'cam_beta': {'type': float, 'default': '-1.0'},
        'cam_mu': {'type': float, 'default': '-1.0'},
        'rad_type': {'type': string, 'default': 'mhl'},
        'rad_npts': {'type': int, 'default': '50'},
        'ang_npts': {'type': int, 'default': '194'},
        'partfun': {'type': string, 'default': 'ssf'},
        'pruned': {'type': string, 'default': 'SG1'},
        'grid_ao_pruned': {'type': bool, 'default': 'True'},
        'grid_ao_threshold': {'type': float, 'default': '1.0e-15'},
        'grid_ao_sparsity_ratio': {'type': float, 'default': '0.9'},
    },
    'tdhf': {
        'type': {'type': string, 'default': 'rpa'},
        'maxit': {'type': int, 'default': '50'},
        'maxit_zv': {'type': int, 'default': '50'},
        'multiplicity': {'type': int, 'default': '1'},
        'conv': {'type': float, 'default': '1.0e-6'},
        'nstate': {'type': int, 'default': '1'},
        'zvconv': {'type': float, 'default': '1.0e-6'},
        'nvdav': {'type': int, 'default': '50'},
        'tlf': {'type': int, 'default': '2'},
        'hfscale': {'type': float, 'default': '-1.0'},
        'cam_alpha': {'type': float, 'default': '-1.0'},
        'cam_beta': {'type': float, 'default': '-1.0'},
        'cam_mu': {'type': float, 'default': '-1.0'},
        'spc_coco': {'type': float, 'default': '-1.0'},
        'spc_ovov': {'type': float, 'default': '-1.0'},
        'spc_coov': {'type': float, 'default': '-1.0'},
        'conf_threshold': {'type': float, 'default': '5.0e-2'},
        'ixcore': {'type' : string, 'default' : '-1'},
    },
    'properties': {
        'scf_prop': {'type': sarray, 'default': 'el_mom,mulliken'},
        'td_prop': {'type': bool, 'default': 'False'},
        'grad': {'type': iarray, 'default': '0'},
        'nac': {'type': str, 'default': ''},
        'soc': {'type': str, 'default': ''},
        'export': {'type': bool, 'default': 'False'},
        'title': {'type': str, 'default': ''},
        'back_door': {'type': bool, 'default': False}
    },
    'optimize': {
        'lib': {'type': str, 'default': 'scipy'},
        'optimizer': {'type': str, 'default': 'bfgs'},
        'step_size': {'type': float, 'default': '0.1'},
        'step_tol': {'type': float, 'default': '1e-2'},
        'maxit': {'type': int, 'default': 30},
        'mep_maxit': {'type': int, 'default': 10},
        'rmsd_grad': {'type': float, 'default': '1e-4'},
        'rmsd_step': {'type': float, 'default': '1e-3'},
        'max_grad': {'type': float, 'default': '3e-4'},
        'max_step': {'type': float, 'default': '2e-3'},
        'istate': {'type': int, 'default': '1'},
        'jstate': {'type': int, 'default': '2'},
        'kstate': {'type': int, 'default': '3'},
        'imult': {'type': int, 'default': '1'},
        'jmult': {'type': int, 'default': '3'},
        'energy_shift': {'type': float, 'default': '1e-6'},
        'energy_gap': {'type': float, 'default': '1e-5'},
        'meci_search': {'type': str, 'default': 'penalty'},
        'pen_sigma': {'type': float, 'default': '1.0'},
        'pen_alpha': {'type': float, 'default': '0.0'},
        'pen_incre': {'type': float, 'default': '1.0'},
        'gap_weight': {'type': float, 'default': '1.0'},
        'init_scf': {'type': bool, 'default': 'False'},
    },
    'dlfind': {
        'printl': {'type': int, 'default': '2'},
        'icoord': {'type': int, 'default': '3'},
        'iopt': {'type': int, 'default': '3'},
        'ims': {'type': int, 'default': '0'},
    },
    'hess': {
        'type': {'type': string, 'default': 'numerical'},
        'state': {'type': int, 'default': '0'},
        'dx': {'type': float, 'default': '0.01'},
        'nproc': {'type': int, 'default': '1'},
        'read': {'type': bool, 'default': 'False'},
        'restart': {'type': bool, 'default': 'False'},
        'temperature': {'type': farray, 'default': '298.15'},
        'clean': {'type': bool, 'default': 'False'},
    },
    'nac': {
        'type': {'type': string, 'default': 'numerical'},
        'dt': {'type': float, 'default': '1'},
        'dx': {'type': float, 'default': '0.0001'},
        'bp': {'type': bool, 'default': 'False'},
        'nproc': {'type': int, 'default': '1'},
        'restart': {'type': bool, 'default': 'False'},
        'clean': {'type': bool, 'default': 'False'},
        'states': {'type': parray, 'default': '1 2'},
        'align': {'type': str, 'default': 'reorder'},

    },
    'json':{
            'scf_type': {'type': string, 'default': ''},
            'basis': {'type': string, 'default': ''},
            'library': {'type': string, 'default': ''},
            'do_init': {'type': string, 'default': 'no'},
            },
    'tests': {
        'exception': {'type': bool, 'default': False},
    },
}

TA_DIMENSIONS_LENGTH = 12


class OQPData:
    """Wrapper for OQP data class"""

    _scftypes = {"rhf": 1, "uhf": 2, "rohf": 3}
    _guesses = {"huckel": 1, "hcore": 2}
    _dft_switch = {False: 10, True: 20}
    _methods = ('hf', 'tdhf')
    _td_types = ('rpa', 'tda', 'sf', 'mrsf')
    _rad_grid_types = {'mhl': 0, 'log3': 1, 'ta': 2, 'becke': 3}
    _diis_types = {'none': 1, 'cdiis': 2, 'ediis': 3, 'adiis': 4, 'vdiis': 5}
    _dftgrid_partition_functions = {'ssf': 0, 'becke': 1, 'erf': 2,
                                    'sstep2': 3, 'sstep3': 4, 'sstep4': 5, 'sstep5': 6}
    _handlers = {
        "input": {
            "charge": "set_mol_charge",
            "functional": "set_dft_functional",
            "system": "set_system",
            "system2": "set_system2",
        },
        "guess": {
        },
        "scf": {
            "type": "set_scf_type",
            "maxit": "set_scf_maxit",
            "maxdiis": "set_scf_maxdiis",
            "diis_reset_mod": "set_scf_diis_reset_mod",
            "diis_reset_conv": "set_scf_diis_reset_conv",
            "diis_method_threshold": "set_scf_diis_method_threshold",
            "diis_type": "set_scf_diis_type",
            "vdiis_cdiis_switch": "set_scf_vdiis_cdiis_switch",
            "vdiis_vshift_switch": "set_scf_vdiis_vshift_switch",
            "vshift_cdiis_switch": "set_scf_vshift_cdiis_switch",
            "ft": "set_scf_vshift",
            "vshift": "set_scf_vshift",
            "mom": "set_scf_mom",
            "mom_switch": "set_scf_mom_switch",
            "pfon": "set_scf_pfon",
            "pfon_start_temp": "set_scf_pfon_start_temp",
            "pfon_cooling_rate": "set_scf_pfon_cooling_rate",
            "pfon_nsmear": "set_scf_pfon_nsmear",
            "multiplicity": "set_mol_multiplicity",
            "conv": "set_scf_conv",
            "incremental": "set_scf_incremental",
            "active_basis": "set_scf_active_basis",
            "rstctmo" : "set_scf_rstctmo",
            "soscf_type": "set_scf_soscf_type",
            "soscf_reset_mod": "set_scf_soscf_reset_mod",
            "soscf_lvl_shift": "set_soscf_lvl_shift",
            "verbose": "set_scf_verbose",
        },
        "dftgrid": {
            "rad_type": "set_dftgrid_rad_type",
            "rad_npts": "set_dftgrid_rad_npts",
            "ang_npts": "set_dftgrid_ang_npts",
            "partfun": "set_dftgrid_partfun",
            "pruned": "set_dftgrid_pruned",
            "grid_ao_pruned": "set_dftgrid_ao_pruned",
            "grid_ao_threshold": "set_dftgrid_ao_threshold",
            "grid_ao_sparsity_ratio": "set_dftgrid_pruned_ao_sparsity_ratio",
            "hfscale": "set_dftgrid_hfscale",
            "cam_flag": "set_dftgrid_cam_flag",
            "cam_alpha": "set_dftgrid_cam_alpha",
            "cam_beta": "set_dftgrid_cam_beta",
            "cam_mu": "set_dftgrid_cam_mu",
        },
        "tdhf": {
            "type": "set_tdhf_type",
            "nstate": "set_tdhf_nstate",
            "multiplicity": "set_tdhf_multiplicity",
            "maxit": "set_tdhf_maxit",
            "maxit_zv": "set_tdhf_maxit_zv",
            "conv": "set_tdhf_conv",
            "zvconv": "set_tdhf_zvconv",
            "tlf": "set_tdhf_tlf",
            "hfscale": "set_tdhf_hfscale",
            "cam_alpha": "set_tdhf_cam_alpha",
            "cam_beta": "set_tdhf_cam_beta",
            "cam_mu": "set_tdhf_cam_mu",
            "spc_coco": "set_tdhf_spc_coco",
            "spc_ovov": "set_tdhf_spc_ovov",
            "spc_coov": "set_tdhf_spc_coov",
            "conf_threshold": "set_conf_threshold",
#            "ixcore" : "set_tdhf_ixcore",
        },
    }
    _typemap = [np.void,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
                np.float32,
                np.float64,
                np.complex64,
                np.complex128,
                np.dtype('S1'),
                ]

    @property
    def typemap(self):
        """Access type map"""
        return OQPData._typemap

    def __init__(self, silent=0):
        self._data = ffi.gc(lib.oqp_init(), lib.oqp_clean)
        self.silent = silent
        self.mol2 = []  # coordinates of the second molecule in 1D in Bohr

    def __getitem__(self, key):
        """Get data from molecule"""

        if key in dir(self._data.mol_prop):
            return getattr(self._data.mol_prop, key)

        if key in dir(self._data.mol_energy):
            return getattr(self._data.mol_energy, key)

        if key in dir(self._data.tddft):
            return getattr(self._data.tddft, key)
        if key in dir(self._data.mpiinfo):
            return getattr(self._data.mpiinfo, key)
        if key in dir(self._data.control):
            return getattr(self._data.control, key)
        if key in dir(self._data.elshell):
            return getattr(self._data.elshell, key)
        if key in dir(self._data):
            return getattr(self._data, key)

        code = bytes(key, 'ascii')
        req = ffi.new('char []', code)
        type_id = ffi.new('int32_t *')
        ndims = ffi.new('int32_t *')
        dims = ffi.new(f'int64_t[{TA_DIMENSIONS_LENGTH}]')
        data_ptr = ffi.new('void **')
        data_size = lib.oqp_get(self._data, req, type_id, ndims, dims, data_ptr)
        if data_size >= 0:
            elem_size = np.dtype(self.typemap[type_id[0]]).itemsize
            shape = np.frombuffer(ffi.buffer(dims, ffi.sizeof('int64_t') * ndims[0]), dtype=np.int64)
            data_type = self.typemap[type_id[0]]
            if np.prod(shape) == 1:
                val = np.frombuffer(ffi.buffer(data_ptr[0], elem_size * data_size),
                                    dtype=data_type).reshape(shape)
                if data_type == np.dtype('S1'):
                    return val.tobytes().decode()

            val = np.frombuffer(
                ffi.buffer(data_ptr[0],
                           elem_size * data_size),
                dtype=data_type,
            ).reshape(shape)

            if data_type == np.dtype('S1'):
                return val.tobytes().decode()
            return val

        raise AttributeError(f"Key `{key}` not found in QOP data")

    def __setitem__(self, key, value):
        """Set data in molecule"""

        if key in dir(self._data.mol_prop):
            self._data.mol_prop[key] = value

        if key in dir(self._data.mol_energy):
            self._data.mol_energy[key] = value

        if key in dir(self._data.mpiinfo):
            setattr(self._data.mpiinfo, key, value)

        if key in dir(self._data.tddft):
            setattr(self._data.tddft, key, value)

        if key in dir(self._data.elshell):
            setattr(self._data.elshell, key, value)
            return

        if isinstance(value, np.ndarray):
            _value = value
        elif isinstance(value, str):
            _value = np.frombuffer(np.bytes_(value), dtype=np.dtype('S1'))
        elif isinstance(value, ffi.CData):
            try:
                _value = np.frombuffer(ffi.buffer(value), dtype=np.int32)
            except Exception as e:
                raise TypeError("CData pointer is not buffer-backed or dtype mismatch") from e            
        else:
            _value = np.array(value)

        try:
            typeid = self.typemap.index(_value.dtype)
            self._setitem_internal(key, _value, typeid)
        except:
            print(f"The type of your data is not allowed: {_value.dtype.name}")

    def _setitem_internal(self, key, value, typeid):
        code = bytes(key, 'ascii')
        req = ffi.new('char []', code)
        type_id = ffi.new('int32_t *', typeid)
        shape = np.shape(value)
        dims = ffi.new(f'int64_t[{TA_DIMENSIONS_LENGTH}]', shape)
        ndims = ffi.new('int32_t *', len(np.shape(value)))
        data_ptr = ffi.new('void **')
        data_size = lib.oqp_alloc(self._data, req, type_id, ndims, dims, data_ptr)
        if data_size < 0:
            raise AttributeError("Cannot allocate memory to store object in QOP data")

        elem_size = np.dtype(self.typemap[typeid]).itemsize
        ffi.memmove(data_ptr[0], value, data_size * elem_size)

    def __delitem__(self, key):
        """Erase data in molecule"""
        code = bytes(key, 'ascii')
        req = ffi.new('char []', code)
        result = lib.oqp_del(self._data, req)
        if result == -1:
            raise ValueError("Handle not initialized")
        if result == -2:
            raise KeyError(f"{key}")
        if result == -3:
            raise KeyError(f"Error when deleting entry {key} in OQP data")

    def set_mol_charge(self, charge):
        """Set charge of a molecule"""
        self._data.mol_prop.charge = charge

    def set_mol_multiplicity(self, multiplicity):
        """Set multiplicity of the molecule"""
        self._data.mol_prop.mult = multiplicity

    def set_scf_type(self, scftype):
        """Set SCF type"""
        self._data.control.scftype = OQPData._scftypes[scftype]

    def set_scf_maxit(self, maxit):
        """Set maximum number of SCF iterations"""
        self._data.control.maxit = maxit

    def set_scf_maxdiis(self, maxdiis):
        """Set maximum number of DIIS Equations"""
        self._data.control.maxdiis = maxdiis

    def set_scf_diis_reset_mod(self, diis_reset_mod):
        """Set reset DIIS Equations for every diis_reset_mod"""
        self._data.control.diis_reset_mod = diis_reset_mod

    def set_scf_diis_reset_conv(self, diis_reset_conv):
        """Set reset DIIS Equations for every diis_reset_mod"""
        self._data.control.diis_reset_conv = diis_reset_conv

    def set_scf_diis_method_threshold(self, diis_method_threshold):
        """Set DIIS threshold to switch DIIS method"""
        self._data.control.diis_method_threshold = diis_method_threshold

    def set_scf_diis_type(self, diistype):
        """Set DIIS method"""
        self._data.control.diis_type = OQPData._diis_types[diistype]

    def set_scf_vdiis_cdiis_switch(self, vdiis_cdiis_switch):
        """Set vdiis_cdiis_switch size for better SCF convergence"""
        self._data.control.vdiis_cdiis_switch = vdiis_cdiis_switch

    def set_scf_vdiis_vshift_switch(self, vdiis_vshift_switch):
        """Set vdiis_vshift_switch size for better SCF convergence"""
        self._data.control.vdiis_vshift_switch = vdiis_vshift_switch

    def set_scf_vshift_cdiis_switch(self, vshift_cdiis_switch):
        """Set vshift_cdiis_switch size for better SCF convergence"""
        self._data.control.vshift_cdiis_switch = vshift_cdiis_switch

    def set_scf_vshift(self, vshift):
        """Set Vshift size for better SCF convergency"""
        self._data.control.vshift = vshift

    def set_scf_mom(self, mom):
        """Set MOM for better SCF convergency"""
        self._data.control.mom = mom

    def set_scf_mom_switch(self, mom_switch):
        """Set MOM turn on criteria of DIIS error """
        self._data.control.mom_switch = mom_switch

    def set_scf_pfon(self, pfon):
        """pfon """
        self._data.control.pfon = pfon

    def set_scf_pfon_start_temp(self, pfon_start_temp):
        """pfon_start_temp """
        self._data.control.pfon_start_temp = pfon_start_temp

    def set_scf_pfon_cooling_rate(self, pfon_cooling_rate):
        """pfon_cooling_rate """
        self._data.control.pfon_cooling_rate = pfon_cooling_rate

    def set_scf_pfon_nsmear(self, pfon_nsmear):
        """pfon_cooling_rate """
        self._data.control.pfon_nsmear = pfon_nsmear

    def set_scf_rstctmo(self, rstctmo): 
        """restrict MO """
        self._data.control.rstctmo = rstctmo

    def set_scf_active_basis(self, active_basis):
        """Select basis set: 0 => info%basis
                             1 => info%alt_basis"""
        self._data.control.active_basis = active_basis

    def set_scf_conv(self, conv):
        """Set SCF convergence threshold"""
        self._data.control.conv = conv

    def set_scf_incremental(self, flag):
        """Set incremental Fock matrix build"""
        self._data.control.scf_incremental = 1 if flag else 0

    def set_scf_soscf_type(self, soscf_type):
        """Set SOSCF type for SCF convergence:
            soscf_type (int): SOSCF algorithm type
                0: SOSCF disabled
                1: SOSCF only
                2: SOSCF+DIIS combined mode
        """
        self._data.control.soscf_type = soscf_type

    def set_soscf_lvl_shift(self, soscf_lvl_shift):
        """Reset the orbital Hessian. If it is zero, we don't reset by default.
        """
        self._data.control.soscf_lvl_shift = soscf_lvl_shift

    def set_scf_soscf_reset_mod(self, soscf_reset_mod):
        """Set the SOSCF Hessian reset mode.
        Parameters:
            soscf_reset_mod (int):
                0      – Disable Hessian reset.
                >0     – Reset the Hessian at the specified SCF iteration.
        """
        self._data.control.soscf_reset_mod = soscf_reset_mod

    def set_scf_verbose(self, verbose):
        """Controls output verbosity"""
        self._data.control.verbose = verbose

    def set_tdhf_type(self, td_type):
        """Handle td-dft calculation type"""
        if td_type.lower() == 'tda':
            self._data.tddft.tda = True

    def set_tdhf_nstate(self, nstate):
        """Set number of states in tdhf calculation"""
        self._data.tddft.nstate = nstate

    def set_tdhf_target(self, target):
        """Set target states in tdhf gradient calculation"""
        self._data.tddft.target_state = target

    def set_tdhf_maxit(self, maxit):
        """Set max number of iterations in Davidson's eigensolver"""
        self._data.control.maxit_dav = maxit

    def set_tdhf_maxit_zv(self, maxit_zv):
        """Set max number of iterations in Davidson's eigensolver"""
        self._data.control.maxit_zv = maxit_zv

    def set_tdhf_conv(self, conv):
        """Set SCF convergence threshold"""
        self._data.tddft.cnvtol = conv

    def set_tdhf_zvconv(self, conv):
        """Set SCF convergence threshold"""
        self._data.tddft.zvconv = conv

    def set_tdhf_multiplicity(self, multiplicity):
        """Set multiplicity in tdhf calculation"""
        self._data.tddft.mult = multiplicity

    def set_tdhf_nvdav(self, nvdav):
        """Set max number of trial vectors in Davidson's eigensolver"""
        self._data.tddft.maxvec = nvdav

    def set_tdhf_tlf(self, tlf):
        """Set TLF in tdhf NAC calculation"""
        self._data.tddft.tlf = tlf

    def set_tdhf_hfscale(self, hfscale):
        """Set HF exact exchange scalar in response calculation"""
        self._data.tddft.hfscale = hfscale

    def set_tdhf_cam_alpha(self, cam_alpha):
        """Set short range HF exact exchange scalar in response calculation"""
        self._data.tddft.cam_alpha = cam_alpha

    def set_tdhf_cam_beta(self, cam_beta):
        """Set long range HF exact exchange scalar in response calculation"""
        self._data.tddft.cam_beta = cam_beta

    def set_tdhf_cam_mu(self, cam_mu):
        """Set range separation parameter mu in response calculation"""
        self._data.tddft.cam_mu = cam_mu

    def set_tdhf_spc_coco(self, spc_coco):
        """Set CO-CO spin-pair coupling parameter (C=closed, O=open, V=virtual MOs) in MRSF calculation"""
        self._data.tddft.spc_coco = spc_coco

    def set_tdhf_spc_ovov(self, spc_ovov):
        """Set OV-OV spin-pair coupling parameter (C=closed, O=open, V=virtual MOs) in MRSF calculation"""
        self._data.tddft.spc_ovov = spc_ovov

    def set_tdhf_spc_coov(self, spc_coov):
        """Set CO-OV spin-pair coupling parameter (C=closed, O=open, V=virtual MOs) in MRSF calculation"""
        self._data.tddft.spc_coov = spc_coov

    def set_conf_threshold(self, conf_threshold):
        """Set configuration printout option"""
        self._data.control.conf_print_threshold = conf_threshold

    def set_dft_functional(self, functional):
        """Set DFT functional"""
        dft = functional != ''
        if dft:

            if self.silent != 1:
                print(f'functional={functional}')

            self._data.dft.XC_functional_name = (
                functional.ljust(20)[:20].upper().encode("ascii")
            )
        self._data.control.hamilton = OQPData._dft_switch[dft]

    def set_dftgrid_rad_type(self, radtype):
        """Set radial grid type in DFT"""
        self._data.dft.rad_grid_type = OQPData._rad_grid_types[radtype]

    def set_dftgrid_rad_npts(self, npts):
        """Set number of radial grid points in DFT"""
        self._data.dft.grid_rad_size = npts

    def set_dftgrid_ang_npts(self, npts):
        """Set number of angular grid points in DFT"""
        self._data.dft.grid_ang_size = npts

    def set_dftgrid_ao_threshold(self, do):
        """Set grid_ao_threshold"""
        self._data.dft.grid_ao_threshold = do

    def set_dftgrid_ao_pruned(self, do):
        """Set grid_ao_pruned"""
        self._data.dft.grid_ao_pruned = do

    def set_dftgrid_pruned_ao_sparsity_ratio(self, do):
        """Set grid_ao_sparsity_ratio """
        self._data.dft.grid_ao_sparsity_ratio = do

    def set_dftgrid_partfun(self, partfun):
        """Set partition function in Becke's fuzzy cell method"""
        self._data.dft.dft_partfun = OQPData._dftgrid_partition_functions[partfun]

    def set_dftgrid_pruned(self, pruned):
        """Set pruned grid"""
        pruned_list = ['SG1', ]
        if pruned != "":
            pruned = pruned.upper()
            if pruned in pruned_list:
                self._data.dft.grid_pruned = True
                self._data.dft.grid_pruned_name = pruned.ljust(16)[:16].upper().encode("ascii")
            else:
                print(f"{pruned} grid is not valid. Available options are: {', '.join(pruned_list)}")

    def set_dftgrid_hfscale(self, hfscale):
        """Set HF exact exchange scalar in DFT calculation"""
        self._data.dft.hfscale = hfscale

    def set_dftgrid_cam_flag(self, cam_flag):
        """Set CAM flag in DFT calculation"""
        self._data.dft.cam_flag = cam_flag

    def set_dftgrid_cam_alpha(self, cam_alpha):
        """Set short range HF exact exchange scalar in DFT calculation"""
        self._data.dft.cam_alpha = cam_alpha

    def set_dftgrid_cam_beta(self, cam_beta):
        """Set long range HF exact exchange scalar in DFT calculation"""
        self._data.dft.cam_beta = cam_beta

    def set_dftgrid_cam_mu(self, cam_mu):
        """Set range separation parameter mu in DFT calculation"""
        self._data.dft.cam_mu = cam_mu

    def set_system(self, system):
        """Set up atomic data"""
        num_atoms, x, y, z, q, mass = read_system(system)
        self._data.mol_prop.natom = num_atoms
        lib.oqp_set_atoms(self._data, num_atoms, x, y, z, q, mass)

    def set_system2(self, system):
        """Set up the second set of atomic data"""
        if system.strip():
            num_atoms, x, y, z, q, mass = read_system(system)
            self.mol2 = np.array(x + y + z).reshape((3, num_atoms)).T.reshape(-1)

    def parse_section(self, config, section):
        cfg_input = config[section]
        for key in cfg_input.keys():
            val = cfg_input[key]
            try:
                handler = getattr(self, OQPData._handlers[section][key])
                handler(val)
            except KeyError:
                continue

    def apply_config(self, config):
        """
        Apply the data from the OQP config
        The latter has to be read from the input file
        """
        for section in config:
            self.parse_section(config, section)

        molecule = self._data
        natom = molecule.mol_prop.natom
        charge = molecule.mol_prop.charge
        nelec = sum(int(molecule.qn[i]) for i in range(natom)) - charge

        molecule.mol_prop.nelec = nelec
        na, nb = compute_alpha_beta_electrons(nelec, molecule.mol_prop.mult)
        molecule.mol_prop.nelec_A = na
        molecule.mol_prop.nelec_B = nb
        molecule.mol_prop.nocc = max(na, nb)
        if molecule.control.scftype == 3 and molecule.mol_prop.mult == 1:
            print("WARNING! ROHF + multiplicity = 1 has bugs!")
            print("Do not trust to these results!")

    def get_basis(self):
        """Get basis set from a molecule"""
        pex = ffi.new('double **')
        pcc = ffi.new('double **')
        pdeg = ffi.new('int64_t **')
        pat = ffi.new('int64_t **')
        pam = ffi.new('int64_t **')
        pnbf = ffi.new('int64_t *')
        pnsh = ffi.new('int64_t *')
        pnprim = ffi.new('int64_t *')

        ret = lib.oqp_get_basis(self._data, pnsh, pnprim, pnbf,
                                pam, pat, pdeg, pex, pcc)

        basis = {}

        if ret == 0:
            nbf = pnbf[0]
            nsh = pnsh[0]
            nprim = pnprim[0]

            centers = np.frombuffer(ffi.buffer(pat[0], ffi.sizeof('int64_t') * nsh), dtype=np.int64)
            angs = np.frombuffer(ffi.buffer(pam[0], ffi.sizeof('int64_t') * nsh), dtype=np.int64)
            ncontr = np.frombuffer(ffi.buffer(pdeg[0], ffi.sizeof('int64_t') * nsh), dtype=np.int64)

            alpha = np.frombuffer(ffi.buffer(pex[0], ffi.sizeof('double') * nprim))
            coef = np.frombuffer(ffi.buffer(pcc[0], ffi.sizeof('double') * nprim))

            basis = {
                'centers': np.copy(centers) - 1,  # make zero-based indexing of atoms
                'angs': np.copy(angs),
                'ncontr': np.copy(ncontr),
                'alpha': np.copy(alpha),
                'coef': np.copy(coef),
                'nbf': np.copy(nbf),
                'nsh': np.copy(nsh),
                'nprim': np.copy(nprim),
            }

        return basis


def compute_alpha_beta_electrons(n_e, mult):
    """
    Compute number of alpha and beta electrons for a given total electron number and multiplicity
    ne - total number of electrons
    mult - multiplicity
    returns (n_alpha, n_beta)
    """
    n_a = n_e + (abs(mult) - 1)
    n_b = n_e - (abs(mult) - 1)
    if n_a % 2 != 0 or n_b % 2 != 0 or mult == 0 or n_a + n_b != 2 * n_e:
        raise ValueError(f"Impossible multiplicity and number of electrons combination: ne={n_e}, mult={mult}")

    n_a //= 2
    n_b //= 2

    return (n_a, n_b) if mult > 0 else (n_b, n_a)


def read_system(system):
    system = system.split("\n")
    if system[0]:
        if not os.path.exists(system[0]):
            raise FileNotFoundError("XYZ file %s is not found!" % system[0])

        with open(system[0], 'r') as xyzfile:
            system = xyzfile.read().splitlines()

        num_atoms = int(system[0])
        system = system[2: 2 + num_atoms]
    else:
        system = system[1:]
        num_atoms = len(system)

    atoms = []
    for i, line in enumerate(system):
        line = line.split()
        if len(line) >= 4:
            atoms.append(line[0: 4])
        else:
            print(f"{system[i]} is not valid line for atom configuration!")

    q = [float(SYMBOL_MAP[atoms[i][0]]) for i in range(0, num_atoms)]
    x = [float(atoms[i][1]) / ANGSTROM_TO_BOHR for i in range(0, num_atoms)]
    y = [float(atoms[i][2]) / ANGSTROM_TO_BOHR for i in range(0, num_atoms)]
    z = [float(atoms[i][3]) / ANGSTROM_TO_BOHR for i in range(0, num_atoms)]
    mass = [MASSES[int(SYMBOL_MAP[atoms[i][0]])] for i in range(0, num_atoms)]

    return num_atoms, x, y, z, q, mass

