"""OQP run type functions"""
import os

import numpy as np

from contextlib import contextmanager, nullcontext

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

@contextmanager
def _single_gradient_target(mol, target_types=("umrsf",)):
    grad_list = mol.config.get("properties", {}).get("grad", [])
    method = mol.config.get("input", {}).get("method")
    td_type = mol.config.get("tdhf", {}).get("type")
    if isinstance(td_type, str):
        td_type = td_type.lower()

    target = None
    if method == "tdhf" and td_type in target_types and len(grad_list) == 1:
        target = int(grad_list[0])

    if target is not None and target > 0:
        # Select the gradient state before the energy solve. Davidson still
        # converges every reported root, since all are published downstream.
        mol.data.set_tdhf_target(target)
    yield


def _soc_requested(mol):
    """True when this run computes spin-orbit coupling as a step of the flow."""
    from oqp.utils.spin_states import soc_requested
    return soc_requested(mol.config)


def _spin_setup(mol):
    """This run's spin labels and the gradient request of each ladder.

    Returns (configured label, other label, {'s': [...], 't': [...]}).
    """
    from oqp.utils.spin_states import configured_spin, grad_requests, other_spin

    label = configured_spin(mol.config)
    return label, other_spin(label), grad_requests(mol.config)


def _configured_manifold(sp, label, per_manifold):
    """Enter the configured ladder when a two-ladder run sizes it.

    response_manifold sizes the ladder it answers from nstate_s / nstate_t.
    A run that solves one ladder, or names neither count, keeps [tdhf] nstate
    and answers whatever multiplicity it configured, so nothing is entered and
    nothing changes -- as before nstate_s / nstate_t meant anything outside
    runtype=soc.
    """
    from oqp.library.single_point import response_manifold

    td = sp.mol.config['tdhf']
    named = (per_manifold and int(td.get('multiplicity', 1)) in (1, 3)
             and int(td.get(f'nstate_{label}', 0)))
    return response_manifold(sp, label) if named else nullcontext()


def _soc_driver(mol, sp=None, ref_energy=None, tb_soc=None):
    """The SOC calculator for this flow.

    Imported at use time, like the native optimizers in get_optimizer(): several
    harnesses import runfunc with a stubbed single_point that carries only the
    names it used to need, and a module-level import would break them all.
    """
    from oqp.library.single_point import SOC
    return SOC(mol, sp=sp, ref_energy=ref_energy, tb_soc=tb_soc)


def _tb_soc(mol):
    """The tight-binding SOC driver for this backend.

    It solves both spin manifolds itself and publishes the same OQP:: tags as
    the native kernel, so there is no per-manifold round to add around it.
    """
    if tb_section_name(mol.config) == 'xtb':
        from oqp.library.openqp_xtb import xtb_soc as tb_soc
    else:
        from oqp.library.openqp_dftb import dftb_soc as tb_soc
    return tb_soc


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


def write_acid_cubes(mol):
    """Write the ACID scalar and the induced-current components as cubes.

    Runs after ``scf_prop=nmr`` with ``nmr_gauge=giao``, which is what publishes
    the magnetic density response the current density is built from.  The
    magnetic field direction only selects which vector field is written; the
    ACID scalar is a tensor invariant and does not depend on it.
    """
    from oqp.export import AcidExporter
    from oqp.utils.mpi_utils import MPIManager

    # Every rank reaches compute_scf_prop, and this writes four shared cube
    # paths.  Left unguarded they race on the same files and each rank rebuilds
    # the whole grid.  world_rank, not rank: inside a split communicator the
    # group roots would all pass.
    if getattr(mol, "usempi", False) and MPIManager().world_rank != 0:
        return []

    properties = mol.config.get("properties", {})
    spacing = float(properties.get("acid_spacing", 0.2))
    padding = float(properties.get("acid_padding", 5.0))
    # nan and inf survive float() and slip past an ordering test, and the grid
    # builder then produces a degenerate box of non-finite coordinates and
    # writes cubes nothing can read.  Reject them as input, not as output.
    if not (np.isfinite(spacing) and np.isfinite(padding)):
        raise ValueError(
            "properties.acid_spacing and acid_padding must be finite; got "
            f"spacing={spacing}, padding={padding}."
        )
    if spacing <= 0.0 or padding < 0.0:
        raise ValueError(
            "properties.acid_spacing must be positive and acid_padding "
            f"non-negative; got spacing={spacing}, padding={padding}."
        )
    try:
        exporter = AcidExporter(mol, padding=padding, spacing=spacing)
    except ValueError as error:
        raise ValueError(
            f"{error} Add 'nmr' before 'acid' in properties.scf_prop and set "
            "properties.nmr_gauge=giao."
        ) from None

    prefix = mol.log[:-4] if mol.log.endswith('.log') else mol.log
    paths = exporter.write_cubes(prefix, bfield=(0.0, 0.0, 1.0))
    with open(mol.log, 'a', encoding='utf-8') as fout:
        fout.write('\n    ACID cubes written: '
                   + ', '.join(os.path.basename(path) for path in paths) + '\n')
    return paths


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
        elif prop == 'acid':
            write_acid_cubes(mol)
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
    from oqp.library.single_point import grad_request, response_manifold

    label, other, requests = _spin_setup(mol)
    # A second ladder is published only when something named it, and only then do
    # the file names carry a spin label.
    per_manifold = bool(requests[other])
    spin_label = label if per_manifold else None

    sp = SinglePoint(mol)
    with _configured_manifold(sp, label, per_manifold):
        # compute energy
        # Set the requested response state while converging the full reported
        # spectrum.  grad_request first: _single_gradient_target reads
        # [properties] grad.
        with grad_request(mol, requests[label]), _single_gradient_target(mol):
            sp.energy()

            # compute gradient, unless only the other ladder asked for one
            if requests[label] or not per_manifold:
                Gradient(mol).gradient()

        # compute properties
        compute_scf_prop(mol)
#        if mol.config['input']['qmmm_flag']:
#           oqp.resp_charges(mol)

        # compute dftd4
        LastStep(mol).compute(mol, grad_list=requests[label], spin_label=spin_label)

    if per_manifold:
        # Record each ladder under its own label, as the properties and data
        # flows do, so a checkpoint of this run can seed a later step's history.
        record = _soc_driver(mol, sp=sp)
        record.snapshot(label)

        # compute the gradient of the other spin manifold
        with response_manifold(sp, other):
            sp.excitation([mol.mol_energy.energy])
            # before LastStep saves the molecule
            record.snapshot(other)
            with grad_request(mol, requests[other]):
                Gradient(mol).gradient()
            LastStep(mol).compute(mol, grad_list=requests[other], spin_label=other)

def compute_md(mol):
    """Run gas-phase ground-state Born--Oppenheimer dynamics.

    Ground-state QM/MM MD is dispatched separately to ``QMMM_MD`` before this
    function is reached.  Keeping the all-QM driver here gives ``runtype=md``
    one unambiguous meaning in the ordinary Runner dispatch.
    """
    from oqp.library.ground_state_md import GroundStateMD

    GroundStateMD(mol).run()

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
    from oqp.library.single_point import response_manifold

    if is_tb_method(mol.config['input']['method']):
        # The tight-binding driver solves both spin manifolds itself.
        _soc_driver(mol, tb_soc=_tb_soc(mol)).coupling()
        LastStep(mol).compute(mol)
        return

    sp = SinglePoint(mol)
    ref_energy = sp.reference()          # one SCF feeds both manifolds
    soc = _soc_driver(mol, sp=sp, ref_energy=ref_energy)

    # compute the two response manifolds the kernel reads, and publish each
    # under its own label.  Gradients belong to runtype=grad, prop and data.
    for label in ('s', 't'):
        with response_manifold(sp, label):
            sp.excitation(ref_energy)
            soc.snapshot(label)
            LastStep(mol).compute(mol, spin_label=label)

    # compute spin-orbit coupling
    soc.coupling()

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
    from oqp.library.single_point import grad_request, response_manifold

    soc_on = _soc_requested(mol)
    label, other, requests = _spin_setup(mol)
    grad_list = requests[label]
    # Two ladders are published side by side whenever a second one was asked
    # for, by its own selector or by spin-orbit coupling; only then does a file
    # name have to say which manifold it holds.
    per_manifold = soc_on or bool(requests[other])
    spin_label = label if per_manifold else None
    nac_type = mol.config['properties']['nac']

    # compute reference energy
    sp = SinglePoint(mol)
    ref_energy = sp.reference()

    with _configured_manifold(sp, label, per_manifold):
        # compute mo overlap and apply alignment
        BasisOverlap(mol).overlap()

        # compute excitation energy; _single_gradient_target reads this
        # ladder's request, so install it first
        with grad_request(mol, grad_list), _single_gradient_target(mol):
            sp.excitation(ref_energy)

        # Transport the central-geometry response-vector gauge before gradients,
        # numerical-NAC workers, and the checkpoint JSON are produced.  PyRAI2MD
        # feeds that JSON into the next MD step; aligning only inside displaced
        # workers leaves the time-series phase undefined.
        nacme = NACME(mol)
        if per_manifold:
            # The unlabelled history tag holds whichever ladder the previous step
            # answered last, so align this ladder against its own record first.
            nacme.load_state_history(label)
        nacme.align_x()

        # compute gradient, unless only the other ladder asked for one
        if grad_list or not per_manifold:
            with grad_request(mol, grad_list):
                Gradient(mol).gradient()

        # compute dftd4
        LastStep(mol).compute(mol, grad_list=grad_list, spin_label=spin_label)

        # compute nac or nacme
        if nac_type == 'nac':
            NAC(mol).nac(spin_label=spin_label)
        elif nac_type == 'nacme':
            nacme.nacme(align=False, spin_label=spin_label)

    if requests[other] or soc_on:
        other_grad_list = requests[other]
        soc = _soc_driver(mol, sp=sp, ref_energy=ref_energy)
        # Record each ladder under its own label whenever two are solved, not
        # only for the coupling: the unlabelled tags end up holding the ladder
        # solved last, so the next step pairs each ladder with its own history
        # through these records (NACME.load_state_history).
        soc.snapshot(label)

        # compute the other spin manifold
        with response_manifold(sp, other):
            sp.excitation(ref_energy)
            # before LastStep saves the molecule
            soc.snapshot(other)

            # compute its gradient, if it requested one
            if other_grad_list:
                with grad_request(mol, other_grad_list):
                    Gradient(mol).gradient()
                    LastStep(mol).compute(mol, grad_list=other_grad_list,
                                          spin_label=other)
            else:
                LastStep(mol).compute(mol, spin_label=other)

            # compute its couplings, against its own previous step
            if soc_on and nac_type == 'nacme':
                other_nacme = NACME(mol)
                if other_nacme.load_state_history(other):
                    other_nacme.align_x()
                    other_nacme.nacme(align=False, spin_label=other)

        if soc_on:
            # publish the phase-fixed amplitudes, then the coupling between them
            soc.snapshot(other)
            soc.coupling()

def compute_data(mol):
    from oqp.library.single_point import grad_request, response_manifold

    soc_on = _soc_requested(mol)
    label, other, requests = _spin_setup(mol)
    grad_list = requests[label]
    # Two ladders are published side by side whenever a second one was asked
    # for, by its own selector or by spin-orbit coupling; only then does a file
    # name have to say which manifold it holds.
    per_manifold = soc_on or bool(requests[other])
    spin_label = label if per_manifold else None

    # compute energy
    sp = SinglePoint(mol)

    with _configured_manifold(sp, label, per_manifold):
        # _single_gradient_target reads this ladder's request, so install it first
        with grad_request(mol, grad_list), _single_gradient_target(mol):
            sp.energy()

        # compute gradient, unless only the other ladder asked for one
        if grad_list or not per_manifold:
            with grad_request(mol, grad_list):
                Gradient(mol).gradient()

        # compute dftd4
        LastStep(mol).compute(mol, grad_list=grad_list, spin_label=spin_label)

        # compute nac
        if mol.config['properties']['nac'] == 'nac':
            NAC(mol).nac(spin_label=spin_label)

    if requests[other] or soc_on:
        other_grad_list = requests[other]
        soc = _soc_driver(mol, sp=sp, ref_energy=[mol.mol_energy.energy])
        # Record each ladder under its own label whenever two are solved, as
        # compute_properties does: a later step reads its history from them.
        soc.snapshot(label)

        # compute the other spin manifold
        with response_manifold(sp, other):
            sp.excitation(soc.ref_energy)
            # before LastStep saves the molecule
            soc.snapshot(other)

            # compute its gradient, if it requested one
            if other_grad_list:
                with grad_request(mol, other_grad_list):
                    Gradient(mol).gradient()
                    LastStep(mol).compute(mol, grad_list=other_grad_list,
                                          spin_label=other)
            else:
                LastStep(mol).compute(mol, spin_label=other)

        if soc_on:
            # compute the coupling between the two recorded ladders
            soc.coupling()


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
