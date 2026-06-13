import openmm.app as app
import openmm as mm
import openmm.unit as unit
import numpy as np

from oqp.openqp import OPENQP
from oqp.library.single_point import (
    SinglePoint, Gradient, Hessian, LastStep,
    BasisOverlap, NACME, NAC
)
from oqp.utils.file_utils import dump_log, dump_data, write_config, write_xyz

import oqp


def unpack_lower_tri_single(packed_atom, nbf):
    """
    Unpack a single lower-triangular packed array (length nbf*(nbf+1)/2)
    into a full symmetric (nbf x nbf) matrix.
    """
    packed_atom = np.asarray(packed_atom)
    nbf_tri = nbf * (nbf + 1) // 2
    if packed_atom.size != nbf_tri:
        raise ValueError(f"Size mismatch: got {packed_atom.size}, expected {nbf_tri}")
    full = np.zeros((nbf, nbf), dtype=packed_atom.dtype)
    idx = 0
    for i in range(nbf):
        for j in range(i + 1):
            val = packed_atom[idx]
            full[i, j] = val
            full[j, i] = val
            idx += 1
    return full

def unpack_lower_tri_multi(packed, nbf, natm):
    packed = np.asarray(packed)
    nbf_tri = nbf * (nbf + 1) // 2
    flat = packed.ravel(order="C")
    if flat.size != natm * nbf_tri:
        raise ValueError(f"Size mismatch: got {flat.size}, expected {natm * nbf_tri}")
    packed_by_atom = flat.reshape((natm, nbf_tri))
    full_all = np.zeros((natm, nbf, nbf), dtype=packed.dtype)
    for a in range(natm):
        full_all[a] = unpack_lower_tri_single(packed_by_atom[a], nbf)
    return full_all

def pack_lower_tri_single(full):
    full = np.asarray(full)
    if full.shape[0] != full.shape[1]:
        raise ValueError("Matrix must be square")
    nbf = full.shape[0]
    nbf_tri = nbf * (nbf + 1) // 2
    packed = np.zeros(nbf_tri, dtype=full.dtype)
    idx = 0
    for i in range(nbf):
        for j in range(i + 1):
            packed[idx] = full[i, j]
            idx += 1
    return packed


def read_xyz(filepath):
    """
    Read a standard XYZ file.

    Returns
    -------
    symbols : list of str
    coords : np.ndarray, shape (natom, 3)  –  Angstroms
    """
    symbols = []
    coords = []
    with open(filepath, "r") as f:
        natom = int(f.readline().strip())
        f.readline()  # comment line
        for _ in range(natom):
            parts = f.readline().split()
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)


def _periodic_nonbonded_cutoff(topology, cutoff_method):
    """Return a safe OpenMM nonbonded cutoff for the current periodic box."""
    if cutoff_method is app.NoCutoff:
        return 1.0 * unit.nanometer
    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return 1.0 * unit.nanometer
    lengths = []
    for vec in vectors:
        xyz = vec.value_in_unit(unit.nanometer)
        lengths.append(float(np.linalg.norm(xyz)))
    min_len = min(lengths)
    return min(1.0, 0.4 * min_len) * unit.nanometer


class OpenQpQMMM:
    """
    Low-level QM/MM driver using OpenMM (MM) + OpenQP (QM).

    Two modes:
      1. **Config mode** – ``oqp_cfg`` (dict).
      2. **Mol mode**    – pre-built ``mol`` object.

    Exactly one of ``oqp_cfg`` / ``mol`` must be provided.
    """

    def __init__(
        self,
        positions,
        topology,
        forcefield,
        qm_atoms,
        oqp_cfg=None,
        mol=None,
        Cutoff=app.NoCutoff,
        Embedding='mechanical',
    ):
        if oqp_cfg is None and mol is None:
            raise ValueError("Either 'oqp_cfg' or 'mol' must be provided.")
        if oqp_cfg is not None and mol is not None:
            raise ValueError(
                "'oqp_cfg' and 'mol' are mutually exclusive – provide only one."
            )

        self.positions = positions
        self.topology = topology
        self.forcefield = forcefield
        self.qm_atoms = np.array(qm_atoms, dtype=int)
        self.Cutoff = Cutoff
        self.Embedding = Embedding

        self.use_mol = mol is not None

        if self.use_mol:
            self.mol = mol
            self.oqp_cfg_base = None
        else:
            self.oqp_cfg_base = oqp_cfg
            self.mol = None

        self.op = None
        self.mm_systems = self.prepare_mm()

    # --- Internal helpers -------------------------------------------------

    def _build_xyz_string(self):
        xyz_atoms = []
        for atom in self.topology.atoms():
            at_index = atom.index
            if at_index in self.qm_atoms:
                sym = atom.element.symbol
                x = self.positions[at_index][0].value_in_unit(unit.angstrom)
                y = self.positions[at_index][1].value_in_unit(unit.angstrom)
                z = self.positions[at_index][2].value_in_unit(unit.angstrom)
                xyz_atoms.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
        return '; '.join(xyz_atoms)

    def _update_mol_positions(self):
        coords = []
        for atom in self.topology.atoms():
            if atom.index in self.qm_atoms:
                x = self.positions[atom.index][0].value_in_unit(unit.angstrom)
                y = self.positions[atom.index][1].value_in_unit(unit.angstrom)
                z = self.positions[atom.index][2].value_in_unit(unit.angstrom)
                coords.append([x, y, z])
        coords = np.array(coords)
        ang2bohr = 1.8897259886
        self.mol.set_atoms2("xyz", (coords * ang2bohr).ravel())

    def forces_qm_openqp(self, potmm=None, potqm=None):

        if self.use_mol:
            # ---- Mol mode ------------------------------------------------
            self._update_mol_positions()
            sp = SinglePoint(self.mol)
            sp._prep_guess()

            self.mol.data["OQP::POTMM"] = potmm
            self.mol.data["OQP::POTQM"] = potqm
            oqp.espf_op_corr(self.mol)
            espf_op_corr = self.mol.data["OQP::ESPF_CORR"]

            basis = self.mol.data.get_basis()
            nat = self.mol.data["natom"]
            nbf = basis["nbf"]
            espf_op_corr_f = unpack_lower_tri_multi(espf_op_corr, nbf, nat)

            if potmm is not None:
                hcore = self.mol.get_hcore()
                hcore_full = unpack_lower_tri_single(hcore, nbf)
                hcore_full += np.einsum("ijk,i->jk", espf_op_corr_f, potmm)
                self.mol.set_hcore(pack_lower_tri_single(hcore_full))

            sp.scf()
            self.eqm = self.mol.get_scf_energy()

            oqp.form_esp_charges(self.mol)
            self.pchg_qm = self.mol.data["OQP::partial_charges"]

            if potqm is not None and potmm is not None:
                potmm -= np.einsum(
                    "ij,j->i", potqm,
                    self.pchg_qm - self.mol.get_atoms2("charge")
                )
                self.mol.data["OQP::POTMM"] = potmm

            if potmm is not None:
                self.eqm -= np.dot(
                    self.pchg_qm - self.mol.get_atoms2("charge"), potmm
                )

            self._sp = sp

        else:
            # ---- Config mode ---------------------------------------------
            xyz_atoms = self._build_xyz_string()
            self.oqp_cfg_base["input.system"] = xyz_atoms
            self.op = OPENQP(self.oqp_cfg_base, True)
            self.op.sp._prep_guess()

            self.op.mol.data["OQP::POTMM"] = potmm
            self.op.mol.data["OQP::POTQM"] = potqm
            oqp.espf_op_corr(self.op.mol)
            espf_op_corr = self.op.mol.data["OQP::ESPF_CORR"]

            basis = self.op.mol.data.get_basis()
            nat = self.op.mol.data["natom"]
            nbf = basis["nbf"]
            espf_op_corr_f = unpack_lower_tri_multi(espf_op_corr, nbf, nat)

            if potmm is not None:
                hcore = self.op.mol.get_hcore()
                hcore_full = unpack_lower_tri_single(hcore, nbf)
                hcore_full += np.einsum("ijk,i->jk", espf_op_corr_f, potmm)
                self.op.mol.set_hcore(pack_lower_tri_single(hcore_full))

            self.op.sp.scf()
            self.eqm = self.op.mol.get_scf_energy()

            oqp.form_esp_charges(self.op.mol)
            self.pchg_qm = self.op.mol.data["OQP::partial_charges"]

            if potqm is not None and potmm is not None:
                potmm -= np.einsum(
                    "ij,j->i", potqm,
                    self.pchg_qm - self.op.mol.get_atoms2("charge")
                )
                self.op.mol.data["OQP::POTMM"] = potmm

            if potmm is not None:
                self.eqm -= np.dot(
                    self.pchg_qm - self.op.mol.get_atoms2("charge"), potmm
                )

            # --- Gradients: pure QM + ESPF contribution -----------------------
            gradient = Gradient(self.op.mol)
            if gradient.method == 'hf':
                oqp.hf_gradient(self.op.mol)
                grad = self.op.mol.get_grad()
                gqm = np.array(grad.copy()).reshape(
                    (1, self.op.mol.get_atoms2("natom"), 3)
                )
                oqp.grad_esp_qmmm(self.op.mol)
                gqm += self.op.mol.data["OQP::ESPF_GRAD"]
                # --- Unit conversion to OpenMM conventions ------------------------
                self.eqm *= 2625.5 * unit.kilojoule_per_mole
                self.gqm = gqm[0]*49578.9  # already dimensionless, factor to kJ/mol/nm
            if gradient.method == 'tdhf':
                energies = self.op.sp.excitation([self.eqm])
                grads = np.zeros(( gradient.nstate + 1,  gradient.natom, 3))
                for i in gradient.grads:
                    dump_log(gradient.mol, title='PyOQP: Gradient of Root %s' % i)
                    gradient.mol.data.set_tdhf_target(i)
                    gradient.zvec_func[gradient.td](gradient.mol)

                    # check convergence
                    z_flag = gradient.mol.mol_energy.Z_Vector_converged

                    if not z_flag:
                        dump_log(gradient.mol, title='PyOQP: TD Z-vector is not converged', section='end')

                        if gradient.exception is True:
                            raise ZVnotConverged()
                        else:
                            exit()

                    gradient.grad_func[gradient.td](gradient.mol)
                    gqm = gradient.mol.get_grad().reshape((gradient.natom, 3))
                    oqp.grad_esp_qmmm_excited(self.op.mol)
                    gqm += self.op.mol.data["OQP::ESPF_GRAD"]
                    grads[i] = gqm.copy()
                    self.gqm = gqm*49578.9#49614.78
                    self.eqm = energies[i] * 2625.5 * unit.kilojoule_per_mole
            self.op.mol.save_data()


            print(f"self.gqm:{self.gqm}")

        return self.eqm, self.gqm, self.pchg_qm

    def _get_mol(self):
        if self.use_mol:
            return self.mol
        else:
            return self.op.mol

    def forces_mm(self, pchg_qm):
        system = self.mm_systems["sys0"]
        simulation = self.mm_systems["sim0"]

        forces = { force.__class__.__name__ : force for force in system.getForces() }
        nonbonded = forces['NonbondedForce']
        for k, iatom in enumerate(self.qm_atoms):
            charge, sigma, epsilon = nonbonded.getParticleParameters(iatom)
            charge = pchg_qm[k]*unit.elementary_charge
            nonbonded.setParticleParameters(iatom, charge, sigma, epsilon)

        for i in range(nonbonded.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = nonbonded.getExceptionParameters(i)
            if p1 in self.qm_atoms and p2 not in self.qm_atoms:
               charge1 = pchg_qm[np.where(self.qm_atoms == p1)[0]]*unit.elementary_charge
               charge2, sigma2, epsilon2 = nonbonded.getParticleParameters(p2)
               chargeProd=charge1*charge2
               nonbonded.setExceptionParameters(i,p1,p2,chargeProd,sigma,epsilon)
            elif p2 in self.qm_atoms and p1 not in self.qm_atoms:
               charge2 = pchg_qm[np.where(self.qm_atoms == p2)[0]]*unit.elementary_charge
               charge1, sigma1, epsilon1 = nonbonded.getParticleParameters(p1)
               chargeProd=charge1*charge2
               nonbonded.setExceptionParameters(i,p1,p2,chargeProd,sigma,epsilon)
        nonbonded.updateParametersInContext(simulation.context)

        state = simulation.context.getState(getEnergy=True,getForces=True)
        return state.getPotentialEnergy(), state.getForces(asNumpy=True)

    def compute_force(self, positions, topology, mm_systems, qm_atoms):
        self.positions = positions
        self.topology = topology
        self.mm_systems = mm_systems
        self.qm_atoms = qm_atoms

        potmm = potqm = None
        if self.Embedding == "electrostatic":
            potmm, potqm = self.electrostatic_potential()

        eqm, gqm, pchg_qm = self.forces_qm_openqp(potmm=potmm, potqm=potqm)
        emm, gmm = self.forces_mm(pchg_qm)

        total_energy = eqm + emm
        total_forces = gmm.copy()

        for k, i in enumerate(self.qm_atoms):
            total_forces[i] = total_forces[i] - self.gqm[k]
        cmm=np.sum(total_forces,axis=0)
        for i in range(len(total_forces)):
            total_forces[i]-=cmm/float(len(total_forces))

        return total_energy, total_forces


    def prepare_mm(self):
        positions=self.positions
        topology=self.topology
        forcefield=self.forcefield
        qm_atoms =self.qm_atoms
        qm_set = set(int(i) for i in qm_atoms)
        Cutoff=self.Cutoff
        nb_cutoff = _periodic_nonbonded_cutoff(topology, Cutoff)

        system=forcefield.createSystem(
            topology, nonbondedMethod=Cutoff, nonbondedCutoff=nb_cutoff,
            constraints=None, rigidWater=False)
        nonbonded = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))

        for i in range(nonbonded.getNumExceptions()):
            p1, p2, chgProd, sigma, epsilon = nonbonded.getExceptionParameters(i)
            if (int(p1) in qm_set) or (int(p2) in qm_set):
               nonbonded.setExceptionParameters(i, p1, p2, chgProd, 0.0, 0.0)

        for p1 in qm_atoms:
           for p2 in qm_atoms:
              if p1 != p2:
                 nonbonded.addException(p1,p2,0,0,0,replace=True)

        int0=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
        simulation=app.Simulation(topology, system, int0)
        simulation.context.setPositions(positions)

        #Deactivate QM-QM interactions
        for f in system.getForces():
        #Deactivating non-bonded terms between QM atoms (but keep QM-MM)
           if isinstance(f, mm.NonbondedForce):
              pass
           elif isinstance(f, mm.CustomNonbondedForce):
              for p1 in qm_atoms:
                for p2 in qm_atoms:
                   f.addExclusion(p1,p2)
              f.updateParametersInContext(simulation.context)
        #Deactivating QM bonded terms
           elif isinstance(f, mm.HarmonicBondForce):
              for i in range(f.getNumBonds()):
                 p1, p2, length, k = f.getBondParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms)
                 if exclude: f.setBondParameters(i, p1, p2, length, 0)
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, mm.CustomBondForce):
              for i in range(f.getNumBonds()):
                 p1, p2, parameters = f.getBondParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms)
                 if exclude: f.setBondParameters(i, p1, p2, (0,0))
        #Deactivating QM angle terms
           elif isinstance(f, mm.HarmonicAngleForce):
              for i in range(f.getNumAngles()):
                 p1, p2, p3, angle, k = f.getAngleParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms)
                 if exclude: f.setAngleParameters(i, p1, p2, p3, angle, 0)
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, mm.CustomAngleForce):
              for i in range(f.getNumAngles()):
                 p1, p2, p3, parameters = f.getAngleParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms)
                 if exclude: f.setAngleParameters(i, p1, p2, p3, (0, 0))
              f.updateParametersInContext(simulation.context)
        #Deactivating QM torsion terms
           elif isinstance(f, mm.PeriodicTorsionForce):
              for i in range(f.getNumTorsions()):
                 p1, p2, p3, p4, periodicity, phase, k = f.getTorsionParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms and p4 in qm_atoms)
                 if exclude: f.setTorsionParameters(i, p1, p2, p3, p4, periodicity, phase, 0)
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, mm.CustomTorsionForce):
              for i in range(f.getNumTorsions()):
                 p1, p2, p3, p4, parameters = f.getTorsionParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms and p4 in qm_atoms)
                 if exclude: f.setTorsionParameters(i, p1, p2, p3, p4, (0,0,0))
              f.updateParametersInContext(simulation.context)
           elif isinstance(f, (mm.CMAPTorsionForce)):
              for i in range(f.getNumTorsions()):
                 cmap, p1, p2, p3, p4, q1, q2, q3, q4 = f.getTorsionParameters(i)
                 exclude = (p1 in qm_atoms and p2 in qm_atoms and p3 in qm_atoms and p4 in qm_atoms)
                 exclude = exclude and (q1 in qm_atoms and q2 in qm_atoms and q3 in qm_atoms and q4 in qm_atoms)
                 if exclude: f.setMapParameters(i,cmap.size,0)
              if f.getNumTorsions() != 0: f.updateParametersInContext(simulation.context)
        #Exception, unless CMMotionRemover
           else:
              if not isinstance(f, mm.CMMotionRemover): exit(f"Force not found")

        if Cutoff is not app.NoCutoff:
           sysew=forcefield.createSystem(
               topology, nonbondedMethod=app.Ewald, nonbondedCutoff=nb_cutoff,
               constraints=None, rigidWater=False)
           intew=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
           simew=app.Simulation(topology, sysew, intew)
           simew.context.setPositions(positions)

           sysor=forcefield.createSystem(topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
           intor=mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.001*unit.picoseconds)
           simor=app.Simulation(topology, sysor, intor)
           simor.context.setPositions(positions)
        else:
           sysew = simew = sysor = simor = None
        return {
         "sys0": system,
         "sim0": simulation,
         "sysew": sysew,
         "simew": simew,
         "sysor": sysor,
         "simor": simor,
        }


    def electrostatic_potential(self):

       syspbc=self.mm_systems["sys0"]
       simpbc=self.mm_systems["sim0"]

    #Focus on non-bonded interactions
       forces = { force.__class__.__name__ : force for force in syspbc.getForces() }
       nonbonded = forces['NonbondedForce']
       if nonbonded is None: ValueError(f"Non-bonded interactions are not present, what shall I do?")

    #######################################################################
    #  1. Compute MM potential (needs QM-QM contributions to be removed)  #
    #######################################################################
       potmm=np.zeros((len(self.qm_atoms)))
       for i in range(len(self.qm_atoms)):
           charge, sigma, epsilon = nonbonded.getParticleParameters(self.qm_atoms[i])
           nonbonded.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, sigma, epsilon)
       nonbonded.updateParametersInContext(simpbc.context)
       state = simpbc.context.getState(getEnergy=True)
       e_pbc_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

       for i in range(len(self.qm_atoms)):
           charge, sigma, epsilon = nonbonded.getParticleParameters(self.qm_atoms[i])
           nonbonded.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, sigma, epsilon)
           nonbonded.updateParametersInContext(simpbc.context)
           state = simpbc.context.getState(getEnergy=True)
           e_pbc_qm_charge=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
           potmm[i]=e_pbc_qm_charge-e_pbc_no_qm_charge
           nonbonded.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, sigma, epsilon)

       # Non-periodic embedding has no Ewald QM-QM self-interaction, so the
       # QM-QM correction potential is identically zero. Return a zero matrix
       # (not None) so the Fortran add_potqm_contributions has a valid record.
       if self.Cutoff == app.NoCutoff:
           return potmm, np.zeros((len(self.qm_atoms), len(self.qm_atoms)))

    #######################################################################
    #                 2. Compute QM pair potential                        #
    #######################################################################
       potqm = np.zeros((len(self.qm_atoms), len(self.qm_atoms)))

    # Create an Ewald system
       sysew=self.mm_systems["sysew"]
       simew=self.mm_systems["simew"]
       forces = { force.__class__.__name__ : force for force in sysew.getForces() }
       nonbondedew = forces['NonbondedForce']

    # Create a non-periodic system
       sysor=self.mm_systems["sysor"]
       simor=self.mm_systems["simor"]
       forcesor = { force.__class__.__name__ : force for force in sysor.getForces() }
       nonbondedor = forcesor['NonbondedForce']

    # Create a system with no charges, both QM and MM (Ewald & Original)
       for i in range(sysew.getNumParticles()):
           nonbondedew.setParticleParameters(i, 0.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(i, 0.0*unit.elementary_charge, 0.0, 0.0)
       nonbondedew.updateParametersInContext(simew.context)
       nonbondedor.updateParametersInContext(simor.context)

       state = simew.context.getState(getEnergy=True)
       e_ew_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

       state = simor.context.getState(getEnergy=True)
       e_or_no_qm_charge = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

    #######################################################################
    # 2.1. Compute the QM-QM diagonal potential and remove QM-QM from MM  #
    #         Note 1: the 1/2 factor needs to be corrected later          #
    #         Note 2: here the PME/Ew method                          #
    #######################################################################
       for i in range(len(self.qm_atoms)):
           nonbondedew.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedew.updateParametersInContext(simew.context)
           state = simew.context.getState(getEnergy=True)
           e_ew_qm_chargei = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
           potqm[i,i] = e_ew_qm_chargei - e_ew_no_qm_charge
           potmm[i] -= potqm[i,i] #Remove QM-QM interactions from MM potential
           potqm[i,i] *= 2.0
           nonbondedew.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)

    #######################################################################
    #        2.2. Compute the QM-QM off-diagonal potential                #
    #######################################################################
       for i in range(len(self.qm_atoms)):

           nonbondedew.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(self.qm_atoms[i], 1.0*unit.elementary_charge, 0.0, 0.0)

           for j in range(i+1,len(self.qm_atoms)):

               nonbondedew.setParticleParameters(self.qm_atoms[j], 1.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedew.updateParametersInContext(simew.context)
               state = simew.context.getState(getEnergy=True)
               e_ew_qm_chargeij = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879

               nonbondedor.setParticleParameters(self.qm_atoms[j], 1.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedor.updateParametersInContext(simor.context)
               state = simor.context.getState(getEnergy=True)
               e_or_qm_chargeij = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*0.000380879
               ecorr = e_or_qm_chargeij - e_or_no_qm_charge

               potij = e_ew_qm_chargeij - e_ew_no_qm_charge - ecorr - 0.5*(potqm[j,j] + potqm[i,i])
               potqm[i,j] = potqm[j,i] = potij

               nonbondedew.setParticleParameters(self.qm_atoms[j], 0.0*unit.elementary_charge, 0.0, 0.0)
               nonbondedor.setParticleParameters(self.qm_atoms[j], 0.0*unit.elementary_charge, 0.0, 0.0)

           nonbondedew.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)
           nonbondedor.setParticleParameters(self.qm_atoms[i], 0.0*unit.elementary_charge, 0.0, 0.0)

       return potmm,potqm
