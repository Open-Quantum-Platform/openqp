"""OQP OpenMM MD class"""
import copy
import oqp
import numpy as np
import scipy as sc
from oqp.library.single_point import SinglePoint, Gradient, LastStep
from oqp.utils.file_utils import dump_log, dump_data
import oqp.utils.qmmm as qmmm 
import openmm as mm
import openmm.app as app
import openmm.unit as unit


class QMMM_MD:
    def __init__(self, mol):
        self.mol = mol
        self.nSteps = mol.config['qmmm']['nsteps']
        self.timeStep = mol.config['qmmm']['timestep']
        self.init_scf = mol.config['optimize']['init_scf']
        self.nstate = mol.config['tdhf']['nstate']
        self.natom = mol.data['natom']
        self.atoms = mol.get_atoms().reshape((self.natom, 1))
        self.sp = SinglePoint(mol)
        self.grad = Gradient(mol)

        dump_log(self.mol, title='PyOQP: Entering QM/MM molecular dynamics')

    def run_md(self):

        pdb = qmmm.pdb0

#Create an empty system and add particles based on PDB topology
        system = mm.System()
        for atom in pdb.topology.atoms():
           system.addParticle(atom.element.mass)

#Add QM/MM force and energy to system as a CustomExternalForce
        oqp_qmmm = mm.CustomExternalForce('fx*x + fy*y + fz*z + qmmm_energy')
        force_index = system.addForce(oqp_qmmm)
        oqp_qmmm.addPerParticleParameter('fx')
        oqp_qmmm.addPerParticleParameter('fy')
        oqp_qmmm.addPerParticleParameter('fz')

#First call to OQP QM.MM
        qmmm_energy, qmmm_grad = self.oqp_qmmm_single_point(0)
        
#Define initial parameters as QM/MM energy and gradient
        number_of_particles = system.getNumParticles()
        qmmm_energy = (qmmm_energy/number_of_particles/qmmm.kj_per_mol_to_hartree)*unit.kilojoules_per_mole
        oqp_qmmm.addGlobalParameter('qmmm_energy', qmmm_energy)
        qmmm_grad = (qmmm_grad/qmmm.kj_per_mol_to_hartree*qmmm.bohr_to_nm) * unit.kilojoules_per_mole/unit.nanometer
        for i in range(number_of_particles):
            oqp_qmmm.addParticle(i, -qmmm_grad[i])
        
#Create simulation integrator
        integrator = mm.VerletIntegrator(self.timeStep*unit.femtoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        
#Add reporters
        simulation.reporters.append(app.PDBReporter('trajectory.pdb', 1))
        oqp_output=open(self.mol.log+'.md','w')
        oqp_output.write(f"\n\nStarting NVE dynamics using Verlet algorithm:\n")
        simulation.reporters.append(app.StateDataReporter(oqp_output, 1, step=True, time=True, separator=" ",potentialEnergy=True, temperature=False, speed=True,elapsedTime=True))
        
        for step in range(self.nSteps):

            simulation.step(1)
        
#Perform OQP calculation and set up QM/MM energy and gradient external parameters
            qmmm_energy, qmmm_grad = self.oqp_qmmm_single_point(step+1)
        
            qmmm_energy = qmmm_energy/number_of_particles*unit.kilojoules_per_mole/qmmm.kj_per_mol_to_hartree
            simulation.context.setParameter('qmmm_energy', qmmm_energy)
        
            qmmm_grad = qmmm_grad /qmmm.kj_per_mol_to_hartree * qmmm.bohr_to_nm * unit.kilojoules_per_mole/unit.nanometer
            for i in qmmm.qm_atoms:
                oqp_qmmm.setParticleParameters(force_index, i, -qmmm_grad[i])
        
#Update OpenMM external force 
            oqp_qmmm.updateParametersInContext(simulation.context)

    def oqp_qmmm_single_point(self,itr):

        dump_log(self.mol, title='PyOQP: QM/MM molecular dynamics (NVE) using Verlet algorithm ')

        num_atoms, x, y, z, q, mass = qmmm.openmm_system()
        coordinates_qm=np.array(x + y + z).reshape((3, num_atoms)).T.reshape(-1)

        self.mol.update_system(coordinates_qm)

        # compute 1e integral
        if itr == 0:
            do_init_scf = True
        else:
            do_init_scf = self.init_scf
            oqp.library.ints_1e(self.mol)

        # compute energy
        energies = self.sp.energy(do_init_scf=do_init_scf)

        # compute QM/MM gradient
        current_xyz = self.mol.get_system().reshape((-1, 3))
        gradient_qm,gradient_mm=qmmm.openmm_gradient(current_xyz,self.mol.data["OQP::partial_charges"])
        self.mol.data["OQP::mm_gradient"]=np.transpose(gradient_qm).tolist()

        self.grad.grads = [qmmm.istate]
        grads = self.grad.gradient()
        self.mol.energies = energies
        self.mol.grads = grads

        gradient_qmmm=qmmm.form_gradient_qmmm(grads,gradient_mm)

        # flatten data
        energy = energies[qmmm.istate]
        grad = gradient_qmmm

        return energy, grad


