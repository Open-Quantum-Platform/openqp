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
        number_of_particles = system.getNumParticles()

#Add QM/MM gradient and energy to system as a CustomExternalForce
        oqp_qmmm = mm.CustomExternalForce('grad_x*x + grad_y*y + grad_z*z + qmmm_energy - ecorr')
        force_index = system.addForce(oqp_qmmm)
        oqp_qmmm.addPerParticleParameter('grad_x')
        oqp_qmmm.addPerParticleParameter('grad_y')
        oqp_qmmm.addPerParticleParameter('grad_z')
        
#First call to OQP and define initial parameters as QM/MM energy and gradient
        qmmm_energy, qmmm_grad = self.oqp_qmmm_single_point(0)
# Energy
        qmmm_energy = (qmmm_energy/number_of_particles/qmmm.kj_per_mol_to_hartree)*unit.kilojoules_per_mole
        oqp_qmmm.addGlobalParameter('qmmm_energy', qmmm_energy)
# Gradient 
        qmmm_grad = (qmmm_grad/qmmm.kj_per_mol_to_hartree/qmmm.bohr_to_nm) * unit.kilojoules_per_mole/unit.nanometer
        ecorr=0.0*unit.kilojoules_per_mole
        for i in range(number_of_particles):
            oqp_qmmm.addParticle(i, qmmm_grad[i])
            ecorr+=pdb.positions[i][0]*qmmm_grad[i][0]
            ecorr+=pdb.positions[i][1]*qmmm_grad[i][1]
            ecorr+=pdb.positions[i][2]*qmmm_grad[i][2]
        ecorr/=number_of_particles
#This is an energy correction to eliminate the grad_x*x + grad_y*y + grad_z*z contribution
        oqp_qmmm.addGlobalParameter('ecorr', ecorr)

#Create simulation integrator
        integrator = mm.VerletIntegrator(self.timeStep*unit.femtoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300.0*unit.kelvin)

#Add reporters
        simulation.reporters.append(app.PDBReporter('trajectory.pdb', 1))
        oqp_output=open(self.mol.log+'.md','w')
        oqp_output.write(f"\n\nStarting NVE dynamics using Verlet algorithm:\n")
        simulation.reporters.append(app.StateDataReporter(oqp_output, 1, step=True, time=True, separator=" ",totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=False, speed=True,elapsedTime=True))

        state = simulation.context.getState(getPositions=True,getForces=True)

        for step in range(self.nSteps):

            simulation.step(1)

#Get the current positions from simulation context and update the pdb
            state = simulation.context.getState(getPositions=True)
            qmmm.pdb0.positions = state.getPositions()[:number_of_particles]
        
#Perform OQP calculation and update parameters in context 
            qmmm_energy, qmmm_grad = self.oqp_qmmm_single_point(step+1)
            qmmm_energy = qmmm_energy/number_of_particles*unit.kilojoules_per_mole/qmmm.kj_per_mol_to_hartree
            simulation.context.setParameter('qmmm_energy', qmmm_energy)

            qmmm_grad = (qmmm_grad/qmmm.kj_per_mol_to_hartree/qmmm.bohr_to_nm) * unit.kilojoules_per_mole/unit.nanometer
            ecorr=0.0*unit.kilojoules_per_mole
            for i in range(number_of_particles):
                oqp_qmmm.setParticleParameters(i, i, qmmm_grad[i])
                ecorr+=pdb.positions[i][0]*qmmm_grad[i][0]
                ecorr+=pdb.positions[i][1]*qmmm_grad[i][1]
                ecorr+=pdb.positions[i][2]*qmmm_grad[i][2]
            ecorr/=number_of_particles
            simulation.context.setParameter('ecorr', ecorr)
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
        grads = self.grad.gradient()

        # flatten data
        energy = energies[qmmm.istate]
        grad = qmmm.gradient_qmmm

        return energy, grad


