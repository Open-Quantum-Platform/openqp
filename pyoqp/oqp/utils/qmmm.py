import openmm as mm
import openmm.app as app
import openmm.unit as unit
import os
import numpy as np
import math

# Information defined in OQP input
force_field = None
nonbondedMethod = None
constraints = None
rigidWater = False

# List of QM atoms (defined in OQP)
qm_atoms = np.array((),dtype=np.uint64)

# Initialize forcefield informations
pdb_file = None
pdb0 = None
forcefield0 = None
system0 = None

# Initialize gradient information
gradient_qmmm = np.array((),dtype=np.float64)

#Bonded info (generated here)
bonded_info = []

# For computing Link atom's g_factor: Define the array of covalent radii
covalent_radii = [
    0.0, 0.354, 0.849, 1.336, 1.010, 0.838,
    0.757, 0.700, 0.658, 0.668, 0.920, 1.539, 1.421,
    1.244, 1.117, 1.101, 1.064, 1.044, 1.032, 1.953,
    1.761, 1.513, 1.412, 1.402, 1.345, 1.382, 1.270,
    1.241, 1.164, 1.302, 1.193, 1.260, 1.197, 1.211,
    1.190, 1.192, 1.147, 2.260, 2.052, 1.698, 1.564,
    1.473, 1.467, 1.322, 1.478, 1.332, 1.338, 1.386,
    1.403, 1.459, 1.398, 1.407, 1.386, 1.382, 1.267,
    2.570, 2.277, 1.943, 1.841, 1.823, 1.816, 1.801,
    1.780, 1.771, 1.735, 1.732, 1.710, 1.696, 1.673,
    1.660, 1.637, 1.671, 1.611, 1.511, 1.392, 1.372,
    1.372, 1.371, 1.364, 1.262, 1.340, 1.518, 1.459,
    1.512, 1.500, 1.545, 1.420, 2.880, 2.512, 1.983,
    1.721, 1.711, 1.684, 1.666, 1.657, 1.660, 1.801,
    1.761, 1.750, 1.724, 1.712, 1.689, 1.679, 1.698,
    1.850
]

#Conversion factor
kj_per_mol_to_hartree = 0.000380879  # 1 kJ/mol = 0.000380879 Hartree
bohr_to_nm = 0.052917721092  # 1 bohr = 0.052917721092 nm

# Define unit charge and zero charge
unit_charge = 1.0 * unit.elementary_charge
zero_charge = 0.0 * unit.elementary_charge
zero_sigma = 0.0 * unit.nanometer
zero_epsilon = 0.0 * unit.kilojoule_per_mole

#!----------------------------------------------------------------------

def _openmm_from_oqp(force_field,nonbondedMethod,constraints,water):

   if nonbondedMethod == 'NoCutoff':
      electrostatics = app.NoCutoff
   else:
      exit(f"The {nonbondedMethod} type of electrostatics is not yet implemented")

   if constraints is None:
      const = None
   elif constraints == 'HBonds':
      const = app.HBonds
   elif constraints == 'AllBonds':
      const = app.AllBonds
   elif constraints == 'HAngles':
      const = app.HAngles
   else:
      exit(f"Unknown type of constraint {constraints}, select a value between None, HBonds, AllBonds or HAngles")

   return force_field,electrostatics,const,water

#!----------------------------------------------------------------------

def _openmm_atoms_bonded(topology, atom_index1, atom_index2):
#
# Purpose: Check if two atoms are bonded (needed for link atom placement)
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#
    for bond in topology.bonds():
        if (bond[0].index == atom_index1 and bond[1].index == atom_index2) or \
           (bond[0].index == atom_index2 and bond[1].index == atom_index1):
            return True
    return False

#!----------------------------------------------------------------------

def _deactivate_bonded_all(system):
#
# Purpose: Deactivate bonded interactions between all atoms
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#
   #1) Deactivate all bonds
   harmonic_bond_force = None
   for force in system.getForces():
      if isinstance(force, mm.HarmonicBondForce):
         harmonic_bond_force = force
         break
   if harmonic_bond_force is None:
      raise ValueError("No HarmonicBondForce found in the system")
   for i in range(harmonic_bond_force.getNumBonds()):
      particle1, particle2, length, k = harmonic_bond_force.getBondParameters(i)
      harmonic_bond_force.setBondParameters(i, particle1, particle2, length, 0.0 * unit.kilojoule_per_mole / unit.nanometer**2)

   #2) Deactivate all angles
   harmonic_angle_force = None
   for force in system.getForces():
       if isinstance(force, mm.HarmonicAngleForce):
          harmonic_angle_force = force
          break
   if harmonic_angle_force is None:
       raise ValueError("No HarmonicAngleForce found in the system")
   for i in range(harmonic_angle_force.getNumAngles()):
       particle1, particle2, particle3, angle, k = harmonic_angle_force.getAngleParameters(i)
       harmonic_angle_force.setAngleParameters(i, particle1, particle2, particle3, angle, 0.0 * unit.kilojoule_per_mole / unit.radian**2)

   #3) Deactivate all torsions
   periodic_torsion_force = None
   for force in system.getForces():
       if isinstance(force, mm.PeriodicTorsionForce):
           periodic_torsion_force = force
           break
   if periodic_torsion_force is None:
       raise ValueError("No PeriodicTorsionForce found in the system")
   for i in range(periodic_torsion_force.getNumTorsions()):
       particle1, particle2, particle3, particle4, periodicity, phase, k = periodic_torsion_force.getTorsionParameters(i)
       periodic_torsion_force.setTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, 0.0 * unit.kilojoule_per_mole)

#!----------------------------------------------------------------------

def _deactivate_bonded_qm(system):
#
# Purpose: Deactivate bonded interactions between QM atoms
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#
   #1) Deactivate all QM-QM bonds
   harmonic_bond_force = None
   for force in system.getForces():
      if isinstance(force, mm.HarmonicBondForce):
         harmonic_bond_force = force
         break
   if harmonic_bond_force is None: raise ValueError("No HarmonicBondForce found in the system")

   for i in range(harmonic_bond_force.getNumBonds()):
      particle1, particle2, length, k = harmonic_bond_force.getBondParameters(i)
      if particle1 in qm_atoms and particle2 in qm_atoms:
         harmonic_bond_force.setBondParameters(i, particle1, particle2, length, 0.0 * unit.kilojoule_per_mole / unit.nanometer**2)

   #2) Deactivate all QM-QM-QM angles
   harmonic_angle_force = None
   for force in system.getForces():
       if isinstance(force, mm.HarmonicAngleForce):
          harmonic_angle_force = force
          break
   if harmonic_angle_force is None: raise ValueError("No HarmonicAngleForce found in the system")

   for i in range(harmonic_angle_force.getNumAngles()):
       particle1, particle2, particle3, angle, k = harmonic_angle_force.getAngleParameters(i)
       if particle1 in qm_atoms and particle2 in qm_atoms and particle3 in qm_atoms:
           harmonic_angle_force.setAngleParameters(i, particle1, particle2, particle3, angle, 0.0 * unit.kilojoule_per_mole / unit.radian**2)

   #3) Deactivate all QM-QM-QM-X or X-QM-QM-QM torsions (X=QM or MM)
   periodic_torsion_force = None
   for force in system.getForces():
       if isinstance(force, mm.PeriodicTorsionForce):
           periodic_torsion_force = force
           break
   if periodic_torsion_force is None: raise ValueError("No PeriodicTorsionForce found in the system")

   for i in range(periodic_torsion_force.getNumTorsions()):
       particle1, particle2, particle3, particle4, periodicity, phase, k = periodic_torsion_force.getTorsionParameters(i)
       if (particle2 in qm_atoms and particle3 in qm_atoms) and (particle1 in qm_atoms or particle4 in qm_atoms):
           periodic_torsion_force.setTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, 0.0 * unit.kilojoule_per_mole)

#!----------------------------------------------------------------------

def openmm_info():
#
# Purpose: Update coordinates for OpenMM
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#

   pdb = pdb0
   ff,electrostatics,const,water = _openmm_from_oqp(force_field,nonbondedMethod,constraints,rigidWater)
   forcefield = app.ForceField(*ff)
   system = forcefield.createSystem(pdb.topology, nonbondedMethod=electrostatics, constraints=const, rigidWater=water)

   xyz = []
   mass = []
   at_num = []
   natoms = system.getNumParticles()

# Iterate over all particles and print their masses
   for atom in pdb.topology.atoms():
      at_num.append(atom.element.atomic_number)
      mass.append(system.getParticleMass(atom.index).value_in_unit(unit.dalton))
      pos=pdb.positions[atom.index].value_in_unit(unit.bohr)
      xyz.append(pos.x)
      xyz.append(pos.y)
      xyz.append(pos.z)

   return xyz,mass,natoms,at_num

#!----------------------------------------------------------------------

def openmm_dump(pdbfile=None):
#
# Purpose: Write PDB in filename
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#

   app.PDBFile.writeFile(pdb0.topology, pdb0.positions, file=pdbfile)
   for atom in pdb0.topology.atoms():
      i=atom.index
      pos=pdb0.positions[i]

#!----------------------------------------------------------------------

def openmm_update_system(coordinates):
#
# Purpose: Update coordinates for OpenMM
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#

############
## Update coordinates
############
   xyz=coordinates.reshape((-1, 3))
   for atom in pdb0.topology.atoms():
      i=atom.index
      pos=mm.Vec3(xyz[i][0],xyz[i][1],xyz[i][2])*unit.bohr
      pdb0.positions[i]=pos

#!----------------------------------------------------------------------

def openmm_init(atom_list):
#
# Purpose: Initialize atom list, forcefield, pdb and system
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#

### Check if atom list is ok
   if len(atom_list) == 0: exit(f"\nError!! Atom list not defined!\n")

   pdb = pdb0

###########
## Read the PDB and load the forcefield information
############
   if pdb is None:
### The pdb file should be defined in the OQP input and it must exist
      if not os.path.exists(pdb_file): exit(f"\nError!! PDB file is required in QM/MM runs!\n")
      pdb = app.PDBFile(pdb_file)
# Attribute forcefield and create system from pdb topology
   ff,electrostatics,const,water = _openmm_from_oqp(force_field,nonbondedMethod,constraints,rigidWater)
   forcefield = app.ForceField(*ff)
   system = forcefield.createSystem(pdb.topology, nonbondedMethod=electrostatics, constraints=const, rigidWater=water)

#Shift atom list and verify it is not out of bounds
   final_atom_list = np.array((),dtype=np.uint32)
   for i in atom_list:
      if int(i)-1 < 0 or int(i)-1 > system.getNumParticles():
         exit(f"\nError!! Atom list out of bounds!\n")
      final_atom_list=np.append(final_atom_list,int(i)-1)

   return final_atom_list,pdb,forcefield,system

#!----------------------------------------------------------------------

def openmm_system():
#
# Purpose: Main function to pass all information OQP needs to know to perform QM/MM calculation
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#
### Reset data to be returned
   num_atoms = 0 #Number of atoms (QM+link atoms)
   q = [] #Atomic numbers (QM+link atoms)
   mass = [] #Masses (QM+link atoms) in Daltons
   x = [] #X positions (QM+link atoms) in Bohr
   y = [] #Y positions (QM+link atoms) in Bohr
   z = [] #Z positions (QM+link atoms) in Bohr

   pdb = pdb0
   forcefield = forcefield0
   system = system0

# Set up integrator
   integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
   simulation = app.Simulation(pdb.topology, system, integrator)
   simulation.context.setPositions(pdb.positions)

# Get the atom's position
   positions = simulation.context.getState(getPositions=True).getPositions()
### Get the information about QM atoms (qm_info) and their connectivity (bonded_info, used to generate link atom list)
   qm_info = []
   for qm_atom in pdb.topology.atoms():
      if qm_atom.index in qm_atoms:
# Retrieve atomic number, mass and position => save it in qm_info
         qm_number = qm_atom.element.atomic_number
         qm_mass = system.getParticleMass(qm_atom.index).value_in_unit(unit.dalton)
         qm_position = positions[qm_atom.index].value_in_unit(unit.bohr)
         qm_info.append((qm_number,qm_mass,[qm_position[0],qm_position[1],qm_position[2]]))
# Check if the specified atoms are bonded => save it in bonded_info
         if not bonded_info:
            for mm_atom in pdb.topology.atoms():
               if mm_atom.index not in qm_atoms:
                  if _openmm_atoms_bonded(pdb.topology, qm_atom.index, mm_atom.index):
                     mm_number = mm_atom.element.atomic_number
                     g_factor=(covalent_radii[1]+covalent_radii[qm_number])/(covalent_radii[qm_number]+covalent_radii[mm_number])
                     bonded_info.append((qm_atom.index,mm_atom.index,g_factor))

### Compute total number of atoms (QM+Link Atoms)
   num_atoms=len(qm_info)+len(bonded_info)

### Adding qm_info to OQP's vectors
   for i in qm_info:
      q.append(i[0])
      mass.append(i[1])
      x.append(i[2][0])
      y.append(i[2][1])
      z.append(i[2][2])

### Adding Link atoms to OQP's vectors (Warning: For the moment, only hydrogen capping is allowed!)
   linkatom_error = False
   for i in range(len(bonded_info)):
      if len(bonded_info[i]) != 0:
         q.append(1)
         mass.append(1.00782503223)
         qm_index,mm_index,g_factor=bonded_info[i]
         qm_position = positions[qm_index].value_in_unit(unit.bohr)
         mm_position = positions[mm_index].value_in_unit(unit.bohr)
         if g_factor >= 1.0 or g_factor <=0:
             if not linkatom_error:
                 print("Error!! You should reconsider your QM/MM partitioning between:")
             print(f"     - QM({qm_atoms[i]}) and MM({mm_index})")
             linkatom_error = True
         x.append(qm_position[0]+g_factor*(mm_position[0]-qm_position[0]))
         y.append(qm_position[1]+g_factor*(mm_position[1]-qm_position[1]))
         z.append(qm_position[2]+g_factor*(mm_position[2]-qm_position[2]))

# If error in link atom partitioning, exit!
   if linkatom_error: exit()

   return num_atoms, x, y, z, q, mass

#!----------------------------------------------------------------------

def openmm_potential(xyz):
#
# Purpose: Compute the MM potential on the xyz coordinates (QM+Link Atoms)
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#

### Reset data to be returned
   mmpot = np.array(()) #MM potential (QM+link atoms) in Hartree/a.u. of charge

   pdb = pdb0
   ff,electrostatics,const,water = _openmm_from_oqp(force_field,nonbondedMethod,constraints,rigidWater)
   forcefield = app.ForceField(*ff)
   system = forcefield.createSystem(pdb.topology, nonbondedMethod=electrostatics, constraints=const, rigidWater=water)

#  forcefield = forcefield0
#  system = system0

#Focus on non-bonded interactions
   forces = { force.__class__.__name__ : force for force in system.getForces() }
   nonbonded_force = forces['NonbondedForce']
   if nonbonded_force is None: raise ValueError("No NonbondedForce found in the system")

#Set all QM charges to zero and the VdW parameters to zero at first, including exceptions (scaling factors)
   for index in range(nonbonded_force.getNumParticles()):
      if index in qm_atoms:
         nonbonded_force.setParticleParameters(index, zero_charge, zero_sigma, zero_epsilon)
   for i in range(nonbonded_force.getNumExceptions()):
      p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
      if p1 in qm_atoms or p2 in qm_atoms:
         nonbonded_force.setExceptionParameters(i, p1, p2, zero_charge*zero_charge, zero_sigma, zero_epsilon)

#Add link atoms to the list with a zero charge. Add exceptions inherited from QM bonded atom
   num_real_particles=system.getNumParticles()
   for ila in range(len(bonded_info)):
       la_index=num_real_particles+ila
       qm_index,mm_index,g_factor=bonded_info[ila]
       system.addParticle(1.00782503223)
       nonbonded_force.addParticle(zero_charge, zero_sigma, zero_epsilon)
       link_atom_position=mm.Vec3(xyz[len(qm_atoms)+ila][0],
                                  xyz[len(qm_atoms)+ila][1],
                                  xyz[len(qm_atoms)+ila][2])*bohr_to_nm
       pdb.positions.append(link_atom_position*unit.nanometer)
       for i in range(nonbonded_force.getNumExceptions()):
          p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
          if p1 == qm_index:
             nonbonded_force.addException(la_index, p2, zero_charge*zero_charge, zero_sigma, zero_epsilon)
          elif p2 == qm_index:
             nonbonded_force.addException(p1, la_index, zero_charge*zero_charge, zero_sigma, zero_epsilon)

# Create a simulation object with link atom
   integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
   simulation = app.Simulation(pdb.topology, system, integrator)
   simulation.context.setPositions(pdb.positions)

# Reset the pdb.positions to contain only real atoms
   pdb.positions=pdb.positions[:num_real_particles]

#Set positions and compute the total energy
   simulation.context.reinitialize(preserveState=True)  # Reinitialize to apply changes
   state_without_charge = simulation.context.getState(getEnergy=True)
   potential_energy_without_charge = state_without_charge.getPotentialEnergy()
#
# Potential on QM Atoms: Loop over all QM particles to get the potential at each QM center
#
   for qm_index in qm_atoms:
#Get the potential with unit charge on target particle
      nonbonded_force.setParticleParameters(qm_index, unit_charge, zero_sigma, zero_epsilon)
      for i in range(nonbonded_force.getNumExceptions()):
         p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
         if p1 == qm_index and p2 not in qm_atoms:
            charge, sigma, epsilon = nonbonded_force.getParticleParameters(p2)
            nonbonded_force.setExceptionParameters(i, p1, p2, charge*unit_charge, zero_sigma, zero_epsilon)
         elif p2 == qm_index and p1 not in qm_atoms:
            charge, sigma, epsilon = nonbonded_force.getParticleParameters(p1)
            nonbonded_force.setExceptionParameters(i, p1, p2, charge*unit_charge, zero_sigma, zero_epsilon)
      simulation.context.reinitialize(preserveState=True)  # Reinitialize to apply changes

# Compute the electrostatic potential energy contribution at the target particle
      state_with_charge = simulation.context.getState(getEnergy=True)
      potential_energy_with_charge = state_with_charge.getPotentialEnergy()
      electrostatic_energy = potential_energy_with_charge - potential_energy_without_charge
      electrostatic_energy_value = electrostatic_energy.value_in_unit(unit.kilojoule_per_mole)

#Get the potential back to zero
      nonbonded_force.setParticleParameters(qm_index, zero_charge, zero_sigma, zero_epsilon)
      for i in range(nonbonded_force.getNumExceptions()):
         p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
         if  p1 == qm_index or p2 == qm_index:
            nonbonded_force.setExceptionParameters(i, p1, p2, zero_charge*zero_charge, zero_sigma, zero_epsilon)
      simulation.context.reinitialize(preserveState=True)  # Reinitialize to apply changes

# Print the potential at each QM center in Hartrees
      mmpot=np.append(mmpot,[electrostatic_energy_value*kj_per_mol_to_hartree],axis=0)

#
# Potential on Link Atoms: Loop over all link atoms to get the potential at each LA center
#
   for ila in range(len(bonded_info)):
      la_index = num_real_particles+ila
      qm_index,mm_index,g_factor=bonded_info[ila]
#Get the potential with unit charge on target particle
      nonbonded_force.setParticleParameters(la_index, unit_charge, zero_sigma, zero_epsilon)
      for i in range(nonbonded_force.getNumExceptions()):
          p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
          if p1 == la_index and p2 not in qm_atoms:
             charge, sigma, epsilon = nonbonded_force.getParticleParameters(p2)
             nonbonded_force.setExceptionParameters(i, la_index, p2, charge*unit_charge, zero_sigma, zero_epsilon)
          elif p2 == la_index and p1 not in qm_atoms:
             charge, sigma, epsilon = nonbonded_force.getParticleParameters(p1)
             nonbonded_force.setExceptionParameters(i, p1, la_index, charge*unit_charge, zero_sigma, zero_epsilon)
      simulation.context.reinitialize(preserveState=True)  # Reinitialize to apply changes

# Compute the electrostatic potential energy contribution at the target particle
      state_with_charge = simulation.context.getState(getEnergy=True)
      potential_energy_with_charge = state_with_charge.getPotentialEnergy()
      electrostatic_energy = potential_energy_with_charge - potential_energy_without_charge
      electrostatic_energy_value = electrostatic_energy.value_in_unit(unit.kilojoule_per_mole)

#Get the potential back to zero
      nonbonded_force.setParticleParameters(la_index, zero_charge, zero_sigma, zero_epsilon)
      for i in range(nonbonded_force.getNumExceptions()):
         p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
         if p1 == la_index and (p2 not in qm_atoms and p2 < num_real_particles):
            nonbonded_force.setExceptionParameters(i, la_index, p2, zero_charge*zero_charge, zero_sigma, zero_epsilon)
         elif p2 == la_index and (p1 not in qm_atoms and p1 < num_real_particles):
            nonbonded_force.setExceptionParameters(i, p1, la_index, zero_charge*zero_charge, zero_sigma, zero_epsilon)
      simulation.context.reinitialize(preserveState=True)  # Reinitialize to apply changes

# Print the potential at each QM center in Hartrees
      mmpot=np.append(mmpot,[electrostatic_energy_value*kj_per_mol_to_hartree],axis=0)

   return mmpot

#!----------------------------------------------------------------------

def openmm_energy():
#
# Purpose: Get the MM energy components after deactivating QM-QM interactions
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#
   pdb = pdb0
   ff,electrostatics,const,water = _openmm_from_oqp(force_field,nonbondedMethod,constraints,rigidWater)
   forcefield = app.ForceField(*ff)
   system = forcefield.createSystem(pdb.topology, nonbondedMethod=electrostatics, constraints=const, rigidWater=water)

# Deactivate all QM-QM bonded interactions and QM bonded terms including exceptions
   _deactivate_bonded_qm(system)
   forces = { force.__class__.__name__ : force for force in system.getForces() }
   nonbonded_force = forces['NonbondedForce']
   if nonbonded_force is None: raise ValueError("No NonbondedForce found in the system")

   for i in range(nonbonded_force.getNumParticles()):
      charge, sigma, epsilon = nonbonded_force.getParticleParameters(i)
      if i in qm_atoms:
         nonbonded_force.setParticleParameters(i, charge*0.0, sigma, epsilon)

   for i in range(nonbonded_force.getNumExceptions()):
      p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
      if p1 in qm_atoms or p2 in qm_atoms:
	      nonbonded_force.setExceptionParameters(i, p1, p2, chargeProd*0.0, sigma*0.0, epsilon*0.0)

# Enable energy decomposition
   for i, f in enumerate(system.getForces()):
       f.setForceGroup(i)

# Set up integrator
   integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
   simulation = app.Simulation(pdb.topology, system, integrator)
   simulation.context.setPositions(pdb.positions)

# Energy decomposition
   energy=np.zeros((24))
   for i, f in enumerate(system.getForces()):
       state = simulation.context.getState(getEnergy=True, groups={i})
       current_energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*kj_per_mol_to_hartree
       key=f.getName()
#
# This code only works for python version >= 3.10
#
#      match key:
#         case "HarmonicBondForce": energy[0]=current_energy
#         case "HarmonicAngleForce": energy[1]=current_energy
#         case "PeriodicTorsionForce": energy[2]=current_energy
#         case "RBTorsionForce": energy[3]=current_energy
#         case "NonbondedForce": energy[4]=current_energy
#         case "CustomBondForce": energy[5]=current_energy
#         case "CustomAngleForce": energy[6]=current_energy
#         case "CustomTorsionForce": energy[7]=current_energy
#         case "CustomNonbondedForce": energy[8]=current_energy
#         case "CustomExternalForce": energy[9]=current_energy
#         case "GBSAOBCForce": energy[10]=current_energy
#         case "CustomGBForce": energy[11]=current_energy
#         case "AmoebaMultipoleForce": energy[12]=current_energy
#         case "AmoebaVdwForce": energy[13]=current_energy
#         case "AmoebaBondForce": energy[14]=current_energy
#         case "AmoebaAngleForce": energy[15]=current_energy
#         case "AmoebaTorsionForce": energy[16]=current_energy
#         case "AmoebaOutOfPlaneBendForce": energy[17]=current_energy
#         case "AmoebaPiTorsionForce": energy[18]=current_energy
#         case "AmoebaStretchBendForce": energy[19]=current_energy
#         case "AmoebaUreyBradleyForce": energy[20]=current_energy
#         case "AmoebaVdw14Force": energy[21]=current_energy
#         case "CMMotionRemover": energy[22]=current_energy
#         case "CustomCompoundBondForce": energy[23]=current_energy
#         case "MonteCarloBarostat": energy[24]=current_energy
#         case _:
#              exit(f"Unknown {key} Force. Please update the qmmm.py code!")
#
# Substituted by if/elif equivalent:
       if key == "HarmonicBondForce": energy[0] = current_energy
       elif key == "HarmonicAngleForce": energy[1] = current_energy
       elif key == "PeriodicTorsionForce": energy[2] = current_energy
       elif key == "RBTorsionForce": energy[3] = current_energy
       elif key == "NonbondedForce": energy[4] = current_energy
       elif key == "CustomBondForce": energy[5] = current_energy
       elif key == "CustomAngleForce": energy[6] = current_energy
       elif key == "CustomTorsionForce": energy[7] = current_energy
       elif key == "CustomNonbondedForce": energy[8] = current_energy
       elif key == "CustomExternalForce": energy[9] = current_energy
       elif key == "GBSAOBCForce": energy[10] = current_energy
       elif key == "CustomGBForce": energy[11] = current_energy
       elif key == "AmoebaMultipoleForce": energy[12] = current_energy
       elif key == "AmoebaVdwForce": energy[13] = current_energy
       elif key == "AmoebaBondForce": energy[14] = current_energy
       elif key == "AmoebaAngleForce": energy[15] = current_energy
       elif key == "AmoebaTorsionForce": energy[16] = current_energy
       elif key == "AmoebaOutOfPlaneBendForce": energy[17] = current_energy
       elif key == "AmoebaPiTorsionForce": energy[18] = current_energy
       elif key == "AmoebaStretchBendForce": energy[19] = current_energy
       elif key == "AmoebaUreyBradleyForce": energy[20] = current_energy
       elif key == "AmoebaVdw14Force": energy[21] = current_energy
       elif key == "CMMotionRemover": energy[22] = current_energy
       elif key == "CustomCompoundBondForce": energy[23] = current_energy
       elif key == "MonteCarloBarostat": energy[24] = current_energy
       else:
           exit(f"Unknown {key} Force. Please update the qmmm.py code!")

# Remove QM-QM electrostatics
#   Deactivate all bonded interactions
   _deactivate_bonded_all(system)

#   Deactivate all QM-MM and MM-MM nonbonded interactions
   for index in range(nonbonded_force.getNumParticles()):
      if index not in qm_atoms:
         nonbonded_force.setParticleParameters(index, 0.0, 0.0, 0.0)
   for i in range(nonbonded_force.getNumExceptions()):
      p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
      if (p1 not in qm_atoms and p2 not in qm_atoms) or (p1 in qm_atoms and p2 not in qm_atoms) or (p1 not in qm_atoms and p2 in qm_atoms):
         nonbonded_force.setExceptionParameters(i, p1, p2, chargeProd*0.0, sigma*0.0, epsilon*0.0)
   simulation.context.reinitialize(preserveState=True)  # Reinitialize to apply changes
   state = simulation.context.getState(getEnergy=True)
   for i, f in enumerate(system.getForces()):
       state = simulation.context.getState(getEnergy=True, groups={i})
       current_energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)*kj_per_mol_to_hartree
       key=f.getName()
#
# This code only works for python version >= 3.10
#
#      match key:
#         case "HarmonicBondForce": energy[0]-=current_energy
#         case "HarmonicAngleForce": energy[1]-=current_energy
#         case "PeriodicTorsionForce": energy[2]-=current_energy
#         case "RBTorsionForce": energy[3]-=current_energy
#         case "NonbondedForce": energy[4]-=current_energy
#         case "CustomBondForce": energy[5]-=current_energy
#         case "CustomAngleForce": energy[6]-=current_energy
#         case "CustomTorsionForce": energy[7]-=current_energy
#         case "CustomNonbondedForce": energy[8]-=current_energy
#         case "CustomExternalForce": energy[9]-=current_energy
#         case "GBSAOBCForce": energy[10]-=current_energy
#         case "CustomGBForce": energy[11]-=current_energy
#         case "AmoebaMultipoleForce": energy[12]-=current_energy
#         case "AmoebaVdwForce": energy[13]-=current_energy
#         case "AmoebaBondForce": energy[14]-=current_energy
#         case "AmoebaAngleForce": energy[15]-=current_energy
#         case "AmoebaTorsionForce": energy[16]-=current_energy
#         case "AmoebaOutOfPlaneBendForce": energy[17]-=current_energy
#         case "AmoebaPiTorsionForce": energy[18]-=current_energy
#         case "AmoebaStretchBendForce": energy[19]-=current_energy
#         case "AmoebaUreyBradleyForce": energy[20]-=current_energy
#         case "AmoebaVdw14Force": energy[21]-=current_energy
#         case "CMMotionRemover": energy[22]-=current_energy
#         case "CustomCompoundBondForce": energy[23]-=current_energy
#         case "MonteCarloBarostat": energy[24]-=current_energy
#         case _:
#              exit(f"Unknown {key} Force. Please update the qmmm.py code!")
#
# Substituted by if/elif equivalent:
       if key == "HarmonicBondForce": energy[0] -= current_energy
       elif key == "HarmonicAngleForce": energy[1] -= current_energy
       elif key == "PeriodicTorsionForce": energy[2] -= current_energy
       elif key == "RBTorsionForce": energy[3] -= current_energy
       elif key == "NonbondedForce": energy[4] -= current_energy
       elif key == "CustomBondForce": energy[5] -= current_energy
       elif key == "CustomAngleForce": energy[6] -= current_energy
       elif key == "CustomTorsionForce": energy[7] -= current_energy
       elif key == "CustomNonbondedForce": energy[8] -= current_energy
       elif key == "CustomExternalForce": energy[9] -= current_energy
       elif key == "GBSAOBCForce": energy[10] -= current_energy
       elif key == "CustomGBForce": energy[11] -= current_energy
       elif key == "AmoebaMultipoleForce": energy[12] -= current_energy
       elif key == "AmoebaVdwForce": energy[13] -= current_energy
       elif key == "AmoebaBondForce": energy[14] -= current_energy
       elif key == "AmoebaAngleForce": energy[15] -= current_energy
       elif key == "AmoebaTorsionForce": energy[16] -= current_energy
       elif key == "AmoebaOutOfPlaneBendForce": energy[17] -= current_energy
       elif key == "AmoebaPiTorsionForce": energy[18] -= current_energy
       elif key == "AmoebaStretchBendForce": energy[19] -= current_energy
       elif key == "AmoebaUreyBradleyForce": energy[20] -= current_energy
       elif key == "AmoebaVdw14Force": energy[21] -= current_energy
       elif key == "CMMotionRemover": energy[22] -= current_energy
       elif key == "CustomCompoundBondForce": energy[23] -= current_energy
       elif key == "MonteCarloBarostat": energy[24] -= current_energy
       else:
            exit(f"Unknown {key} Force. Please update the qmmm.py code!")

   return energy

#!----------------------------------------------------------------------

def openmm_gradient(xyz,oqp_charges):
#
# Purpose: Get the MM gradient components after deactivating QM-QM interactions
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#
   pdb = pdb0
   ff,electrostatics,const,water = _openmm_from_oqp(force_field,nonbondedMethod,constraints,rigidWater)
   forcefield = app.ForceField(*ff)
   system = forcefield.createSystem(pdb.topology, nonbondedMethod=electrostatics, constraints=const, rigidWater=water)

# Remove QM-QM bonded interactions
   _deactivate_bonded_qm(system)

# Add QM charges to the atom list
   forces = { force.__class__.__name__ : force for force in system.getForces() }
   nonbonded_force = forces['NonbondedForce']
   if nonbonded_force is None: raise ValueError("No NonbondedForce found in the system")
   for index in range(nonbonded_force.getNumParticles()):
      if index in qm_atoms:
         charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
         charge = oqp_charges[np.where(qm_atoms == index)[0][0]]*unit.elementary_charge
         nonbonded_force.setParticleParameters(index, charge, sigma, epsilon)

#Add link atoms to the list with a zero charge. Add exceptions inherited from QM bonded atom
   num_real_particles=system.getNumParticles()
   for ila in range(len(bonded_info)):
       la_index=num_real_particles+ila
       qm_index,mm_index,g_factor=bonded_info[ila]
       system.addParticle(1.00782503223)
       nonbonded_force.addParticle(oqp_charges[len(qm_atoms)+ila],0.0,0.0)
       link_atom_position=mm.Vec3(xyz[len(qm_atoms)+ila][0],
                                  xyz[len(qm_atoms)+ila][1],
                                  xyz[len(qm_atoms)+ila][2])*bohr_to_nm
       pdb.positions.append(link_atom_position*unit.nanometer)

# Create a simulation object with link atom
   integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
   simulation = app.Simulation(pdb.topology, system, integrator)
   simulation.context.setPositions(pdb.positions)

# Reset the pdb.positions to contain only real atoms
   pdb.positions=pdb.positions[:num_real_particles]

# Compute gradient
   state = simulation.context.getState(getForces=True)
   gradient = -1.0*state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.bohr)*kj_per_mol_to_hartree

# Remove QM-QM non-bonded

#   Deactivate all bonded interactions
   _deactivate_bonded_all(system)

#   Deactivate all QM-MM and MM-MM nonbonded interactions, including exceptions
   for index in range(nonbonded_force.getNumParticles()):
      if index not in qm_atoms:
         nonbonded_force.setParticleParameters(index, 0.0, 0.0, 0.0)

   for i in range(nonbonded_force.getNumExceptions()):
      p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
      if (p1 not in qm_atoms and p2 not in qm_atoms) or (p1 in qm_atoms and p2 not in qm_atoms) or (p1 not in qm_atoms and p2 in qm_atoms):
         nonbonded_force.setExceptionParameters(i, p1, p2, chargeProd*0.0, sigma*0.0, epsilon*0.0)

   simulation.context.reinitialize(preserveState=True)  # Reinitialize to apply changes
   state = simulation.context.getState(getForces=True)
   gradient2 = -1.0*state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.bohr)*kj_per_mol_to_hartree

   gradient_mm = gradient - gradient2

   grad_qm=np.zeros((len(qm_atoms)+len(bonded_info),3))

   for i in qm_atoms:
         grad_qm[np.where(i == qm_atoms)[0][0],:]=gradient_mm[i,:]

   for ila in range(len(bonded_info)):
         grad_qm[len(qm_atoms)+ila,:]=gradient_mm[num_real_particles+ila,:]

   return grad_qm,gradient_mm

#!----------------------------------------------------------------------

def form_gradient_qmmm(gradient_qm,gradient_mm):
#
# Purpose: Form the final QM/MM gradient
# Author: Miquel Huix-Rotllant
# Date: 31/07/2024
#
   pdb = pdb0
   ff,electrostatics,const,water = _openmm_from_oqp(force_field,nonbondedMethod,constraints,rigidWater)
   forcefield = app.ForceField(*ff)
   system = forcefield.createSystem(pdb.topology, nonbondedMethod=electrostatics, constraints=const, rigidWater=water)

# For link atoms, first project back the contributions
   num_real_particles=system.getNumParticles()

# Project back link atom contributions to real atoms
   if len(bonded_info) != 0:
      for i in range(len(bonded_info)):
         qm_index,mm_index,g_factor=bonded_info[i]
         gradient_qm[0,np.where(qm_index == qm_atoms)[0][0],:]+=(1.0-g_factor)*gradient_qm[0,len(qm_atoms)+i,:]
         gradient_mm[mm_index,:]+=g_factor*gradient_qm[0,len(qm_atoms)+i,:]

# Form the final QM/MM gradient
   qmmm_gradient=np.zeros((num_real_particles,3),dtype=np.float64)
#  for i in qm_atoms:
#     qmmm_gradient[i,:]=gradient_qm[0,np.where(i == qm_atoms)[0][0],:]
   for i in range(num_real_particles):
      if i not in qm_atoms:
         qmmm_gradient[i,:]=gradient_mm[i,:]
      elif i in qm_atoms:
         qmmm_gradient[i,:]=gradient_qm[0,np.where(i == qm_atoms)[0][0],:]
      else:
         exit(f"Error in form_gradient_qmmm, unknown atom index {i}")

   return qmmm_gradient
