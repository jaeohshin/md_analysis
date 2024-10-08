#!/usr/bin/env python
# coding: utf-8

# # Molecular dynamics simulation

# ### 1. Import modules

# In[1]:


import copy
import sys
from pathlib import Path
import requests
from IPython.display import display

import numpy as np
import openmm as mm
import openmm.app as app
from openmm import unit
import mdtraj as md
import pdbfixer

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

from openff.toolkit.topology import Molecule, Topology
from openmmforcefields.generators import GAFFTemplateGenerator


# ### 2. Input Files

# In[2]:


# create data directory if not exists
#HERE = Path(_dh[-1])
HERE = Path.cwd()
DATA = HERE / "test"
pdb_path = DATA / "cutoff_protein.pdb"

print("PDB Path:")
print(pdb_path)


"""
I plan to use the amber99sb-ildn force field for the protein, TIP3P for water molecules, 
and the Joung-Cheatham (JC) model for ions. Since the "amber99sbildn.xml" file already includes ion parameters, 
I have commented out the ion potentials (Na+ and Cl-, since I will use only them anyway) in this file to avoid conflicts. The JC parameters are provided in the "amber14/tip3p.xml" file.

The path to the "amber99sbildn.xml" file on the cluster (under my account) is:
/Users/jaeohshin/miniconda3/envs/md/lib/python3.12/site-packages/openmm/app/data/amber99sbildn.xml
I have removed the parameters for Na+ and Cl- from "amber99sbildn.xml" to ensure that the system uses the ion parameters defined in 
"amber14/tip3p.xml" (as of October 1st, 2024, by Jaeoh Shin).

---
Below is comment from OpenMM website for your infomation.
The solvent model XML files included under the amber14/ directory include both water and ions compatible with that water model, 
so if you mistakenly specify tip3p.xml instead of amber14/tip3p.xml, you run the risk of having ForceField throw an exception 
since tip3p.xml will be missing parameters for ions in your system.
---

"""
protein_ff="amber99sbildn.xml"
solvent_ff="amber14/tip3p.xml"
forcefield = app.ForceField(protein_ff, solvent_ff)


# ### 3. System Configuration

# In[3]:


nonbondedMethod = app.PME
rigidWater = True
hydrogenMass = 2*unit.amu


# Integration Options

dt = 2.0 * unit.femtoseconds
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picoseconds
pressure = 1.0 * unit.atmosphere
barostatInterval = 25 # default value


# Simulation Options
steps = 1e5 #5e8                 ## 1e8= 200 ns, 1e6=2ns
write_interval = 25000        ## 25,000 = 50 ps
log_interval = 5e4          ## 5e4 = 100 ps
equilibrationSteps = 1e5

platform = mm.Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}



dataReporter= app.StateDataReporter(
        sys.stdout, log_interval, step=True,
        potentialEnergy=True, kineticEnergy=True, 
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=steps, separator="\t"
    )


# ### 3. Prepare the protein

# #### Protein preparation
# 
# A crucial part for successful simulation is a correct and complete system. Crystallographic structures retrieved from the Protein Data Bank often miss atoms, mainly hydrogens, and may contain non-standard residues. In this talktorial, we will use the Python package [PDBFixer](https://github.com/openmm/pdbfixer) to prepare the protein structure. However, co-crystallized ligands are not handled well by [PDBFixer](https://github.com/openmm/pdbfixer) and will thus be prepared separately.

# In[4]:


def prepare_protein(
    pdb_file, ignore_missing_residues=True, ignore_terminal_missing_residues=True, ph=7.0
):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens
    return fixer


# In[5]:


# prepare protein and build only missing non-terminal residues
prepared_protein = prepare_protein(pdb_path, ignore_missing_residues=False)


# ### 4. Merge protein and ligand
# 
# In the next step, we want to merge the prepared protein and ligand structures using the Python package [MDTraj](https://github.com/mdtraj/mdtraj). [MDTraj](https://github.com/mdtraj/mdtraj) can handle the prepared protein, which is currently a [PDBFixer](https://github.com/openmm/pdbfixer) molecule, a format that has a topology and atom positions similar to and usually interchangeable with [OpenMM Modeller](http://docs.openmm.org/latest/userguide/application.html#model-building-and-editing) topologies and positions. For the ligand however, we need to do several conversions, since it is currently an [RDKit](https://github.com/rdkit/rdkit) molecule.

# Now protein and ligand are both in [OpenMM](https://github.com/openmm/openmm) like formats and can be merged with [MDTraj](https://github.com/mdtraj/mdtraj).

# In[6]:


def merge_protein(protein):
    """
    Merge two OpenMM objects.

    Parameters
    ----------
    protein: pdbfixer.pdbfixer.PDBFixer
        Protein to merge.
    ligand: openmm.app.Modeller
        Ligand to merge.

    Returns
    -------
    complex_topology: openmm.app.topology.Topology
        The merged topology.
    complex_positions: openmm.unit.quantity.Quantity
        The merged positions.
    """
    # combine topologies
    md_protein_topology = md.Topology.from_openmm(protein.topology)  # using mdtraj for protein top
    complex_topology = md_protein_topology.to_openmm()

    # combine positions
    #total_atoms = len(protein.positions) + len(ligand.positions)
    total_atoms = len(protein.positions) 

    # create an array for storing all atom positions as tupels containing a value and a unit
    # called OpenMM Quantities
    complex_positions = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)
    complex_positions = protein.positions  # add protein positions

    return complex_topology, complex_positions


# In[7]:


complex_topology, complex_positions = merge_protein(prepared_protein)


# In[8]:


print("Complex topology has", complex_topology.getNumAtoms(), "atoms.")


# ### 7. System setup

# With our configured force field we can now  use the  [OpenMM Modeller](http://docs.openmm.org/latest/userguide/application.html#model-building-and-editing) class to create the MD environment, a simulation box which contains the complex and is filled with a solvent. The standard solvent is water with a specified amount of ions. The size of the box can be determined in various ways. We define it with a padding, which results in a cubic box with dimensions dependent on the largest dimension of the complex.
# 
# > Note this step can take a long time, in the order of minutes, depending on your hardware.

# In[9]:


modeller = app.Modeller(complex_topology, complex_positions)
modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers, boxShape='dodecahedron', ionicStrength=0.15 * unit.molar)


# ### Building system
# With our solvated system and force field, we can finally create an [OpenMM System](http://docs.openmm.org/development/api-python/generated/openmm.openmm.System.html#openmm.openmm.System) and set up the simulation.
# Additionally to the system the simulation needs an integrator. An [OpenMM Integrator](http://docs.openmm.org/development/api-python/library.html#integrators) defines a method for simulating a system by integrating the equations of motion. The chosen **Langevin Integrator** uses Langevin equations. A list of all different kinds of integrators can be found in the [OpenMM Docs](http://docs.openmm.org/development/api-python/library.html#integrators). For further insight into the **Langevin Integrator**, we recommend reading about Langevin equations, e.g. on [Wikipedia](https://en.wikipedia.org/wiki/Langevin_equation).

# In[10]:


print('Building system...')

system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbondedMethod, rigidWater=rigidWater, hydrogenMass=hydrogenMass)
system.addForce(mm.MonteCarloBarostat(pressure,temperature, barostatInterval))

integrator = mm.LangevinMiddleIntegrator(temperature, friction, dt)
simulation = app.Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(modeller.positions)


# ### 8. Minimize and Equilibrate
# Now that everything is set up, we can perform the simulation. We need to set starting positions and minimize the energy of the system to get a low energy starting configuration, which is important to decrease the chance of simulation failures due to severe atom clashes. The energy minimized system is saved.

# In[11]:


print("Minimizing energy...")
simulation.minimizeEnergy()

print("Equilibrating...")
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)


# In[12]:


#pdb = ap.PDBFile(pdb_)file)
#print(simulation.topology)
#print(type(simulation))
protein_atom_indices = []

for atom in simulation.topology.atoms():
    # Check if the atom belongs to a protein residue
    if atom.residue.name not in ['HOH', 'NA', 'CL']:  # Exclude common non-protein residues
    #if False:
    #if atom.atom.name in ['ATOM']:  # Exclude common non-protein residues
        protein_atom_indices.append(atom.index)


# In[13]:


print(len(protein_atom_indices))


# In[14]:


# Energy minimized system is saved.

with open(DATA / "top.pdb", "w") as pdb_file:
    #positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(protein_atom_indices)
    positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions=positions, file=pdb_file, keepIds=True)


# In[15]:


# Save only protein system
traj = md.load(DATA / "top.pdb")
protein_atoms = traj.topology.select("protein")
protein_traj = traj.atom_slice(protein_atoms)
protein_traj.save_pdb(DATA / "top_protein.pdb")


# In[ ]:





# In[16]:


# Create a filtered topology 
"""
with open(DATA / "topology_subset.pdb", "w") as pdb_file:
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    app.PDBFile.writeFile(
        simulation.topology,
        state.getPositions(protein_atom_indices),
        file=pdb_file,
        keepIds=True,
    )
"""


# In[17]:


xtcReporter = md.reporters.XTCReporter(file=str(DATA / "traj_protein.xtc"), 
                                       reportInterval=write_interval, atomSubset=protein_atoms)


# In[18]:


print(len(protein_atom_indices))


# In[19]:


pwd


# #### Simulating

# In[20]:


print('Simulating...')
simulation.reporters.append(xtcReporter)
simulation.reporters.append(dataReporter)
simulation.currentStep = 0
simulation.step(steps)

