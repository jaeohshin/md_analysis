#!/usr/bin/env python
# coding: utf-8

# # Molecular dynamics simulation

# ### 1. Import dependencies

# In[140]:


import copy
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


# In[141]:


# Simulation Options
steps = 2e8                 ## 1e8= 200 ns
write_interval = 1e4        ## 1e4 = 20 ps
log_interval = 1e4          ## 1e4 = 20 ps


# 
# ### 2. Load PDB file

# In[142]:


# create data directory if not exists
#HERE = Path(_dh[-1])
HERE = Path.cwd()
DATA = HERE / "test"
pdb_path = DATA / "cutoff_protein.pdb"

print("PDB Path:")
print(pdb_path)


# ### 3. Prepare the protein

# #### Protein preparation
# 
# A crucial part for successful simulation is a correct and complete system. Crystallographic structures retrieved from the Protein Data Bank often miss atoms, mainly hydrogens, and may contain non-standard residues. In this talktorial, we will use the Python package [PDBFixer](https://github.com/openmm/pdbfixer) to prepare the protein structure. However, co-crystallized ligands are not handled well by [PDBFixer](https://github.com/openmm/pdbfixer) and will thus be prepared separately.

# In[143]:


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


# In[144]:


# prepare protein and build only missing non-terminal residues
prepared_protein = prepare_protein(pdb_path, ignore_missing_residues=False)


# ### Check what has changed after fixing

# ### 4. Merge protein and ligand
# 
# In the next step, we want to merge the prepared protein and ligand structures using the Python package [MDTraj](https://github.com/mdtraj/mdtraj). [MDTraj](https://github.com/mdtraj/mdtraj) can handle the prepared protein, which is currently a [PDBFixer](https://github.com/openmm/pdbfixer) molecule, a format that has a topology and atom positions similar to and usually interchangeable with [OpenMM Modeller](http://docs.openmm.org/latest/userguide/application.html#model-building-and-editing) topologies and positions. For the ligand however, we need to do several conversions, since it is currently an [RDKit](https://github.com/rdkit/rdkit) molecule.

# Now protein and ligand are both in [OpenMM](https://github.com/openmm/openmm) like formats and can be merged with [MDTraj](https://github.com/mdtraj/mdtraj).

# In[145]:


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


# In[146]:


complex_topology, complex_positions = merge_protein(prepared_protein)


# In[147]:


print("Complex topology has", complex_topology.getNumAtoms(), "atoms.")
# NBVAL_CHECK_OUTPUT


# ### 5. MD simulation set up
# 
# 

# #### Force field
# 
# Common force fields like AMBER have parameters for amino acids, nucleic acids, water and ions and usually offer several options to choose from depending on your aim. We use the `amber99sbildn.xml` force field file. For solvation we use the standard three-site [water model](https://en.wikipedia.org/wiki/Water_model) [**TIP3P**](https://aip.scitation.org/doi/10.1063/1.445869).

# In[148]:


"""
I plan to use the amber99sb-ildn force field for the protein, TIP3P for water molecules, 
and the Joung-Cheatham (JC) model for ions. Since the "amber99sbildn.xml" file already includes ion parameters, 
I have commented out the ion potentials in this file to avoid conflicts. The JC parameters are already provided in the "amber14/tip3p.xml" file.

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


# ### 6. Define the parameters

# In[149]:


# Integration Options
dt = 2.0 * unit.femtoseconds
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picoseconds
pressure = 1.0 * unit.atmosphere
barostatInterval = 25 # default value




#platform = mm.Platform.getPlatformByName('CUDA')
#platformProperties = {'Precision': 'mixed'}
hydrogenMass = 1.5*unit.amu


# ### 7. System setup

# With our configured force field we can now  use the  [OpenMM Modeller](http://docs.openmm.org/latest/userguide/application.html#model-building-and-editing) class to create the MD environment, a simulation box which contains the complex and is filled with a solvent. The standard solvent is water with a specified amount of ions. The size of the box can be determined in various ways. We define it with a padding, which results in a cubic box with dimensions dependent on the largest dimension of the complex.
# 
# > Note this step can take a long time, in the order of minutes, depending on your hardware.

# In[150]:


modeller = app.Modeller(complex_topology, complex_positions)
modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers, boxShape='dodecahedron', ionicStrength=0.15 * unit.molar)

print('Save the output file...')
output_file = 'protein_with_solvent2.pdb'
with open(output_file, 'w') as f:
    app.PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(), f)


# With our solvated system and force field, we can finally create an [OpenMM System](http://docs.openmm.org/development/api-python/generated/openmm.openmm.System.html#openmm.openmm.System) and set up the simulation.
# Additionally to the system the simulation needs an integrator. An [OpenMM Integrator](http://docs.openmm.org/development/api-python/library.html#integrators) defines a method for simulating a system by integrating the equations of motion. The chosen **Langevin Integrator** uses Langevin equations. A list of all different kinds of integrators can be found in the [OpenMM Docs](http://docs.openmm.org/development/api-python/library.html#integrators). For further insight into the **Langevin Integrator**, we recommend reading about Langevin equations, e.g. on [Wikipedia](https://en.wikipedia.org/wiki/Langevin_equation).

# In[151]:


print('Building system...')


system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME)
integrator = mm.LangevinMiddleIntegrator(temperature, friction, dt)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)


# ### 8. Minimize energy
# Now that everything is set up, we can perform the simulation. We need to set starting positions and minimize the energy of the system to get a low energy starting configuration, which is important to decrease the chance of simulation failures due to severe atom clashes. The energy minimized system is saved.

# In[152]:


print("Minimizing energy...")
simulation.minimizeEnergy()

#simulation.minimizeEnergy(tolerance=10.0 * unit.kilojoule_per_mole/unit.nanometer, maxIterations=1000)


# ### 9. Equilibration: NVT

# ### NPT equilibration

# In[153]:


system.addForce(mm.MonteCarloBarostat(pressure,temperature, barostatInterval))
simulation.context.reinitialize(preserveState=True)
print("Running NPT ...")
simulation.step(10000)


# In[154]:


# Create a filtered topology excluding water molecules
"""
filtered_topology = mm.app.Topology()
filtered_topology.setPeriodicBoxVectors(simulation.topology.getPeriodicBoxVectors())

# Copy over only non-water chains/residues
filtered_topology = simulation.topology.toPDBTopology(lambda res: res.name not in ["HOH", "WAT", "TIP3"])

###
"""


with open(DATA / "topology_temp.pdb", "w") as pdb_file:
    app.PDBFile.writeFile(
        simulation.topology,
        simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
        file=pdb_file,
        keepIds=True,
    )


# ### Once the minimization has finished, we can perform the MD simulation. 

# In[155]:


import sys

# output settings
simulation.reporters.append(
    md.reporters.XTCReporter(file=str(DATA / "trajectory_temp.xtc"), reportInterval=write_interval)
)
simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout,
        log_interval,
        step=True,
        potentialEnergy=True,
        temperature=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=steps,
        separator="\t",
    )
)


# The velocities for all particles in the system are randomly chosen from a distribution at the given temperature. We chose 300 Kelvin, which is some degrees above room temperature.
# A random seed is generated, but could be explicitly given to reproduce results.
# 
# Then the simulation is performed by taking the steps defined before.

# In[156]:


simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(steps)  # perform the simulation


# In[55]:


# Check the trajectory exists and is not empty
(DATA / "trajectory.xtc").stat().st_size > 0
# NBVAL_CHECK_OUTPUT

