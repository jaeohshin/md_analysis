import argparse
from pathlib import Path
import pdbfixer
from openmm.app import PDBFile  # Needed to save the fixed PDB

# Set up argument parser
parser = argparse.ArgumentParser(description="Fix a PDB file in the current directory.")
parser.add_argument("filename", type=str, help="PDB file name")
parser.add_argument("--ph", type=float, default=7.0, help="pH level for adding hydrogens (default: 7.0)")

# Parse the arguments
args = parser.parse_args()

# Define PDB path in the current directory
pdb_path = Path(args.filename)

# Load PDB file with PDBFixer
fixer = pdbfixer.PDBFixer(str(pdb_path))
fixer.removeHeterogens()  # Remove co-crystallized ligands
fixer.findMissingResidues()  # Identify missing residues
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(args.ph)

# Save the fixed PDB file
output_path = Path(f"fixed_{args.filename}")
with open(output_path, 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
print(f"Fixed PDB saved as: {output_path}")

