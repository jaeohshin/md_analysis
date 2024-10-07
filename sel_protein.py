import mdtraj as md
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Select protein from PDB file and save to a new PDB file.')
parser.add_argument('input_pdb', type=str, help='Input PDB file')
parser.add_argument('output_pdb', type=str, help='Output PDB file name')

# Parse the command-line arguments
args = parser.parse_args()

# Load the PDB file
traj = md.load(args.input_pdb)

# Select only the protein atoms
protein_atoms = traj.topology.select('protein')

# Slice the trajectory to only include the protein atoms
protein_traj = traj.atom_slice(protein_atoms)

# Save the new PDB file containing only the protein
protein_traj.save_pdb(args.output_pdb)

print(f'Protein-only PDB saved as {args.output_pdb}')

