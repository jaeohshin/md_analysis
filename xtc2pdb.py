""" Make pdb files extracted from a xtc file.
    The residue index start from 1 (which is necessary for the flowpacker) """

import mdtraj as md
import os

# Load trajectory and topology
traj = md.load('samples.xtc', top='topology.pdb')

# Save each frame as a PDB file
for i, frame in enumerate(traj):
    frame.save(f'frame_{i+1}.pdb')

# Post-process: Renumber residues using a simple script
def renumber_pdb(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        res_map = {}  # Maps old residue numbers to new ones
        new_res_idx = 1  # Start numbering from 1
        
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                old_resnum = int(line[22:26])  # Extract original residue number
                if old_resnum not in res_map:
                    res_map[old_resnum] = new_res_idx
                    new_res_idx += 1
                
                new_resnum = res_map[old_resnum]
                new_line = line[:22] + f"{new_resnum:4d}" + line[26:]  # Replace residue number
                outfile.write(new_line)
            else:
                outfile.write(line)

# Apply renumbering and overwrite the original files
for i in range(1, len(traj) + 1):
    original_file = f'frame_{i}.pdb'
    temp_file = f'frame_{i}_temp.pdb'  # Temporary file
    
    # Renumber the PDB and write to a temporary file
    renumber_pdb(original_file, temp_file)
    
    # Replace the original file with the renumbered one
    os.replace(temp_file, original_file)
