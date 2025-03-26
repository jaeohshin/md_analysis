""" 
This simple code will add chain id (A).
"""
import os

def add_chain_id(input_file, chain_id="A"):
    # Read the input file and write back with the chain ID added if necessary
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
    
    with open(input_file, 'w') as outfile:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Check if chain identifier is already present (22nd column)
                if line[21] != chain_id:
                    # Insert the chain identifier at the 22nd column without shifting other columns
                    modified_line = line[:21] + chain_id + line[22:]
                    outfile.write(modified_line)
                else:
                    # If the chain ID is already correct, write the line as is
                    outfile.write(line)
            else:
                # Write non-atom lines as is
                outfile.write(line)

# Loop through all the PDB files (frame_re_1.pdb to frame_re_9.pdb)
for i in range(1, 101):  # Adjust the range if needed
    input_pdb = f'frame_re_{i}.pdb'
    if os.path.exists(input_pdb):
        add_chain_id(input_pdb)
        print(f'Processed: {input_pdb}')
    else:
        print(f'File {input_pdb} does not exist.')

