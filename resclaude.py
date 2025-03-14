def renumber_pdb_residues(input_pdb, output_pdb, start_from=0):
    """
    Renumber residues in a PDB file starting from a specific value,
    automatically determining the offset from the file.
    
    Parameters:
    -----------
    input_pdb : str
        Path to the input PDB file
    output_pdb : str
        Path to the output PDB file
    start_from : int, optional
        The new starting residue number (default: 0)
    """
    # First pass: determine the minimum residue number in the file
    min_residue = float('inf')
    
    with open(input_pdb, 'r') as infile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    resnum = int(line[22:26])
                    min_residue = min(min_residue, resnum)
                except ValueError:
                    continue
    
    if min_residue == float('inf'):
        print("No valid residue numbers found in the file.")
        return
    
    # Calculate the offset
    offset = min_residue - start_from
    print(f"Minimum residue number in file: {min_residue}")
    print(f"Applying offset of {offset} to start from {start_from}")
    
    # Second pass: renumber the residues
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    resnum = int(line[22:26])
                    new_resnum = resnum - offset
                    new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                    outfile.write(new_line)
                except ValueError:
                    outfile.write(line)
            else:
                outfile.write(line)
    
    print(f"Renumbered PDB saved to {output_pdb}")

# Example usage
if __name__ == "__main__":
    input_pdb = "init.pdb"
    output_pdb = "renumbered.pdb"
    
    # Renumber residues starting from 0
    renumber_pdb_residues(input_pdb, output_pdb, start_from=0)
