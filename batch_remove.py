import os
import sys

def remove_and_renumber(pdb_file, n):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # Get all unique residue numbers by (chain, resnum)
    residue_keys = []
    for line in lines:
        if line.startswith("ATOM"):
            key = (line[21], int(line[22:26]))
            if key not in residue_keys:
                residue_keys.append(key)

    if len(residue_keys) <= n:
        print(f"Skipping {pdb_file} (only {len(residue_keys)} residues)")
        return

    renumber_map = {
        residue_keys[i]: i + 1
        for i in range(n, len(residue_keys))
    }

    output_lines = []
    for line in lines:
        if line.startswith("ATOM"):
            chain_id = line[21]
            resnum = int(line[22:26])
            key = (chain_id, resnum)

            if key in renumber_map:
                new_resnum = renumber_map[key]
                new_line = line[:22] + f"{new_resnum:>4}" + line[26:]
                output_lines.append(new_line)
        else:
            output_lines.append(line)

    with open(pdb_file, 'w') as f:
        f.writelines(output_lines)

# === Main ===
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python batch_remove.py N")
        sys.exit(1)

    N_REMOVE = int(sys.argv[1])

    for i in range(1, 1001):
        fname = f"frame_{i}.pdb"
        if os.path.exists(fname):
            print(f"Processing {fname}...")
            remove_and_renumber(fname, N_REMOVE)
        else:
            print(f"{fname} not found, skipping.")

