##
## Combine pdb files into one xtc file
import mdtraj as md
import glob
import numpy as np

# Get all PDB files and sort them
pdb_files = sorted(glob.glob('*.pdb'))

# Load the first PDB to use as reference topology
reference = md.load(pdb_files[0])
ref_topology = reference.topology

# Initialize list to store coordinates
all_coords = [reference.xyz]

# Load each subsequent PDB and store coordinates
for pdb_file in pdb_files[1:]:
    try:
        # Try to load with the reference topology
        traj = md.load(pdb_file, top=ref_topology)
        all_coords.append(traj.xyz)
    except Exception as e:
        print(f"Skipping {pdb_file}: {e}")

# Combine all coordinates
combined_xyz = np.concatenate(all_coords)

# Create new trajectory with reference topology
combined_traj = md.Trajectory(
    xyz=combined_xyz,
    topology=ref_topology,
    time=np.arange(len(combined_xyz)) * 0.002  # Adding time info, 2 ps intervals
)

# Save as XTC
combined_traj.save_xtc('combined_trajectory.xtc')
# Save reference topology
reference.save_pdb('reference.pdb')
