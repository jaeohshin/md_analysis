import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for cluster

# === Config ===
input_dir = "./kinase/abl1"  # update this
output_plot = "ca_ca_histogram.png"

files = sorted(glob(f"{input_dir}/frame_*.pdb"))
all_distances = []

for f in files:
    u = mda.Universe(f)
    calphas = u.select_atoms("name CA")
    pos = calphas.positions
    dists = np.linalg.norm(pos[1:] - pos[:-1], axis=1)
    all_distances.extend(dists)

# Save histogram
plt.figure(figsize=(6, 4))
plt.hist(all_distances, bins=np.linspace(3.5, 4.1, 150), edgecolor='black')
plt.xlabel("Cα–Cα distance (Å)")
plt.ylabel("Frequency")
plt.title("Histogram of consecutive Cα distances")
plt.xlim(3.4, 4.2)  # zoom in around expected range
plt.tight_layout()
plt.savefig(output_plot)
print(f"[INFO] Histogram saved to: {output_plot}")

