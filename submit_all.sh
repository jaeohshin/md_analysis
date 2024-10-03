#!/bin/bash

# Define an array of directories
directories=("ABL" "AKT1" "BRAF")

# Loop through each directory and submit the SLURM job
for dir in "${directories[@]}"; do
    echo "Submitting job in directory: $dir"
    cd "$dir" || exit  # Change to the directory or exit if it fails
    sbatch ./run.sh  # Submit the job
    cd ..  # Go back to the parent directory
done

