#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=rome
#SBATCH --time 4:00
#SBATCH --mem=2G

# Parse argument
CMD="$@"

# Run
$CMD

