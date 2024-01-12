#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=rome
#SBATCH --time 8:00
#SBATCH --mem=2G

# Parse argument
CMD="$@"

# Run
$CMD

# # Parse arguments
# # CMD="$@"
# bed_file=$1
# cojo_file=$2
# chr=$3
# snp_file=$4
# out_pre=$5
# 
# # Set static arguments
# gcta_binary=/projects/0/prjs0817/software/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
# 
# # Print arguments
# echo "Binary:     ${gcta_binary}"
# echo "BED file:   ${bed_file}"
# echo "COJO file:  ${cojo_file}"
# echo "Chromosome: ${chr}"
# echo "SNP file:   ${snp_file}"
# echo "Out prefix: ${out_pre}"
# 
# # Run
# # $CMD
# ${gcta_binary} \
#   --bfile ${bed_file} \
#   --cojo-file ${cojo_file} \
#   --chr ${chr} \
#   --extract ${snp_file} \
#   --cojo-slct \
#   --out ${out_pre}
    
