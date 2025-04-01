#!/bin/sh

#SBATCH --job-name=QuantumClone
#SBATCH --ntasks-per-node=6
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH -o RunClustering_1.out
#SBATCH -e RunClustering_1.err
#SBATCH --chdir=/masterthesis_marina/DREAM_benchmarking/QuantumClone
#SBATCH --mem=150gb

ml R/4.4.1-gfbf-2023b
module load R-bundle-CRAN/2024.06-foss-2023b
cd /masterthesis_marina/DREAM_benchmarking/QuantumClone
Rscript "/masterthesis_marina/DREAM_benchmarking/QuantumClone/clustering.R"