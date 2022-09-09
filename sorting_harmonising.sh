#!/bin/bash

#SBATCH --job-name=harmonising
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=15
#SBATCH --partition=mrcieu
#SBATCH --time=100:00:00
#SBATCH --mem=10000M

WORK_DIR="/mnt/storage/home/dc15053/PhenoSPD"


cd $WORK_DIR


R CMD BATCH /mnt/storage/home/dc15053/PhenoSPD/sorting_harmonised.R
