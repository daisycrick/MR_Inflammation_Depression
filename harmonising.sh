#!/bin/bash

#SBATCH --job-name=harmonising
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=15
#SBATCH --time=100:00:00
#SBATCH --mem=100000M

WORK_DIR="/mnt/storage/home/dc15053/PhenoSPD"


cd $WORK_DIR


R CMD BATCH /mnt/storage/home/dc15053/PhenoSPD/harmonising_BMI_MDD.R
