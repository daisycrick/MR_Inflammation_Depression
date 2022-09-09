#!/bin/bash

#SBATCH --job-name=harmonising
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=15
#SBATCH --time=100:00:00
#SBATCH --mem=100000M

WORK_DIR="/"


cd $WORK_DIR


R CMD BATCH //harmonising_BMI_MDD.R
