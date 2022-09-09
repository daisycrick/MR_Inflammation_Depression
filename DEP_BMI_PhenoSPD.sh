#!/bin/bash

#SBATCH --job-name=PhenoSpD
#SBATCH --output=/user/home/dc15053/shell_logs/PhenoSpD.o
#SBATCH --error=/user/home/dc150053/shell_logs/PhenoSpD.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=05:00:00
#SBATCH --mem=5000MB

# PhenoSpD
# mkdir /user/home/dc15053/PhenoSPD
MYDIR="/user/home/dc15053/PhenoSPD
cd $MYDIR

# In order to download PhenoSpD, you should clone this repository via the command
# git clone https://github.com/MRCIEU/PhenoSpD.git
# update to the newest version
#git pull

module load languages/r/4.1.0 

# Need following packages:
#install.packages("optparse")
#library(optparse)
#install.packages("installr")
#library(installr)
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("metaCCA"))


cd PhenoSpD
Rscript ./script/phenospd.r --sumstats ./data/harmonised_BMI_DEP.txt --out BMI_DEP_correlation_test

