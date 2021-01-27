#!/bin/bash -l
#SBATCH -A snic2021-22-59
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J 01_cellranger_mkfastq
#SBATCH --mail-type=ALL
#SBATCH --mail-user alva.annett.4036@student.uu.se

#Load modules
module load bioinfo-tools
module load cellranger

#--- mkfastq -----------
cellranger mkfastq \
 --id cellranger_test \
 --run data/cellranger-tiny-bcl-1.2.0 \
 --csv data/cellranger-tiny-bcl-simple-1.2.0.csv

