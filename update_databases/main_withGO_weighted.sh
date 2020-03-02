#!/bin/sh

#PBS -N weighted_GO_katana
#PBS -o weighted_GO_katana.out
#PBS -e weighted_GO_katana.err
#PBS -l select=ncpus=8:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae


module load R/3.5.3

Rscript /srv/scratch/z3526914/DrugRepo/Scripts/main_withGO_weighted.R \
	
