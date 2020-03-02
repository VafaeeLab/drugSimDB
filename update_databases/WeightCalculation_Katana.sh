#!/bin/sh

#PBS -N DR_WeightCal_katana
#PBS -o DR_WeightCal.out
#PBS -e DR_WeightCal.err
#PBS -l select=ncpus=8:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae


module load R/3.5.3

Rscript /srv/scratch/z3526914/DrugRepo/Scripts/WeightCalculation_Katana.R \
	
