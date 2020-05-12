#!/bin/sh

#PBS -N Rev_target_sim_FV
#PBS -o Rev_target_sim_FV.out
#PBS -e Rev_target_sim_FV.err
#PBS -l select=ncpus=32:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae



module load R/3.5.3

Rscript /srv/scratch/z3526914/DrugRepo/Scripts/Rev_target_sim_FV.R >  output-targetJIRev 