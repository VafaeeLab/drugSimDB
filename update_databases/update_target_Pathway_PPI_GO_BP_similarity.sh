#!/bin/sh

#PBS -N PPI_GO_sim_BP
#PBS -o PPI_GO_sim_BP.out
#PBS -e PPI_GO_sim_BP.err
#PBS -l select=ncpus=32:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae



module load R/3.5.3

Rscript /srv/scratch/z3526914/DrugRepo/Scripts/update_target_Pathway_PPI_GO_BP_similarity.R