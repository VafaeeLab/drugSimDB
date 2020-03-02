#!/bin/sh

#PBS -N PPI_GO_sim
#PBS -o PPI_GO_sim.out
#PBS -e PPI_GO_sim.err
#PBS -l select=ncpus=8:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae



module load R/3.5.3

Rscript /srv/scratch/z3526914/Host_Path_Interaction/Regulon_Enrichment/update_target_Pathway_PPI_GO_similarity.R