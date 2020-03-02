#!/bin/sh

#PBS -N update_pathway_similarity_v2
#PBS -o update_pathway_similarity_v2.out
#PBS -e update_pathway_similarity_v2.err
#PBS -l select=ncpus=32:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae



module load R/3.5.3

Rscript /srv/scratch/z3526914/DrugRepo/Scripts/update_pathway_similarity_v2.R