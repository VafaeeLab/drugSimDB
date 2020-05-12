#!/bin/sh

#PBS -N noPPI_GO_sim_JI_CC
#PBS -o noPPI_GO_sim_JI_CC.out
#PBS -e noPPI_GO_sim_JI_CC.err
#PBS -l select=ncpus=32:mem=90G
#PBS -l walltime=12:00:00
#PBS -M akm.azad@unsw.edu.au
#PBS -m ae



module load R/3.5.3

Rscript /srv/scratch/z3526914/DrugRepo/Scripts/Drug_GO_terms_without_PPI.R \
	CC \
	/srv/scratch/z3526914/DrugRepo/Data/external/Enrichr/GO_Cellular_Component_2018.txt \
	0.05 \
	 >  output-CC 