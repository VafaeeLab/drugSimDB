#!/bin/bash

#PBS -l select=ncpus=32:mem=100G
#PBS -l walltime=12:00:00
#PBS -m ae


module load R/3.5.3
Rscript /srv/scratch/z3526914/DrugRepo/Scripts/parallel_GO_Sim_Katana.R \
	/srv/scratch/z3526914/DrugRepo/Data/ \
	result.$onto.rds \
	$start \
	$offset \
	$onto \
	32 \
	GO_sim \
	 >  output-$start 