#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l mem=10GB
#PBS -l ncpus=32
#PBS -P yr31
#PBS -l walltime=02:00:00

module load R/3.5.1

R -f /short/yr31/aa7970/azData/DD/Scripts/newTargetSim.R /short/yr31/aa7970/azData/DD/Data/ protein_similarity.csv uniprot_links.csv $start $offset tarhetSimNew >  output_newTargetSim.txt
