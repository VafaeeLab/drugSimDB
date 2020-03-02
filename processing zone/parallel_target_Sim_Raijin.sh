#!/bin/bash
#PBS -l wd
#PBS -q normal
#PBS -l mem=150GB
#PBS -l ncpus=32
#PBS -P yr31
#PBS -l walltime=05:00:00

module load R/3.5.1

R -f /short/yr31/aa7970/azData/DD/Scripts/parallel_target_Sim_Raijin.R /short/yr31/aa7970/azData/DD/Data/ protein.fasta $start $offset 122 protSim  >  output-$start
