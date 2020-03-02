#!/bin/bash
#PBS -l wd
#PBS -o DB_xml_parse.out
#PBS -e DB_xml_parse.err
#PBS -q normal
#PBS -l mem=100GB
#PBS -l ncpus=16
#PBS -P yr31
#PBS -l walltime=10:30:00

module load R/3.5.1

R -f /short/yr31/aa7970/azData/DD/Scripts/drugs_xml_new.R >  output
