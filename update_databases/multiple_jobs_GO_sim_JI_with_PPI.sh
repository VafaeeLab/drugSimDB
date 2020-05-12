#!/bin/sh

#PBS -N multiple_jobs_GO_sim_JI_with_PPI


offset=100
onto=CC

possibleLength=10317
for start in $(seq 1 $((offset)) $((possibleLength))); do
	qsub -N GOSim_withPPI-$onto-$start -o GOSim_withPPI-$onto-$start.out -e GOSim_withPPI-$onto-$start.err  -v start=$start,offset=$offset,onto=$onto /srv/scratch/z3526914/DrugRepo/Scripts/parallel_GO_sim_JI_with_PPI.sh
done
