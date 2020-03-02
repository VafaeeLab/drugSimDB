#!/bin/sh

#PBS -N multi_jobs_GO_KATANA


offset=100
onto=BP

possibleLength=10317
for start in $(seq 1 $((offset)) $((possibleLength))); do
	qsub -N GOSim-$onto-$start -o GOSim-$onto-$start.out -e GOSim-$onto-$start.err  -v start=$start,offset=$offset,onto=$onto /srv/scratch/z3526914/DrugRepo/Scripts/parallel_GO_Sim_Katana.sh
done
