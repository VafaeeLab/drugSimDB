#!/bin/sh

offset=100
possibleLength=4986
for start in $(seq 1 $((offset)) $((possibleLength))); do
	qsub -N targetSim-$start -o targetSim-$start.out -e targetSim-$start.err  -v start=$start,offset=$offset /short/yr31/aa7970/azData/DD/Scripts/parallel_target_Sim_Raijin.sh
done
