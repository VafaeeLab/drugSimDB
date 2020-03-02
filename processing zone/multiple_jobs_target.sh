#!/bin/sh

offset=200
possibleLength=7026
for start in $(seq 1 $((offset)) $((possibleLength))); do
	qsub -N targetSimNew-$start -o targetSimNew-$start.out -e targetSimNew-$start.err  -v start=$start,offset=$offset /short/yr31/aa7970/azData/DD/Scripts/newTargetSim.sh
done
