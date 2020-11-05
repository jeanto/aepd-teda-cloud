#!/bin/bash

islands=5
funi=18
numruns=0

for fun in $(seq 18 $funi); do
for dim in 100; do
for run in $(seq $numruns $numruns); do
	echo "run: $run; fun: $fun; dim: $dim; islands: $islands" 
	mpirun -np $islands main $fun $run $dim
done
done
done
