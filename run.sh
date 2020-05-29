#!/bin/bash

islands=5
fun=9
numruns=1

for dim in 10; do
for run in $(seq 1 $numruns); do
	echo "run: $run; fun: $fun; dim: $dim; islands: $islands" 
	mpirun -np $islands main $fun $run $dim
done
done
