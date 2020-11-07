#!/bin/bash

# MIT License
# 
# Copyright (c) 2020 Jean Nunes Ribeiro Araujo
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

###
### This script runs the same experiments described in the PhD qualification.
###

# number of islands
if [ -z "$1" ]
  then
    islands=5
else
	islands="$1"
fi


# dimension
if [ -z "$2" ]
  then
    dimi=(50, 100)
else
	dimi="$2"
fi


# function
if [ -z "$3" ]
  then
    funi="$(seq 1 30)"
else
	funi="$3"
fi


# number of runs
if [ -z "$4" ]
  then
    runi="$(seq 1 1)"
else
	runi="$(seq 1 "$4")"
fi


# strategies
# 0 - S1 
# 1 = S2
# 2 = FIX
# 3 = PROB
if [ -z "$5" ]
  then
    meti=(0,1,2,3)
else
	meti="$5"
fi

if [ $islands = '-h' ]
	then
		echo "./run arg1 arg2 arg3 arg4 arg5"
		echo "arg1: Number of islands."
		echo "arg2: Dimension: number of variables. Valid values are: (10, 30, 50, 100)."
		echo "arg3: Function: benchmark function to be solved. Valid values are: 1-30."
		echo "arg4: Number of runs: number of times that function is called." 
		echo "arg5: Migration strategy: Valid values are: 0: S1, 1: S2, 2: FIX, 3: PROB."

else
	# call main function with mpirun
	for met in $meti; do
	for fun in $funi; do
	for dim in ${dimi[@]}; do
	for run in $runi; do
		echo "number of islands: $islands; dimension: $dim; function: $fun; dimension: $dim; method: $met" 
		mpirun -np $islands main $fun $run $dim $met
	done
	done
	done
	done
fi