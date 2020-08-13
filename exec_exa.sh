#!/bin/bash

remainder=$(( $PMIX_RANK % 2 ))

if [ $remainder -eq 0 ]; then
 export OMP_PLACES="{0},{4},{8},{12},{28:56},{16},{20},{24}"
else
 export OMP_PLACES="{88},{92},{96},{100},{116:56},{104},{108},{112}"
fi

export OMP_PROC_BIND="close,spread,spread"

if [ $PMIX_RANK -eq 0 ]; then
 nvprof --openmp-profiling off -o trace.%q{OMPI_COMM_WORLD_RANK} ./Qforce.x
elif [ $PMIX_RANK -eq 127 ]; then
 nvprof --openmp-profiling off -o trace.%q{OMPI_COMM_WORLD_RANK} ./Qforce.x
elif [ $PMIX_RANK -eq 255 ]; then
 nvprof --openmp-profiling off -o trace.%q{OMPI_COMM_WORLD_RANK} ./Qforce.x
elif [ $PMIX_RANK -eq 256 ]; then
 nvprof --openmp-profiling off -o trace.%q{OMPI_COMM_WORLD_RANK} ./Qforce.x
else
 exec ./Qforce.x
fi
