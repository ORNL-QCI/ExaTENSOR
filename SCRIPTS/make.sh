#!/bin/bash
rm *.o *.x
gfortran -c -O3 stsubs.F90
gfortran -c -O3 combinatoric.F90
gfortran -c -O3 log_perf_analysis.F90
gfortran log_perf_analysis.o combinatoric.o stsubs.o -o log_perf_analysis.x
