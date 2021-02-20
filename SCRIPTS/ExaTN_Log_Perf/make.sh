#!/bin/bash
rm *.o *.mod *.x
gfortran -c -O3 stsubs.F90
gfortran -c -O3 combinatoric.F90
gfortran -c -O3 parse_prim.F90
gfortran -c -O3 log_perf_analysis.F90
gfortran log_perf_analysis.o parse_prim.o combinatoric.o stsubs.o -o log_perf_analysis.x
