NAME = gfc

#ADJUST THE FOLLOWING ACCORDINGLY:

#Cross-compiling wrappers: [WRAP|NOWRAP]:
export WRAP ?= NOWRAP
#Compiler: [GNU|PGI|INTEL|CRAY|IBM]:
export TOOLKIT ?= GNU
#Optimization: [DEV|OPT|PRF]:
export BUILD_TYPE ?= OPT
#MPI Library: [NONE|MPICH|OPENMPI]:
export MPILIB ?= NONE
#Operating system: [LINUX|NO_LINUX]:
export EXA_OS ?= LINUX


#SET YOUR LOCAL PATHS (for unwrapped builds):

#MPI library (whichever you have, set one):
# Set this if you have MPICH or its derivative:
export PATH_MPICH ?= /usr/local/mpi/mpich/3.2
#  Only reset these if MPI files are spread in the system directories:
 export PATH_MPICH_INC ?= $(PATH_MPICH)/include
 export PATH_MPICH_LIB ?= $(PATH_MPICH)/lib
 export PATH_MPICH_BIN ?= $(PATH_MPICH)/bin
# Set this if you have OPENMPI or its derivative:
export PATH_OPENMPI ?= /usr/local/mpi/openmpi/2.0.1
#  Only reset these if MPI files are spread in the system directories:
 export PATH_OPENMPI_INC ?= $(PATH_OPENMPI)/include
 export PATH_OPENMPI_LIB ?= $(PATH_OPENMPI)/lib
 export PATH_OPENMPI_BIN ?= $(PATH_OPENMPI)/bin

#YOU ARE DONE!


#=======================
ifeq ($(BUILD_TYPE),PRF)
 COMP_PREF = scorep --thread=omp --cuda
else
 COMP_PREF =
endif
#Fortran compiler:
FC_GNU = gfortran
FC_PGI = pgf90
FC_INTEL = ifort
FC_CRAY = ftn
FC_IBM = xlf2008_r
FC_MPICH = $(PATH_MPICH_BIN)/mpif90
FC_OPENMPI = $(PATH_OPENMPI_BIN)/mpifort
ifeq ($(MPILIB),NONE)
FC_NOWRAP = $(FC_$(TOOLKIT))
else
FC_NOWRAP = $(FC_$(MPILIB))
endif
FC_WRAP = ftn
FCOMP = $(COMP_PREF) $(FC_$(WRAP))
#C compiler:
CC_GNU = gcc
CC_PGI = pgcc
CC_INTEL = icc
CC_CRAY = cc
CC_IBM = xlc_r
CC_MPICH = $(PATH_MPICH_BIN)/mpicc
CC_OPENMPI = $(PATH_OPENMPI_BIN)/mpicc
ifeq ($(MPILIB),NONE)
CC_NOWRAP = $(CC_$(TOOLKIT))
else
CC_NOWRAP = $(CC_$(MPILIB))
endif
CC_WRAP = cc
CCOMP = $(COMP_PREF) $(CC_$(WRAP))
#C++ compiler:
CPP_GNU = g++
CPP_PGI = pgc++
CPP_INTEL = icc
CPP_CRAY = CC
CPP_IBM = xlC_r
CPP_MPICH = $(PATH_MPICH_BIN)/mpic++
ifeq ($(EXA_OS),LINUX)
CPP_OPENMPI = $(PATH_OPENMPI_BIN)/mpic++
else
CPP_OPENMPI = $(PATH_OPENMPI_BIN)/mpicxx
endif
ifeq ($(MPILIB),NONE)
CPP_NOWRAP = $(CPP_$(TOOLKIT))
else
CPP_NOWRAP = $(CPP_$(MPILIB))
endif
CPP_WRAP = CC
CPPCOMP = $(COMP_PREF) $(CPP_$(WRAP))

#COMPILER INCLUDES:
INC_GNU = -I.
INC_PGI = -I.
INC_INTEL = -I.
INC_CRAY = -I.
INC_IBM = -I.
INC_NOWRAP = $(INC_$(TOOLKIT))
INC_WRAP = -I.
INC = $(INC_$(WRAP))

#COMPILER LIBS:
LIB_GNU = -L.
LIB_PGI = -L.
LIB_INTEL = -L.
LIB_CRAY = -L.
LIB_IBM = -L. -lxlsmp
LIB_NOWRAP = $(LIB_$(TOOLKIT))
LIB_WRAP = -L.
ifeq ($(EXA_OS),LINUX)
 LIB = $(LIB_$(WRAP)) -lstdc++
else
ifeq ($(TOOLKIT),INTEL)
 LIB = $(LIB_$(WRAP)) -lc++
else
 LIB = $(LIB_$(WRAP)) -lstdc++
endif
endif

#MPI INCLUDES:
MPI_INC_MPICH = -I$(PATH_MPICH_INC)
MPI_INC_OPENMPI = -I$(PATH_OPENMPI_INC)
ifeq ($(MPILIB),NONE)
MPI_INC_NOWRAP = -I.
else
MPI_INC_NOWRAP = $(MPI_INC_$(MPILIB))
endif
MPI_INC_WRAP = -I.
MPI_INC = $(MPI_INC_$(WRAP))

#MPI LIBS:
MPI_LINK_MPICH = -L$(PATH_MPICH_LIB)
MPI_LINK_OPENMPI = -L$(PATH_OPENMPI_LIB)
ifeq ($(MPILIB),NONE)
MPI_LINK_NOWRAP = -L.
else
MPI_LINK_NOWRAP = $(MPI_LINK_$(MPILIB))
endif
MPI_LINK_WRAP = -L.
MPI_LINK = $(MPI_LINK_$(WRAP))

#Accelerator support:
ifeq ($(TOOLKIT),IBM)
DF := -WF,
else
DF :=
endif
ifeq ($(TOOLKIT),PGI)
NO_PGI :=
else
NO_PGI := -DNO_PGI
endif
NO_GPU = -DNO_GPU
NO_AMD = -DNO_AMD
NO_PHI = -DNO_PHI

#Platform independence:
PIC_FLAG_GNU = -fPIC
PIC_FLAG_PGI = -fpic
PIC_FLAG_INTEL = -fpic
PIC_FLAG_IBM = -qpic=large
PIC_FLAG_CRAY = -fpic
PIC_FLAG = $(PIC_FLAG_$(TOOLKIT))

#C FLAGS:
CFLAGS_INTEL_DEV = -c -g -O0 -qopenmp -D_DEBUG
CFLAGS_INTEL_OPT = -c -O3 -qopenmp
CFLAGS_INTEL_PRF = -c -g -O3 -qopenmp
CFLAGS_CRAY_DEV = -c -g -O0 -fopenmp -D_DEBUG
CFLAGS_CRAY_OPT = -c -O3 -fopenmp
CFLAGS_CRAY_PRF = -c -g -O3 -fopenmp
CFLAGS_GNU_DEV = -c -g -O0 -fopenmp -D_DEBUG
CFLAGS_GNU_OPT = -c -O3 -fopenmp
CFLAGS_GNU_PRF = -c -g -O3 -fopenmp
CFLAGS_PGI_DEV = -c -g -O0 -D_DEBUG -silent -w
CFLAGS_PGI_OPT = -c -O3 -silent -w -Mnovect
CFLAGS_PGI_PRF = -c -g -O3 -silent -w -Mnovect
CFLAGS_IBM_DEV = -c -g -O0 -qsmp=omp -D_DEBUG
CFLAGS_IBM_OPT = -c -O3 -qsmp=omp
CFLAGS_IBM_PRF = -c -g -O3 -qsmp=omp
CFLAGS = $(CFLAGS_$(TOOLKIT)_$(BUILD_TYPE)) $(NO_GPU) $(NO_AMD) $(NO_PHI) $(NO_BLAS) -D$(EXA_OS) $(PIC_FLAG)

#CPP FLAGS:
CPPFLAGS = $(CFLAGS) -std=c++11

#FORTRAN FLAGS:
FFLAGS_INTEL_DEV = -c -g -O0 -fpp -vec-threshold4 -traceback -qopenmp -mkl=parallel
FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -traceback -qopenmp -mkl=parallel
FFLAGS_INTEL_PRF = -c -g -O3 -fpp -vec-threshold4 -traceback -qopenmp -mkl=parallel
FFLAGS_CRAY_DEV = -c -g -fopenmp -J OBJ
FFLAGS_CRAY_OPT = -c -O3 -fopenmp -J OBJ
FFLAGS_CRAY_PRF = -c -g -O3 -fopenmp -J OBJ
FFLAGS_GNU_DEV = -c -fopenmp -g -Og -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -ffpe-trap=invalid,zero,overflow
FFLAGS_GNU_OPT = -c -fopenmp -O3
FFLAGS_GNU_PRF = -c -fopenmp -g -O3
FFLAGS_PGI_DEV = -c -mp -Mcache_align -Mbounds -Mchkptr -Mstandard -Mallocatable=03 -g -O0
FFLAGS_PGI_OPT = -c -mp -Mcache_align -Mstandard -Mallocatable=03 -O3
FFLAGS_PGI_PRF = -c -mp -Mcache_align -Mstandard -Mallocatable=03 -g -O3
FFLAGS_IBM_DEV = -c -qsmp=noopt -g9 -O0 -qfullpath -qkeepparm -qcheck -qsigtrap -qstackprotect=all
FFLAGS_IBM_OPT = -c -qsmp=omp -O3
FFLAGS_IBM_PRF = -c -qsmp=omp -g -O3
FFLAGS = $(FFLAGS_$(TOOLKIT)_$(BUILD_TYPE)) $(DF)$(NO_GPU) $(DF)$(NO_AMD) $(DF)$(NO_PHI) $(DF)-D$(EXA_OS) $(DF)$(NO_PGI) $(PIC_FLAG)

#THREADS:
LTHREAD_GNU   = -lgomp
LTHREAD_PGI   = -mp -lpthread
LTHREAD_INTEL = -liomp5
LTHREAD_CRAY  = -fopenmp
LTHREAD_IBM   = -lxlsmp
LTHREAD = $(LTHREAD_$(TOOLKIT))

#LINKING:
LFLAGS = $(LTHREAD) $(MPI_LINK) $(LIB)

OBJS =  ./OBJ/dil_basic.o ./OBJ/timers.o ./OBJ/stsubs.o ./OBJ/multords.o ./OBJ/combinatoric.o ./OBJ/stack.o ./OBJ/lists.o ./OBJ/dictionary.o \
	./OBJ/gfc_base.o ./OBJ/gfc_range.o ./OBJ/gfc_bank.o ./OBJ/gfc_vector.o ./OBJ/gfc_list.o ./OBJ/gfc_tree.o ./OBJ/gfc_stack.o \
	./OBJ/gfc_vec_tree.o ./OBJ/gfc_queue.o ./OBJ/gfc_pri_queue.o ./OBJ/gfc_dictionary.o ./OBJ/gfc_hash_map.o ./OBJ/gfc_graph.o

$(NAME): lib$(NAME).a ./OBJ/main.o
	$(FCOMP) ./OBJ/main.o lib$(NAME).a $(LFLAGS) -o test_$(NAME).x

lib$(NAME).a: $(OBJS)
	ar cr lib$(NAME).a $(OBJS)

./OBJ/dil_basic.o: dil_basic.F90
	mkdir -p ./OBJ
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) dil_basic.F90 -o ./OBJ/dil_basic.o

./OBJ/stsubs.o: stsubs.F90
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) stsubs.F90 -o ./OBJ/stsubs.o

./OBJ/multords.o: multords.F90
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) multords.F90 -o ./OBJ/multords.o

./OBJ/combinatoric.o: combinatoric.F90
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) combinatoric.F90 -o ./OBJ/combinatoric.o

./OBJ/timers.o: timers.F90
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) timers.F90 -o ./OBJ/timers.o

./OBJ/stack.o: stack.F90 ./OBJ/dil_basic.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) stack.F90 -o ./OBJ/stack.o

./OBJ/lists.o: lists.F90 ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) lists.F90 -o ./OBJ/lists.o

./OBJ/dictionary.o: dictionary.F90 ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) dictionary.F90 -o ./OBJ/dictionary.o

./OBJ/gfc_base.o: gfc_base.F90 ./OBJ/dil_basic.o ./OBJ/stsubs.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_base.F90 -o ./OBJ/gfc_base.o

./OBJ/gfc_range.o: gfc_range.F90 ./OBJ/gfc_base.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_range.F90 -o ./OBJ/gfc_range.o

./OBJ/gfc_bank.o: gfc_bank.F90 ./OBJ/gfc_base.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_bank.F90 -o ./OBJ/gfc_bank.o

./OBJ/gfc_list.o: gfc_list.F90 ./OBJ/gfc_base.o ./OBJ/gfc_vector.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_list.F90 -o ./OBJ/gfc_list.o

./OBJ/gfc_tree.o: gfc_tree.F90 ./OBJ/gfc_base.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_tree.F90 -o ./OBJ/gfc_tree.o

./OBJ/gfc_stack.o: gfc_stack.F90 ./OBJ/gfc_base.o ./OBJ/gfc_list.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_stack.F90 -o ./OBJ/gfc_stack.o

./OBJ/gfc_vector.o: gfc_vector.F90 ./OBJ/gfc_base.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_vector.F90 -o ./OBJ/gfc_vector.o

./OBJ/gfc_vec_tree.o: gfc_vec_tree.F90 ./OBJ/gfc_base.o ./OBJ/gfc_vector.o ./OBJ/gfc_tree.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_vec_tree.F90 -o ./OBJ/gfc_vec_tree.o

./OBJ/gfc_queue.o: gfc_queue.F90 ./OBJ/gfc_base.o ./OBJ/gfc_list.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_queue.F90 -o ./OBJ/gfc_queue.o

./OBJ/gfc_pri_queue.o: gfc_pri_queue.F90 ./OBJ/gfc_base.o ./OBJ/gfc_tree.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_pri_queue.F90 -o ./OBJ/gfc_pri_queue.o

./OBJ/gfc_dictionary.o: gfc_dictionary.F90 ./OBJ/gfc_base.o ./OBJ/gfc_list.o ./OBJ/gfc_vector.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_dictionary.F90 -o ./OBJ/gfc_dictionary.o

./OBJ/gfc_hash_map.o: gfc_hash_map.F90 ./OBJ/gfc_base.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_hash_map.F90 -o ./OBJ/gfc_hash_map.o

./OBJ/gfc_graph.o: gfc_graph.F90 ./OBJ/gfc_base.o ./OBJ/gfc_vector.o ./OBJ/gfc_list.o ./OBJ/gfc_dictionary.o ./OBJ/combinatoric.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) gfc_graph.F90 -o ./OBJ/gfc_graph.o

./OBJ/main.o: main.F90 ./OBJ/dil_basic.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) main.F90 -o ./OBJ/main.o


.PHONY: clean
clean:
	rm -f *.x *.a ./OBJ/* *.mod *.modmic *.ptx *.log
