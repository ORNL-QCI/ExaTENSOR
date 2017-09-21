NAME = utility

#ADJUST THE FOLLOWING ACCORDINGLY:
#Cross-compiling wrappers: [WRAP|NOWRAP]:
export WRAP ?= NOWRAP
#Compiler: [GNU|PGI|INTEL|CRAY|IBM]:
export TOOLKIT ?= GNU
#Optimization: [DEV|OPT]:
export BUILD_TYPE ?= OPT
#MPI Library: [NONE|MPICH|OPENMPI]:
export MPILIB ?= NONE
#BLAS: [ATLAS|MKL|ACML|ESSL|NONE]:
export BLASLIB ?= NONE
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

#BLAS library (whichever you have, set one):
# Set this if you do not have a vendor provided BLAS:
export PATH_BLAS_ATLAS ?= /usr/lib
# Set this if you have vendor provided BLAS (choose):
#  MKL BLAS:
export PATH_BLAS_MKL ?= /opt/intel/mkl/lib/intel64
export PATH_BLAS_MKL_DEP ?= /opt/intel/compilers_and_libraries/linux/lib/intel64_lin
#  ACML BLAS:
export PATH_BLAS_ACML ?= /opt/acml/5.3.1/gfortran64_fma4_mp/lib
#  ESSL BLAS:
export PATH_BLAS_ESSL ?= /sw/summitdev/essl/5.5.0/lib64
export PATH_BLAS_ESSL_DEP ?= /sw/summitdev/xl/20161123/xlf/15.1.5/lib

#YOU ARE DONE!


#=================
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
FCOMP = $(FC_$(WRAP))
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
CCOMP = $(CC_$(WRAP))
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
CPPCOMP = $(CPP_$(WRAP))

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
LIB_IBM = -L.
LIB_NOWRAP = $(LIB_$(TOOLKIT))
LIB_WRAP = -L.
ifeq ($(TOOLKIT),PGI)
 LIB = $(LIB_$(WRAP)) -lstdc++
else
 LIB = $(LIB_$(WRAP)) -lstdc++
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

#LINEAR ALGEBRA FLAGS:
LA_LINK_ATLAS = -L$(PATH_BLAS_ATLAS) -lblas -llapack
ifeq ($(TOOLKIT),GNU)
LA_LINK_MKL = -L$(PATH_BLAS_MKL) -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
else
LA_LINK_MKL = -L$(PATH_BLAS_MKL) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl
endif
LA_LINK_ACML = -L$(PATH_BLAS_ACML) -lacml_mp
ifeq ($(BLASLIB),NONE)
LA_LINK_NOWRAP = -L.
else
LA_LINK_NOWRAP = $(LA_LINK_$(BLASLIB))
endif
LA_LINK_WRAP = -L.
LA_LINK = $(LA_LINK_$(WRAP))

#Accelerator support:
ifeq ($(TOOLKIT),IBM)
DF := -WF,
else
DF :=
endif
ifeq ($(BLASLIB),NONE)
NO_BLAS = -DNO_BLAS
else
NO_BLAS :=
endif
NO_GPU = -DNO_GPU
NO_AMD = -DNO_AMD
NO_PHI = -DNO_PHI

#Platform independence:
PIC_FLAG_GNU = -fPIC
PIC_FLAG_PGI = -fpic
PIC_FLAG_INTEL = -fpic
PIC_FLAG_IBM = -qpic=large
PIC_FLAG_CRAY = -hpic
PIC_FLAG = $(PIC_FLAG_$(TOOLKIT))

#C FLAGS:
CFLAGS_DEV = -c -g -O0 -D_DEBUG
CFLAGS_OPT = -c -O3
CFLAGS = $(CFLAGS_$(BUILD_TYPE)) $(NO_GPU) $(NO_AMD) $(NO_PHI) $(NO_BLAS) -D$(EXA_OS) $(PIC_FLAG)

#FORTRAN FLAGS:
FFLAGS_INTEL_DEV = -c -g -O0 -fpp -vec-threshold4 -qopenmp -mkl=parallel
#FFLAGS_INTEL_DEV = -c -g -fpp -vec-threshold4 -openmp
FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -qopenmp -mkl=parallel
#FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -openmp
FFLAGS_CRAY_DEV = -c -g
FFLAGS_CRAY_OPT = -c -O3
FFLAGS_GNU_DEV = -c -fopenmp -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -g -O0
FFLAGS_GNU_OPT = -c -fopenmp -O3
FFLAGS_PGI_DEV = -c -mp -Mcache_align -Mbounds -Mchkptr -Mstandard -g -O0
FFLAGS_PGI_OPT = -c -mp -Mcache_align -Mstandard -O3
FFLAGS_IBM_DEV = -c -qsmp=omp -g -qkeepparm -qcheck -qsigtrap
FFLAGS_IBM_OPT = -c -qsmp=omp -O3
FFLAGS = $(FFLAGS_$(TOOLKIT)_$(BUILD_TYPE)) $(DF)$(NO_GPU) $(DF)$(NO_AMD) $(DF)$(NO_PHI) $(DF)$(NO_BLAS) $(DF)-D$(EXA_OS) $(PIC_FLAG)

#THREADS:
LTHREAD_GNU   = -lgomp
LTHREAD_PGI   = -mp -lpthread
LTHREAD_INTEL = -liomp5
LTHREAD_CRAY  = -L.
LTHREAD_IBM   = -lxlsmp
LTHREAD = $(LTHREAD_$(TOOLKIT))

#LINKING:
LFLAGS = $(LTHREAD) $(MPI_LINK) $(LA_LINK) $(LIB)

OBJS = ./OBJ/dil_basic.o ./OBJ/stsubs.o ./OBJ/timers.o ./OBJ/combinatoric.o ./OBJ/parse_prim.o

$(NAME): $(OBJS)
	ar cr lib$(NAME).a $(OBJS)

./OBJ/dil_basic.o: dil_basic.F90
	mkdir -p ./OBJ
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) dil_basic.F90 -o ./OBJ/dil_basic.o

./OBJ/stsubs.o: stsubs.F90
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) stsubs.F90 -o ./OBJ/stsubs.o

./OBJ/timers.o: timers.F90
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) timers.F90 -o ./OBJ/timers.o

./OBJ/combinatoric.o: combinatoric.F90
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) combinatoric.F90 -o ./OBJ/combinatoric.o

./OBJ/parse_prim.o: parse_prim.F90 ./OBJ/stsubs.o
	$(FCOMP) $(INC) $(MPI_INC) $(FFLAGS) parse_prim.F90 -o ./OBJ/parse_prim.o


.PHONY: clean
clean:
	rm -f *.x *.a ./OBJ/* *.mod *.modmic *.ptx *.log
