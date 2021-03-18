NAME = ddss

#ADJUST THE FOLLOWING ACCORDINGLY:

#Cray cross-compiling wrappers (only for Cray): [WRAP|NOWRAP]:
export WRAP ?= NOWRAP
#Compiler: [GNU|INTEL|CRAY|IBM|PGI]:
export TOOLKIT ?= GNU
#Optimization: [DEV|OPT|PRF]:
export BUILD_TYPE ?= OPT
#MPI library base: [MPICH|OPENMPI]:
export MPILIB ?= MPICH
#NVIDIA GPU via CUDA: [CUDA|NOCUDA]:
export GPU_CUDA ?= NOCUDA
#NVIDIA GPU architecture (two digits):
export GPU_SM_ARCH ?= 35
#Operating system: [LINUX|NO_LINUX]:
export EXA_OS ?= LINUX
#Only for Linux DEV builds with GNU: [YES|NO]:
export LINUX_GNU_ASAN ?= NO


#SET YOUR LOCAL PATHS (for direct build without Cray compiler wrappers):

#MPI library base (whichever you have, set one):
# Set this if you use MPICH or its derivative (e.g. Cray-MPICH):
export PATH_MPICH ?= /usr/local/mpi/mpich/3.2.1
#  Only reset these if MPICH files are spread in system directories:
 export PATH_MPICH_INC ?= $(PATH_MPICH)/include
 export PATH_MPICH_LIB ?= $(PATH_MPICH)/lib
 export PATH_MPICH_BIN ?= $(PATH_MPICH)/bin
# Set this if you use OPENMPI or its derivative (e.g. IBM Spectrum MPI):
export PATH_OPENMPI ?= /usr/local/mpi/openmpi/3.1.0
#  Only reset these if OPENMPI files are spread in system directories:
 export PATH_OPENMPI_INC ?= $(PATH_OPENMPI)/include
 export PATH_OPENMPI_LIB ?= $(PATH_OPENMPI)/lib
 export PATH_OPENMPI_BIN ?= $(PATH_OPENMPI)/bin

#IBM XL (only set these if you use IBM XL compiler):
export PATH_IBM_XL_CPP ?= /sw/summit/xl/16.1.1-3/xlC/16.1.1/lib
export PATH_IBM_XL_FOR ?= /sw/summit/xl/16.1.1-3/xlf/16.1.1/lib
export PATH_IBM_XL_SMP ?= /sw/summit/xl/16.1.1-3/xlsmp/5.1.1/lib

#CUDA (set these only if you build with CUDA):
export PATH_CUDA ?= /usr/local/cuda
# Only reset these if CUDA files are spread in system directories:
 export PATH_CUDA_INC ?= $(PATH_CUDA)/include
 export PATH_CUDA_LIB ?= $(PATH_CUDA)/lib64
 export PATH_CUDA_BIN ?= $(PATH_CUDA)/bin
# Reset your CUDA Host compiler if needed:
 export CUDA_HOST_COMPILER ?= /usr/bin/g++

#YOU ARE DONE! MAKE IT!


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
CPP_MPICH = $(PATH_MPICH_BIN)/mpicxx
CPP_OPENMPI = $(PATH_OPENMPI_BIN)/mpicxx
ifeq ($(MPILIB),NONE)
CPP_NOWRAP = $(CPP_$(TOOLKIT))
else
CPP_NOWRAP = $(CPP_$(MPILIB))
endif
CPP_WRAP = CC
CPPCOMP = $(COMP_PREF) $(CPP_$(WRAP))
#CUDA compiler:
CUDA_COMP = nvcc

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
ifeq ($(TOOLKIT),IBM)
 LIB = $(LIB_$(WRAP)) -L$(PATH_IBM_XL_CPP) -libmc++ -lstdc++
else
ifeq ($(EXA_OS),LINUX)
 LIB = $(LIB_$(WRAP)) -lstdc++
else
ifeq ($(TOOLKIT),INTEL)
 LIB = $(LIB_$(WRAP)) -lc++
else
 LIB = $(LIB_$(WRAP)) -lstdc++
endif
endif
endif

#SANITIZERS:
ifeq ($(LINUX_GNU_ASAN),YES)
 ASAN_COMPILE = -fsanitize=address -fno-omit-frame-pointer
 ASAN_CXX = -Wno-maybe-uninitialized -Wno-unused-result
else
 ASAN_COMPILE = -I.
 ASAN_CXX = -I.
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

#CUDA INCLUDES:
ifeq ($(GPU_CUDA),CUDA)
CUDA_INC_NOWRAP = -I$(PATH_CUDA_INC)
CUDA_INC_WRAP = -I.
CUDA_INC = $(CUDA_INC_$(WRAP))
else
CUDA_INC = -I.
endif

#CUDA LIBS:
CUDA_LINK_NOWRAP = -L$(PATH_CUDA_LIB) -lcublas -lcudart -lnvToolsExt
CUDA_LINK_WRAP = -lcublas -lcudart -lnvToolsExt
CUDA_LINK_CUDA = $(CUDA_LINK_$(WRAP))
CUDA_LINK_NOCUDA = -L.
CUDA_LINK = $(CUDA_LINK_$(GPU_CUDA))

#Platform independence:
PIC_FLAG_GNU = -fPIC
PIC_FLAG_PGI = -fpic
PIC_FLAG_INTEL = -fpic
PIC_FLAG_IBM = -qpic=large
PIC_FLAG_CRAY = -fpic
PIC_FLAG = $(PIC_FLAG_$(TOOLKIT))
PIC_FLAG_CUDA = $(PIC_FLAG_GNU)

#CUDA FLAGS:
ifeq ($(GPU_CUDA),CUDA)
GPU_SM = sm_$(GPU_SM_ARCH)
GPU_ARCH = $(GPU_SM_ARCH)0
CUDA_HOST_COMPILER ?= /usr/bin/g++
CUDA_HOST_NOWRAP = -I.
CUDA_HOST_WRAP = -I.
CUDA_HOST = $(CUDA_HOST_$(WRAP))
CUDA_FLAGS_DEV = --compile -arch=$(GPU_SM) -g -G -DDEBUG_GPU
CUDA_FLAGS_OPT = --compile -arch=$(GPU_SM) -O3 -lineinfo
CUDA_FLAGS_PRF = --compile -arch=$(GPU_SM) -g -G -O3
CUDA_FLAGS_CUDA = $(CUDA_HOST) $(CUDA_FLAGS_$(BUILD_TYPE)) -D_FORCE_INLINES -D$(EXA_OS) -Xcompiler $(PIC_FLAG_CUDA)
ifeq ($(GPU_FINE_TIMING),YES)
CUDA_FLAGS = $(CUDA_FLAGS_CUDA) -DGPU_FINE_TIMING
else
CUDA_FLAGS = $(CUDA_FLAGS_CUDA)
endif
else
CUDA_FLAGS = -D_FORCE_INLINES
endif

#Accelerator support:
ifeq ($(TOOLKIT),IBM)
DF := -WF,
else
DF :=
endif
ifeq ($(GPU_CUDA),CUDA)
NO_GPU = -DCUDA_ARCH=$(GPU_ARCH)
else
NO_GPU = -DNO_GPU
endif
NO_AMD = -DNO_AMD
NO_PHI = -DNO_PHI

#C FLAGS:
CFLAGS_INTEL_DEV = -c -g -O0 -qopenmp -D_DEBUG
CFLAGS_INTEL_OPT = -c -O3 -qopenmp
CFLAGS_INTEL_PRF = -c -g -O3 -qopenmp
CFLAGS_CRAY_DEV = -c -g -O0 -fopenmp -D_DEBUG
CFLAGS_CRAY_OPT = -c -O3 -fopenmp
CFLAGS_CRAY_PRF = -c -g -O3 -fopenmp
CFLAGS_GNU_DEV = -c -g -O0 -fopenmp -D_DEBUG $(ASAN_COMPILE)
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
CPPFLAGS = $(CFLAGS) -std=c++11 $(ASAN_CXX)

#FORTRAN FLAGS:
FFLAGS_INTEL_DEV = -c -g -O0 -fpp -vec-threshold4 -traceback -qopenmp -mkl=parallel
FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -traceback -qopenmp -mkl=parallel
FFLAGS_INTEL_PRF = -c -g -O3 -fpp -vec-threshold4 -traceback -qopenmp -mkl=parallel
FFLAGS_CRAY_DEV = -c -g -fopenmp -J OBJ
FFLAGS_CRAY_OPT = -c -O3 -fopenmp -J OBJ
FFLAGS_CRAY_PRF = -c -g -O3 -fopenmp -J OBJ
FFLAGS_GNU_DEV = -c -fopenmp -g -Og -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -ffpe-trap=invalid,zero,overflow $(ASAN_COMPILE)
FFLAGS_GNU_OPT = -c -fopenmp -O3
FFLAGS_GNU_PRF = -c -fopenmp -g -O3
FFLAGS_PGI_DEV = -c -Mcache_align -Mbounds -Mchkptr -Mstandard -Mallocatable=03 -g -O0
FFLAGS_PGI_OPT = -c -Mcache_align -Mstandard -Mallocatable=03 -O3
FFLAGS_PGI_PRF = -c -Mcache_align -Mstandard -Mallocatable=03 -g -O3
FFLAGS_IBM_DEV = -c -qsmp=noopt -g9 -O0 -qfullpath -qkeepparm -qcheck -qsigtrap -qstackprotect=all
FFLAGS_IBM_OPT = -c -qsmp=omp -O3
FFLAGS_IBM_PRF = -c -qsmp=omp -g -O3
FFLAGS = $(FFLAGS_$(TOOLKIT)_$(BUILD_TYPE)) $(DF)$(NO_GPU) $(DF)$(NO_AMD) $(DF)$(NO_PHI) $(DF)-D$(EXA_OS) $(PIC_FLAG)

#THREADS:
LTHREAD_GNU   = -fopenmp $(ASAN_COMPILE)
LTHREAD_PGI   = -mp -lpthread
LTHREAD_INTEL = -liomp5
LTHREAD_CRAY  = -fopenmp
LTHREAD_IBM   = -lxlsmp
LTHREAD = $(LTHREAD_$(TOOLKIT))

#LINKING:
LFLAGS = $(LTHREAD) $(MPI_LINK) $(CUDA_LINK) $(LIB)

OBJS =  ./OBJ/dil_basic.o ./OBJ/sys_service.o ./OBJ/c2fortran.o ./OBJ/c2f_ifc.o ./OBJ/stsubs.o ./OBJ/timer.o ./OBJ/timers.o \
	./OBJ/nvtx_profile.o ./OBJ/mpi_fort.o ./OBJ/service_mpi.o ./OBJ/distributed.o ./OBJ/pack_prim.o ./OBJ/test_pack_prim.o

$(NAME): lib$(NAME).a ./OBJ/main.o ./OBJ/test_pack_prim.o
	$(FCOMP) ./OBJ/main.o lib$(NAME).a $(LFLAGS) -o test_$(NAME).x
	$(FCOMP) ./OBJ/test_pack_prim.o ./OBJ/pack_prim.o ./OBJ/stsubs.o ./OBJ/dil_basic.o $(LFLAGS) -o test_pack_prim.x

lib$(NAME).a: $(OBJS)
	ar cr lib$(NAME).a $(OBJS)

./OBJ/dil_basic.o: dil_basic.F90
	mkdir -p ./OBJ
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dil_basic.F90 -o ./OBJ/dil_basic.o

./OBJ/sys_service.o: sys_service.c sys_service.h
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) sys_service.c -o ./OBJ/sys_service.o

./OBJ/c2fortran.o: c2fortran.cu
ifeq ($(GPU_CUDA),CUDA)
	$(CUDA_COMP) -ccbin $(CUDA_HOST_COMPILER) $(INC) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) c2fortran.cu -o ./OBJ/c2fortran.o
else
	cp c2fortran.cu c2fortran.c
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) c2fortran.c -o ./OBJ/c2fortran.o
	rm -f c2fortran.c
endif

./OBJ/c2f_ifc.o: c2f_ifc.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) c2f_ifc.F90 -o ./OBJ/c2f_ifc.o

./OBJ/stsubs.o: stsubs.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stsubs.F90 -o ./OBJ/stsubs.o

./OBJ/timer.o: timer.cpp timer.h
	$(CPPCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CPPFLAGS) timer.cpp -o ./OBJ/timer.o

./OBJ/timers.o: timers.F90 ./OBJ/timer.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) timers.F90 -o ./OBJ/timers.o

./OBJ/nvtx_profile.o: nvtx_profile.c nvtx_profile.h
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) nvtx_profile.c -o ./OBJ/nvtx_profile.o

./OBJ/mpi_fort.o: mpi_fort.c
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) mpi_fort.c -o ./OBJ/mpi_fort.o

./OBJ/service_mpi.o: service_mpi.F90 ./OBJ/mpi_fort.o ./OBJ/stsubs.o ./OBJ/c2f_ifc.o ./OBJ/dil_basic.o ./OBJ/timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) service_mpi.F90 -o ./OBJ/service_mpi.o

./OBJ/distributed.o: distributed.F90 ./OBJ/service_mpi.o ./OBJ/stsubs.o ./OBJ/c2f_ifc.o ./OBJ/timers.o ./OBJ/pack_prim.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) distributed.F90 -o ./OBJ/distributed.o

./OBJ/pack_prim.o: pack_prim.F90 ./OBJ/stsubs.o ./OBJ/dil_basic.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) pack_prim.F90 -o ./OBJ/pack_prim.o

./OBJ/test_pack_prim.o: test_pack_prim.F90 ./OBJ/pack_prim.o ./OBJ/dil_basic.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) test_pack_prim.F90 -o ./OBJ/test_pack_prim.o

./OBJ/main.o: main.F90 ./OBJ/dil_basic.o ./OBJ/timers.o ./OBJ/service_mpi.o ./OBJ/distributed.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) main.F90 -o ./OBJ/main.o


.PHONY: clean
clean:
	rm -f *.x *.a ./OBJ/* *.mod *.modmic *.ptx *.log
