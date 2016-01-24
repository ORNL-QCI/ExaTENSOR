NAME = qforce.x

#Cross-compiling wrappers: [WRAP|NOWRAP]:
export WRAP = NOWRAP
#Compiler: [GNU|PGI|INTEL|CRAY]:
export TOOLKIT = GNU
#Optimization: [DEV|OPT]:
export BUILD_TYPE = DEV
#MPI Library: [MPICH|OPENMPI]:
export MPILIB = MPICH
#BLAS: [ATLAS|MKL|ACML]:
export BLASLIB = ATLAS
#Nvidia GPU via CUDA: [CUDA|NOCUDA]:
export GPU_CUDA = NOCUDA

#Local paths (for unwrapped compilation):
# MPI:
export PATH_MPICH = /usr/local/mpich3.2
export PATH_OPENMPI = /usr/local/openmpi1.10.1
# BLAS:
export PATH_BLAS_ATLAS = /usr/lib
export PATH_BLAS_INTEL = /usr/lib
export PATH_BLAS_ACML = /usr/lib
PATH_BLAS = $(PATH_BLAS_$(BLASLIB))
# CUDA:
export PATH_CUDA=/usr/local/cuda
#DONE.

#=================
#Fortran compiler:
FC_GNU = gfortran
FC_PGI = pgf90
FC_INTEL = ifort
FC_CRAY = ftn
FC_MPICH = $(PATH_MPICH)/bin/mpif90
FC_OPENMPI = $(PATH_OPENMPI)/bin/mpifort
FC_NOWRAP = $(FC_$(MPILIB))
FC_WRAP = ftn
FCOMP = $(FC_$(WRAP))
#C compiler:
CC_GNU = gcc
CC_PGI = pgcc
CC_INTEL = icc
CC_CRAY = cc
CC_MPICH = $(PATH_MPICH)/bin/mpicc
CC_OPENMPI = $(PATH_OPENMPI)/bin/mpicc
CC_NOWRAP = $(CC_$(MPILIB))
CC_WRAP = cc
CCOMP = $(CC_$(WRAP))
#C++ compiler:
CPP_GNU = g++
CPP_PGI = pgc++
CPP_INTEL = icc
CPP_CRAY = CC
CPP_MPICH = $(PATH_MPICH)/bin/mpic++
CPP_OPENMPI = $(PATH_OPENMPI)/bin/mpic++
CPP_NOWRAP = $(CPP_$(MPILIB))
CPP_WRAP = CC
CPPCOMP = $(CPP_$(WRAP))
#CUDA compiler:
CUDA_COMP = nvcc

#COMPILER INCLUDES:
INC_GNU = -I.
INC_PGI = -I.
INC_INTEL = -I.
INC_CRAY = -I.
INC_NOWRAP = $(INC_$(TOOLKIT))
INC_WRAP = -I.
INC = $(INC_$(WRAP))

#COMPILER LIBS:
LIB_GNU = -L.
LIB_PGI = -L.
LIB_INTEL = -L.
LIB_CRAY = -L.
LIB_NOWRAP = $(LIB_$(TOOLKIT))
LIB_WRAP = -L.
LIB = $(LIB_$(WRAP))

#MPI INCLUDES:
MPI_INC_MPICH = -I$(PATH_MPICH)/include
MPI_INC_OPENMPI = -I$(PATH_OPENMPI)/include
MPI_INC_NOWRAP = $(MPI_INC_$(MPILIB))
MPI_INC_WRAP = -I.
MPI_INC = $(MPI_INC_$(WRAP))

#MPI LIBS:
MPI_LINK_MPICH = -L$(PATH_MPICH)/lib
MPI_LINK_OPENMPI = -L$(PATH_OPENMPI)/lib
MPI_LINK_NOWRAP = $(MPI_LINK_$(MPILIB))
MPI_LINK_WRAP = -L.
MPI_LINK = $(MPI_LINK_$(WRAP))

#LINEAR ALGEBRA FLAGS:
LA_LINK_ATLAS = -L$(PATH_BLAS) -lblas -llapack
LA_LINK_MKL = -L$(PATH_BLAS) -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lrt
LA_LINK_ACML = -L$(PATH_BLAS) -lacml_mp
LA_LINK_WRAP = -L.
LA_LINK_NOWRAP = $(LA_LINK_$(BLASLIB))
LA_LINK = $(LA_LINK_$(WRAP))

#CUDA INCLUDES:
CUDA_INC_CUDA = -I$(PATH_CUDA)/include
CUDA_INC_NOCUDA = -I.
CUDA_INC_NOWRAP = $(CUDA_INC_$(GPU_CUDA))
CUDA_INC_WRAP = -I.
CUDA_INC = $(CUDA_INC_$(WRAP))

#CUDA LIBS:
CUDA_LINK_NOWRAP = -L$(PATH_CUDA)/lib64 -lcudart -lcublas
CUDA_LINK_WRAP = -lcudart -lcublas
CUDA_LINK_CUDA = $(CUDA_LINK_$(WRAP))
CUDA_LINK_NOCUDA = -L.
CUDA_LINK = $(CUDA_LINK_$(GPU_CUDA))

#CUDA FLAGS:
CUDA_HOST_NOWRAP = --compiler-bindir /usr/bin
CUDA_HOST_WRAP = -I.
CUDA_HOST = $(CUDA_HOST_$(WRAP))
CUDA_FLAGS_DEV = --compile -arch=sm_35 -g -G -D CUDA_ARCH=350 -D DEBUG_GPU
CUDA_FLAGS_OPT = --compile -arch=sm_35 -O3 -D CUDA_ARCH=350
CUDA_FLAGS_CUDA = $(CUDA_HOST) $(CUDA_FLAGS_$(BUILD_TYPE))
CUDA_FLAGS_NOCUDA = -I.
CUDA_FLAGS = $(CUDA_FLAGS_$(GPU_CUDA))

#Accelerator support:
NO_ACCEL_CUDA = -D NO_AMD -D NO_PHI -D CUDA_ARCH=350
NO_ACCEL_NOCUDA = -D NO_AMD -D NO_PHI -D NO_GPU
NO_ACCEL = $(NO_ACCEL_$(GPU_CUDA))

#C FLAGS:
CFLAGS_DEV = -c -g $(NO_ACCEL)
CFLAGS_OPT = -c -O3 $(NO_ACCEL)
CFLAGS = $(CFLAGS_$(BUILD_TYPE))

#FORTRAN FLAGS:
FFLAGS_INTEL_DEV = -c -g -fpp -vec-threshold4 -openmp $(NO_ACCEL)
FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -openmp $(NO_ACCEL)
FFLAGS_CRAY_DEV = -c -g $(NO_ACCEL)
FFLAGS_CRAY_OPT = -c -O3 $(NO_ACCEL)
FFLAGS_GNU_DEV = -c -fopenmp -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -pg $(NO_ACCEL)
FFLAGS_GNU_OPT = -c -fopenmp -O3 $(NO_ACCEL)
FFLAGS_PGI_DEV = -c -mp -Mcache_align -Mbounds -Mchkptr -Mstandard -pg $(NO_ACCEL)
FFLAGS_PGI_OPT = -c -mp -Mcache_align -Mstandard -O3 $(NO_ACCEL)
FFLAGS = $(FFLAGS_$(TOOLKIT)_$(BUILD_TYPE))

#THREADS:
LTHREAD_GNU   = -lgomp
LTHREAD_PGI   = -lpthread
LTHREAD_INTEL = -liomp5
LTHREAD_CRAY  = -L.
LTHREAD = $(LTHREAD_$(TOOLKIT))

#LINKING:
LFLAGS = $(LIB) $(LTHREAD) $(MPI_LINK) $(LA_LINK) $(CUDA_LINK) -o $(NAME)

OBJS =  dil_basic.o sys_service.o stsubs.o multords.o combinatoric.o symm_index.o timers.o stack.o lists.o dictionary.o \
	c2f_ifc.o tensor_algebra.o tensor_algebra_cpu.o tensor_algebra_cpu_phi.o tensor_dil_omp.o \
	c2fortran.o mem_manager.o tensor_algebra_gpu_nvidia.o talshf.o talshc.o \
	mpi_fort.o service_mpi.o distributed.o subspaces.o virta.o c_process.o exatensor.o \
	qforce.o main.o

$(NAME): $(OBJS)
	$(FCOMP) $(OBJS) $(LFLAGS)

dil_basic.o: dil_basic.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dil_basic.F90

sys_service.o: sys_service.c
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) sys_service.c

stsubs.o: stsubs.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stsubs.F90

multords.o: multords.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) multords.F90

combinatoric.o: combinatoric.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) combinatoric.F90

symm_index.o: symm_index.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) symm_index.F90

timers.o: timers.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) timers.F90

stack.o: stack.F90 timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stack.F90

lists.o: lists.F90 timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) lists.F90

dictionary.o: dictionary.F90 timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dictionary.F90

c2f_ifc.o: c2f_ifc.F90
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) c2f_ifc.F90

tensor_algebra.o: tensor_algebra.F90 dil_basic.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra.F90

tensor_algebra_cpu.o: tensor_algebra_cpu.F90 tensor_algebra.o stsubs.o combinatoric.o symm_index.o timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu.F90

tensor_algebra_cpu_phi.o: tensor_algebra_cpu_phi.F90 tensor_algebra_cpu.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu_phi.F90

tensor_dil_omp.o: tensor_dil_omp.F90 timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_dil_omp.F90

c2fortran.o: c2fortran.c
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) c2fortran.c

mem_manager.o: mem_manager.c tensor_algebra.h
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) mem_manager.c

tensor_algebra_gpu_nvidia.o: tensor_algebra_gpu_nvidia.cu tensor_algebra.h
	$(CUDA_COMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) -ptx tensor_algebra_gpu_nvidia.cu
	$(CUDA_COMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) tensor_algebra_gpu_nvidia.cu

talshf.o: talshf.F90 tensor_algebra_cpu_phi.o tensor_algebra_gpu_nvidia.o mem_manager.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) talshf.F90

talshc.o: talshc.c talsh.h tensor_algebra.h tensor_algebra_cpu_phi.o tensor_algebra_gpu_nvidia.o mem_manager.o
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) talshc.c

mpi_fort.o: mpi_fort.c
	$(CCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) mpi_fort.c

service_mpi.o: service_mpi.F90 mpi_fort.o stsubs.o c2f_ifc.o dil_basic.o timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) service_mpi.F90

distributed.o: distributed.F90 service_mpi.o stsubs.o c2f_ifc.o timers.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) distributed.F90

subspaces.o: subspaces.F90 dil_basic.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) subspaces.F90

virta.o: virta.F90 talshf.o talshc.o distributed.o subspaces.o stack.o lists.o dictionary.o multords.o c2f_ifc.o service_mpi.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) virta.F90

c_process.o: c_process.F90 virta.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) c_process.F90

exatensor.o: exatensor.F90 virta.o c_process.o service_mpi.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) exatensor.F90

qforce.o: qforce.F90 exatensor.o dil_basic.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) qforce.F90

main.o: main.F90 exatensor.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) main.F90

#%.o: %.F90
#	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) $?

.PHONY: clean
clean:
	rm *.x *.o *.a *.mod *.modmic *.ptx
