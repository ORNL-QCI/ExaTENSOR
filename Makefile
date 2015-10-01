NAME = qforce.v13.01.x
#Cross-compiling wrappers: [YES|NO]:
WRAP = NO
#Compiler: [GNU|PGI|INTEL|CRAY]:
TOOLKIT = GNU
#Optimization: [DEV|OPT]:
TYPE = OPT
#MPI Library: [MPICH|OPENMPI]:
MPILIB = MPICH
#DONE.
#-----------------------------------------------
#Fortran compiler:
FC_GNU = /home/dima/build/gcc_5_2_0/bin/gfortran
FC_PGI = pgf90
FC_INTEL = ifort
FC_CRAY = ftn
FC_MPICH = /usr/bin/mpif90
FC_OPENMPI = /usr/local/bin/mpifort
FC_NO = $(FC_$(MPILIB))
FC_YES = ftn
FC = $(FC_$(WRAP))
#C compiler:
CC_GNU = /home/dima/build/gcc_5_2_0/bin/gcc
CC_PGI = pgcc
CC_INTEL = icc
CC_CRAY = cc
CC_MPICH = /usr/bin/mpicc
CC_OPENMPI = /usr/local/bin/mpicc
CC_NO = $(CC_$(MPILIB))
CC_YES = cc
CC = $(CC_$(WRAP))
#C++ compiler:
CPP_GNU = /home/dima/build/gcc_5_2_0/bin/g++
CPP_PGI = pgc++
CPP_INTEL = icc
CPP_CRAY = CC
CPP_MPICH = /usr/bin/mpic++
CPP_OPENMPI = /usr/local/bin/mpic++
CPP_NO = $(CPP_$(MPILIB))
CPP_YES = CC
CPP = $(CPP_$(WRAP))
#CUDA compiler:
CUDA_C = nvcc

#COMPILER INCLUDES:
INC_GNU = -I/home/dima/build/gcc_5_2_0/include/c++/5.2.0
INC_PGI = -I.
INC_INTEL = -I.
INC_CRAY = -I.
INC_NO = $(INC_$(TOOLKIT))
INC_YES = -I.
INC = $(INC_$(WRAP))

#MPI INCLUDES:
MPI_INC_MPICH = -I/usr/lib/mpich/include
MPI_INC_OPENMPI = -I/usr/local/include
MPI_INC_NO = $(MPI_INC_$(MPILIB))
MPI_INC_YES = -I.
MPI_INC = $(MPI_INC_$(WRAP))

#CUDA INCLUDES:
CUDA_INC = -I/usr/local/cuda/include

#COMPILER LIBS:
LIB_GNU = -L/home/dima/build/gcc_5_2_0/lib64
LIB_PGI = -L.
LIB_INTEL = -L.
LIB_CRAY = -L.
LIB_NO = $(LIB_$(TOOLKIT))
LIB_YES = -L.
LIB = $(LIB_$(WRAP))

#MPI LIBS:
MPI_LINK_MPICH = -L/usr/lib
MPI_LINK_OPENMPI = -L/usr/local/lib
MPI_LINK_NO = $(MPI_LINK_$(MPILIB))
MPI_LINK_YES = -L.
MPI_LINK = $(MPI_LINK_$(WRAP))

#CUDA LIBS:
CUDA_LINK_NO = -L/usr/local/cuda/lib64 -lcudart -lcublas
CUDA_LINK_YES = -L. -lcudart -lcublas
CUDA_LINK = $(CUDA_LINK_$(WRAP))

#CUDA FLAGS:
CUDA_HOST_NO = --compiler-bindir /home/dima/build/gcc_5_2_0/bin
CUDA_HOST_YES = " "
CUDA_HOST = $(CUDA_HOST_$(WRAP))
CUDA_FLAGS_DEV = --compile -arch=sm_35 -D CUDA_ARCH=350 -g -G -D DEBUG_GPU
CUDA_FLAGS_OPT = --compile -arch=sm_35 -D CUDA_ARCH=350 -O3
CUDA_FLAGS = $(CUDA_HOST) $(CUDA_FLAGS_$(TYPE))
#LINEAR ALGEBRA FLAGS:
LA_LINK_MKL = -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lrt
LA_LINK_ACML = -lacml_mp -L/opt/acml/5.3.1/gfortran64_fma4_mp/lib
LA_LINK_DEFAULT_YES = -L.
LA_LINK_DEFAULT_NO = -L/usr/lib/atlas-base/atlas -lblas -llapack
LA_LINK_DEFAULT = $(LA_LINK_DEFAULT_$(WRAP))
LA_LINK_INTEL = $(LA_LINK_DEFAULT)
LA_LINK_CRAY = $(LA_LINK_DEFAULT)
LA_LINK_GNU = $(LA_LINK_DEFAULT)
LA_LINK_PGI = $(LA_LINK_DEFAULT)
LA_LINK = $(LA_LINK_$(TOOLKIT))
#C FLAGS:
CFLAGS_DEV = -c -D CUDA_ARCH=350 -g
CFLAGS_OPT = -c -D CUDA_ARCH=350 -O3
CFLAGS = $(CFLAGS_$(TYPE)) -D NO_PHI -D NO_AMD
#FORTRAN FLAGS:
FFLAGS_INTEL_DEV = -c -g -fpp -vec-threshold4 -openmp
FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -openmp
FFLAGS_CRAY_DEV = -c -D CUDA_ARCH=350 -g
FFLAGS_CRAY_OPT = -c -D CUDA_ARCH=350 -O3
FFLAGS_GNU_DEV = -c -fopenmp -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -pg
FFLAGS_GNU_OPT = -c -fopenmp -O3
FFLAGS_PGI_DEV = -c -mp -Mcache_align -Mbounds -Mchkptr -Mstandard -pg
FFLAGS_PGI_OPT = -c -mp -Mcache_align -Mstandard -O3
FFLAGS = $(FFLAGS_$(TOOLKIT)_$(TYPE)) -D NO_PHI -D NO_AMD
#THREADS:
LTHREAD_INTEL = -liomp5
LTHREAD_CRAY  = -L.
LTHREAD_GNU   = -lgomp
LTHREAD_PGI   = -lpthread
LTHREAD = $(LTHREAD_$(TOOLKIT))
LFLAGS = $(LIB) $(LTHREAD) $(MPI_LINK) $(LA_LINK) $(CUDA_LINK) -o

OBJS = dil_kinds.o sys_service.o stsubs.o multords.o combinatoric.o symm_index.o timers.o stack.o lists.o dictionary.o \
	extern_names.o tensor_algebra.o tensor_algebra_cpu.o tensor_algebra_cpu_phi.o tensor_dil_omp.o \
	c2fortran.o mem_manager.o tensor_algebra_gpu_nvidia.o talshf.o talshc.o \
	mpi_fort.o service_mpi.o distributed.o subspaces.o virta.o c_process.o exatensor.o \
	qforce.o main.o

$(NAME): $(OBJS)
	$(FC) $(OBJS) $(LFLAGS) $(NAME)

dil_kinds.o: dil_kinds.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dil_kinds.F90

sys_service.o: sys_service.c
	$(CC) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) sys_service.c

stsubs.o: stsubs.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stsubs.F90

multords.o: multords.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) multords.F90

combinatoric.o: combinatoric.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) combinatoric.F90

symm_index.o: symm_index.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) symm_index.F90

timers.o: timers.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) timers.F90

stack.o: stack.F90 timers.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stack.F90

lists.o: lists.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) lists.F90

dictionary.o: dictionary.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dictionary.F90

extern_names.o: extern_names.F90
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) extern_names.F90

tensor_algebra.o: tensor_algebra.F90 dil_kinds.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra.F90

tensor_algebra_cpu.o: tensor_algebra_cpu.F90 tensor_algebra.o stsubs.o combinatoric.o timers.o symm_index.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu.F90

tensor_algebra_cpu_phi.o: tensor_algebra_cpu_phi.F90 tensor_algebra_cpu.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu_phi.F90

tensor_dil_omp.o: tensor_dil_omp.F90 timers.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_dil_omp.F90

c2fortran.o: c2fortran.cu
	$(CUDA_C) $(INC) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) c2fortran.cu

mem_manager.o: mem_manager.cu tensor_algebra.h
	$(CUDA_C) $(INC) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) mem_manager.cu

tensor_algebra_gpu_nvidia.o: tensor_algebra_gpu_nvidia.cu tensor_algebra.h
	$(CUDA_C) $(INC) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) -ptx tensor_algebra_gpu_nvidia.cu
	$(CUDA_C) $(INC) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) tensor_algebra_gpu_nvidia.cu

talshf.o: talshf.F90 tensor_algebra_cpu_phi.o tensor_algebra_gpu_nvidia.o mem_manager.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) talshf.F90

talshc.o: talshc.c talsh.h tensor_algebra.h tensor_algebra_cpu_phi.o tensor_algebra_gpu_nvidia.o mem_manager.o
	$(CC) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) talshc.c

mpi_fort.o: mpi_fort.c
	$(CC) $(INC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) mpi_fort.c

service_mpi.o: service_mpi.F90 mpi_fort.o stsubs.o extern_names.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) service_mpi.F90

distributed.o: distributed.F90 service_mpi.o tensor_algebra.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) distributed.F90

subspaces.o: subspaces.F90 dil_kinds.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) subspaces.F90

virta.o: virta.F90 talshf.o talshc.o distributed.o subspaces.o stack.o lists.o dictionary.o multords.o extern_names.o service_mpi.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) virta.F90

c_process.o: c_process.F90 virta.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) c_process.F90

exatensor.o: exatensor.F90 c_process.o service_mpi.o virta.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) exatensor.F90

qforce.o: qforce.F90 dil_kinds.o exatensor.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) qforce.F90

main.o: main.F90 exatensor.o
	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) main.F90

#%.o: %.F90
#	$(FC) $(INC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) $?

clean:
	rm *.o *.mod *.modmic *.ptx *.x
