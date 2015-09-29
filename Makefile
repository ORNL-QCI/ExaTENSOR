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
FC_YES = ftn
FC_NO = /usr/bin/mpif90
FC = $(FC_$(WRAP))
CC_YES = cc
CC_NO = mpicc
CC = $(CC_$(WRAP))
CPP_YES = CC
CPP_NO = mpic++
CPP = $(CPP_$(WRAP))
CUDA_C = nvcc
MPI_INC_MPICH = -I/usr/lib/mpich/include
MPI_INC_OPENMPI = -I/usr/local/include
MPI_INC_YES = -I.
MPI_INC_NO = $(MPI_INC_$(MPILIB))
MPI_INC = $(MPI_INC_$(WRAP))
MPI_LINK_MPICH = -L/usr/lib
MPI_LINK_OPENMPI = -L/usr/local/lib
MPI_LINK_YES = -L.
MPI_LINK_NO = $(MPI_LINK_$(MPILIB))
MPI_LINK = $(MPI_LINK_$(WRAP))
CUDA_INC_YES = -I.
CUDA_INC_NO = -I/usr/local/cuda/include
CUDA_INC = $(CUDA_INC_$(WRAP))
CUDA_LINK_YES = -L. -lcudart -lcublas
CUDA_LINK_NO = -L/usr/local/cuda/lib64 -lcudart -lcublas
CUDA_LINK = $(CUDA_LINK_$(WRAP))

CUDA_FLAGS_DEV = --compile -arch=sm_35 -D CUDA_ARCH=350 -g -G -D DEBUG_GPU
CUDA_FLAGS_OPT = --compile -arch=sm_35 -D CUDA_ARCH=350 -O3
CUDA_FLAGS = $(CUDA_FLAGS_$(TYPE))
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
CFLAGS_DEV = -c -D CUDA_ARCH=350 -g
CFLAGS_OPT = -c -D CUDA_ARCH=350 -O3
CFLAGS = $(CFLAGS_$(TYPE)) -D NO_PHI -D NO_AMD
FFLAGS_INTEL_DEV = -c -g -fpp -vec-threshold4 -openmp
FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -openmp
FFLAGS_CRAY_DEV = -c -D CUDA_ARCH=350 -g
FFLAGS_CRAY_OPT = -c -D CUDA_ARCH=350 -O3
FFLAGS_GNU_DEV = -c -fopenmp -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -pg
FFLAGS_GNU_OPT = -c -fopenmp -O3
FFLAGS_PGI_DEV = -c -mp -Mcache_align -Mbounds -Mchkptr -Mstandard -pg
FFLAGS_PGI_OPT = -c -mp -Mcache_align -Mstandard -O3
FFLAGS = $(FFLAGS_$(TOOLKIT)_$(TYPE)) -D NO_PHI -D NO_AMD
LTHREAD_INTEL = -liomp5
LTHREAD_CRAY  = -L.
LTHREAD_GNU   = -lgomp
LTHREAD_PGI   = -lpthread
LTHREAD = $(LTHREAD_$(TOOLKIT))
LFLAGS = $(LTHREAD) $(MPI_LINK) $(LA_LINK) $(CUDA_LINK) -o

OBJS = dil_kinds.o sys_service.o stsubs.o multords.o combinatoric.o symm_index.o timers.o stack.o lists.o dictionary.o \
	extern_names.o tensor_algebra.o tensor_algebra_cpu.o tensor_algebra_cpu_phi.o tensor_dil_omp.o \
	c2fortran.o mem_manager.o tensor_algebra_gpu_nvidia.o talshf.o talshc.o \
	mpi_fort.o service_mpi.o distributed.o subspaces.o virta.o c_process.o exatensor.o \
	qforce.o main.o

$(NAME): $(OBJS)
	$(FC) $(OBJS) $(LFLAGS) $(NAME)

dil_kinds.o: dil_kinds.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dil_kinds.F90

sys_service.o: sys_service.c
	$(CC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) sys_service.c

stsubs.o: stsubs.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stsubs.F90

multords.o: multords.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) multords.F90

combinatoric.o: combinatoric.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) combinatoric.F90

symm_index.o: symm_index.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) symm_index.F90

timers.o: timers.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) timers.F90

stack.o: stack.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stack.F90

lists.o: lists.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) lists.F90

dictionary.o: dictionary.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dictionary.F90

extern_names.o: extern_names.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) extern_names.F90

tensor_algebra.o: tensor_algebra.F90 dil_kinds.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra.F90

tensor_algebra_cpu.o: tensor_algebra_cpu.F90 tensor_algebra.o stsubs.o combinatoric.o timers.o symm_index.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu.F90

tensor_algebra_cpu_phi.o: tensor_algebra_cpu_phi.F90 tensor_algebra_cpu.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu_phi.F90

tensor_dil_omp.o: tensor_dil_omp.F90 timers.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_dil_omp.F90

c2fortran.o: c2fortran.cu
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) c2fortran.cu

mem_manager.o: mem_manager.cu tensor_algebra.h
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) mem_manager.cu

tensor_algebra_gpu_nvidia.o: tensor_algebra_gpu_nvidia.cu tensor_algebra.h
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) -ptx tensor_algebra_gpu_nvidia.cu
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) tensor_algebra_gpu_nvidia.cu

talshf.o: talshf.F90 tensor_algebra_cpu_phi.o tensor_algebra_gpu_nvidia.o mem_manager.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) talshf.F90

talshc.o: talshc.c talsh.h tensor_algebra.h tensor_algebra_cpu_phi.o tensor_algebra_gpu_nvidia.o mem_manager.o
	$(CC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) talshc.c

mpi_fort.o: mpi_fort.c
	$(CC) $(MPI_INC) $(CUDA_INC) $(CFLAGS) mpi_fort.c

service_mpi.o: service_mpi.F90 mpi_fort.o stsubs.o extern_names.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) service_mpi.F90

distributed.o: distributed.F90 service_mpi.o tensor_algebra.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) distributed.F90

subspaces.o: subspaces.F90 dil_kinds.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) subspaces.F90

virta.o: virta.F90 talshf.o talshc.o distributed.o subspaces.o stack.o lists.o dictionary.o multords.o extern_names.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) virta.F90

c_process.o: c_process.F90 virta.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) c_process.F90

exatensor.o: exatensor.F90 c_process.o service_mpi.o virta.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) exatensor.F90

qforce.o: qforce.F90 dil_kinds.o exatensor.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) qforce.F90

main.o: main.F90 exatensor.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) main.F90

#%.o: %.F90
#	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) $?

clean:
	rm *.o *.mod *.modmic *.ptx *.x
