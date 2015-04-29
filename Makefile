NAME = qforce.v13.01.x
FC  = ftn
CC  = cc
CPP = CC
CUDA_C = nvcc
MPI_INC = -I.
CUDA_INC = -I.
CUDA_LINK = -lcudart -lcublas -L.
CUDA_FLAGS_DEV = --compile -arch=sm_35 -D CUDA_ARCH=350 -g -G -D DEBUG_GPU
CUDA_FLAGS_OPT = --compile -arch=sm_35 -D CUDA_ARCH=350 -O3
CUDA_FLAGS = $(CUDA_FLAGS_OPT)
LA_LINK_INTEL = -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lrt
LA_LINK_AMD = -lacml_mp -L/opt/acml/5.3.1/gfortran64_mp/lib
LA_LINK_CRAY = " "
LA_LINK = $(LA_LINK_CRAY)
CFLAGS_DEV = -c -D CUDA_ARCH=350 -g
CFLAGS_OPT = -c -D CUDA_ARCH=350 -O3
CFLAGS = $(CFLAGS_OPT)
FFLAGS_DEV = -c -D CUDA_ARCH=350 -g
FFLAGS_OPT = -c -D CUDA_ARCH=350 -O3
FFLAGS_DEV_GNU = $(FFLAGS_DEV) -fopenmp -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -pg
FFLAGS_OPT_GNU = $(FFLAGS_OPT) -fopenmp
FFLAGS_DEV_PGI = $(FFLAGS_DEV) -mp -Mcache_align -Mbounds -Mchkptr
FFLAGS_OPT_PGI = $(FFLAGS_OPT) -mp -Mcache_align
FFLAGS_DEV_INTEL = $(FFLAGS_DEV) -fpp -vec-threshold4 -vec-report2 -openmp -openmp-report2 -D USE_MKL
FFLAGS_OPT_INTEL = $(FFLAGS_OPT) -fpp -vec-threshold4 -vec-report2 -openmp -openmp-report2 -D USE_MKL
FFLAGS = $(FFLAGS_OPT) -D NO_PHI -D NO_AMD
LFLAGS_GNU = -lgomp
LFLAGS_PGI = -lpthread
LFLAGS = $(LA_LINK) $(CUDA_LINK) -o

OBJS =  stsubs.o combinatoric.o timers.o extern_names.o lists.o dictionary.o symm_index.o \
	tensor_algebra.o tensor_algebra_cpu.o tensor_algebra_cpu_phi.o tensor_dil_omp.o \
	service_mpi.o cuda2fortran.o c_proc_bufs.o tensor_algebra_gpu_nvidia.o sys_service.o \
	distributed.o subspaces.o exatensor.o c_process.o qforce.o proceed.o main.o

$(NAME): $(OBJS)
	$(FC) $(OBJS) $(LFLAGS) $(NAME)

stsubs.o: stsubs.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) stsubs.F90

combinatoric.o: combinatoric.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) combinatoric.F90

timers.o: timers.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) timers.F90

extern_names.o: extern_names.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) extern_names.F90

lists.o: lists.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) lists.F90

dictionary.o: dictionary.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) dictionary.F90

symm_index.o: symm_index.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) symm_index.F90

tensor_algebra.o: tensor_algebra.F90
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra.F90

tensor_algebra_cpu.o: tensor_algebra_cpu.F90 tensor_algebra.o stsubs.o combinatoric.o timers.o symm_index.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu.F90

tensor_algebra_cpu_phi.o: tensor_algebra_cpu_phi.F90 tensor_algebra_cpu.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_algebra_cpu_phi.F90

tensor_dil_omp.o: tensor_dil_omp.F90 timers.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) tensor_dil_omp.F90

service_mpi.o: service_mpi.F90 stsubs.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) service_mpi.F90

cuda2fortran.o: cuda2fortran.cu
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) cuda2fortran.cu

c_proc_bufs.o: c_proc_bufs.cu tensor_algebra.h
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) c_proc_bufs.cu

tensor_algebra_gpu_nvidia.o: tensor_algebra_gpu_nvidia.cu tensor_algebra.h
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) -ptx tensor_algebra_gpu_nvidia.cu
	$(CUDA_C) $(MPI_INC) $(CUDA_INC) $(CUDA_FLAGS) tensor_algebra_gpu_nvidia.cu

sys_service.o: sys_service.c
	$(CC) $(CFLAGS) sys_service.c

distributed.o: distributed.F90 service_mpi.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) distributed.F90

subspaces.o: subspaces.F90 tensor_algebra.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) subspaces.F90

exatensor.o: exatensor.F90 tensor_algebra_cpu_phi.o distributed.o subspaces.o lists.o dictionary.o extern_names.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) exatensor.F90

c_process.o: c_process.F90 exatensor.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) c_process.F90

qforce.o: qforce.F90 exatensor.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) qforce.F90

proceed.o: proceed.F90 exatensor.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) proceed.F90

main.o: main.F90 exatensor.o qforce.o
	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) main.F90

#%.o: %.F90
#	$(FC) $(MPI_INC) $(CUDA_INC) $(FFLAGS) $?

clean:
	rm *.o *.mod *.modmic *.ptx *.x
