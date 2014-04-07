NAME = qforce.v13.01.x
FC = mpif90
CC = gcc
CUDA_C = /usr/local/cuda-5.5/bin/nvcc
MPI_INC = /usr/include/mpich2
CUDA_INC = /usr/local/cuda-5.5/include
CUDA_LIB = /usr/local/cuda-5.5/lib
CUDA_LINK = -lcublas -lcudart
CUDA_FLAGS = --compile -O3 -arch=sm_11 -DDEBUG
LA_LINK = -lblas -llapack
C_FLAGS = -c -O3
FFLAGS = -c -O3 --free-line-length-none -x f95-cpp-input -fopenmp
LFLAGS = -lgomp

OBJS = cuda2fortran.o extern_names.o stsubs.o service.o combinatoric.o tensor_algebra.o qforce.o \
	tensor_dil_omp.o tensor_algebra_gpu_nvidia.o c_process.o c_proc_bufs.o \
	main.o proceed.o \

$(NAME): $(OBJS)
	$(FC) $(OBJS) -I$(MPI_INC) -I$(CUDA_INC) -L$(CUDA_LIB) $(CUDA_LINK) $(LA_LINK) $(LFLAGS) -o $(NAME)

%.o: %.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) $?

tensor_dil_omp.mod: tensor_dil_omp.o
tensor_dil_omp.o: tensor_dil_omp.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) tensor_dil_omp.F90

qforce.mod: qforce.o extern_names.mod service.mod stsubs.mod combinatoric.mod tensor_algebra.mod
qforce.o: qforce.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) qforce.F90

tensor_algebra.mod: tensor_algebra.o stsubs.mod combinatoric.mod
tensor_algebra.o: tensor_algebra.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) tensor_algebra.F90

combinatoric.mod: combinatoric.o
combinatoric.o: combinatoric.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) combinatoric.F90

service.mod: service.o stsubs.mod
service.o: service.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) service.F90

stsubs.mod: stsubs.o
stsubs.o: stsubs.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) stsubs.F90

extern_names.mod: extern_names.o
extern_names.o: extern_names.F90
	$(FC) -I$(MPI_INC) -I$(CUDA_INC) $(FFLAGS) extern_names.F90

c_proc_bufs.o: c_proc_bufs.cu tensor_algebra_gpu_nvidia.h
	$(CUDA_C) $(CUDA_FLAGS) c_proc_bufs.cu

cuda2fortran.o: cuda2fortran.cu
	$(CUDA_C) $(CUDA_FLAGS) -ptx cuda2fortran.cu
	$(CUDA_C) $(CUDA_FLAGS) cuda2fortran.cu

tensor_algebra_gpu_nvidia.o: tensor_algebra_gpu_nvidia.cu tensor_algebra_gpu_nvidia.h
	$(CUDA_C) $(CUDA_FLAGS) -ptx tensor_algebra_gpu_nvidia.cu
	$(CUDA_C) $(CUDA_FLAGS) tensor_algebra_gpu_nvidia.cu

clean:
	rm *.o *.mod *.x *.ptx
