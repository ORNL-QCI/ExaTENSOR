NAME = intravirt

#ADJUST THE FOLLOWING ACCORDINGLY:
#Cross-compiling wrappers: [WRAP|NOWRAP]:
export WRAP ?= NOWRAP
#Compiler: [GNU|PGI|INTEL|CRAY|IBM]:
export TOOLKIT ?= GNU
#Optimization: [DEV|OPT]:
export BUILD_TYPE ?= OPT
#MPI Library: [MPICH|OPENMPI]:
export MPILIB ?= MPICH
#BLAS: [ATLAS|MKL|ACML|ESSL|NONE]:
export BLASLIB ?= ATLAS
#Nvidia GPU via CUDA: [CUDA|NOCUDA]:
export GPU_CUDA ?= NOCUDA
#Nvidia GPU architecture (two digits):
export GPU_SM_ARCH ?= 35
#Operating system: [LINUX|NO_LINUX]:
export EXA_OS ?= LINUX

#ADJUST EXTRAS (optional):
#Fast GPU tensor transpose (cuTT library): [YES|NO]:
export WITH_CUTT ?= NO

#WORKAROUNDS (ignore if you do not experience problems):
#Fool CUDA 7.0 with GCC > 4.9: [YES|NO]:
export FOOL_CUDA ?= NO

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

# CUDA (only if you build with CUDA):
export PATH_CUDA ?= /usr/local/cuda
#  Only reset these if CUDA files are spread in the system directories:
 export PATH_CUDA_INC ?= $(PATH_CUDA)/include
 export PATH_CUDA_LIB ?= $(PATH_CUDA)/lib64
 export PATH_CUDA_BIN ?= $(PATH_CUDA)/bin
# cuTT path (if you use cuTT library):
export PATH_CUTT ?= /home/div/src/cutt

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
LA_LINK_MKL = -L$(PATH_BLAS_MKL) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lm -ldl
else
LA_LINK_MKL = -L$(PATH_BLAS_MKL) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl -L$(PATH_BLAS_MKL_DEP) -liomp5
endif
LA_LINK_ACML = -L$(PATH_BLAS_ACML) -lacml_mp
LA_LINK_ESSL = -L$(PATH_BLAS_ESSL) -lesslsmp -L$(PATH_BLAS_ESSL_DEP) -lxlf90_r -lxlfmath
ifeq ($(BLASLIB),NONE)
LA_LINK_NOWRAP = -L.
else
LA_LINK_NOWRAP = $(LA_LINK_$(BLASLIB))
endif
ifeq ($(BLASLIB),MKL)
LA_LINK_WRAP = $(LA_LINK_MKL)
else
LA_LINK_WRAP = -L.
endif
LA_LINK = $(LA_LINK_$(WRAP))

#CUDA INCLUDES:
ifeq ($(GPU_CUDA),CUDA)
CUDA_INC_NOWRAP = -I$(PATH_CUDA_INC)
CUDA_INC_WRAP = -I.
ifeq ($(WITH_CUTT),YES)
CUDA_INC = $(CUDA_INC_$(WRAP)) -I$(PATH_CUTT)/include
else
CUDA_INC = $(CUDA_INC_$(WRAP))
endif
else
CUDA_INC = -I.
endif

#CUDA LIBS:
CUDA_LINK_NOWRAP = -L$(PATH_CUDA_LIB) -lcudart -lcublas
CUDA_LINK_WRAP = -lcudart -lcublas
CUDA_LINK_CUDA = $(CUDA_LINK_$(WRAP))
CUDA_LINK_NOCUDA = -L.
CUDA_LINK = $(CUDA_LINK_$(GPU_CUDA))

#CUDA FLAGS:
ifeq ($(GPU_CUDA),CUDA)
GPU_SM = sm_$(GPU_SM_ARCH)
GPU_ARCH = $(GPU_SM_ARCH)0
CUDA_HOST_NOWRAP = --compiler-bindir /usr/bin
CUDA_HOST_WRAP = -I.
CUDA_HOST = $(CUDA_HOST_$(WRAP))
CUDA_FLAGS_DEV = --compile -arch=$(GPU_SM) -g -G -lineinfo -DDEBUG_GPU
CUDA_FLAGS_OPT = --compile -arch=$(GPU_SM) -O3 -lineinfo
CUDA_FLAGS_CUDA = $(CUDA_HOST) $(CUDA_FLAGS_$(BUILD_TYPE)) -D_FORCE_INLINES
ifeq ($(FOOL_CUDA),NO)
CUDA_FLAGS_PRE1 = $(CUDA_FLAGS_CUDA) -D$(EXA_OS)
else
CUDA_FLAGS_PRE1 = $(CUDA_FLAGS_CUDA) -D$(EXA_OS) -D__GNUC__=4
endif
ifeq ($(WITH_CUTT),YES)
CUDA_FLAGS_PRE2 = $(CUDA_FLAGS_PRE1) -DUSE_CUTT
else
CUDA_FLAGS_PRE2 = $(CUDA_FLAGS_PRE1)
endif
ifeq ($(GPU_FINE_TIMING),YES)
CUDA_FLAGS = $(CUDA_FLAGS_PRE2) -DGPU_FINE_TIMING
else
CUDA_FLAGS = $(CUDA_FLAGS_PRE2)
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
ifeq ($(TOOLKIT),GNU)
NO_GNU :=
else
NO_GNU := -DNO_GNU
endif
ifeq ($(BLASLIB),NONE)
NO_BLAS = -DNO_BLAS
else
NO_BLAS :=
endif
ifeq ($(GPU_CUDA),CUDA)
NO_GPU = -DCUDA_ARCH=$(GPU_ARCH)
else
NO_GPU = -DNO_GPU
endif
NO_AMD = -DNO_AMD
NO_PHI = -DNO_PHI

#C FLAGS:
CFLAGS_DEV = -c -g -O0 -D_DEBUG
CFLAGS_OPT = -c -O3
CFLAGS = $(CFLAGS_$(BUILD_TYPE)) $(NO_GPU) $(NO_AMD) $(NO_PHI) $(NO_BLAS) -D$(EXA_OS) -fPIC

#FORTRAN FLAGS:
FFLAGS_INTEL_DEV = -c -g -O0 -fpp -vec-threshold4 -qopenmp -mkl=parallel
#FFLAGS_INTEL_DEV = -c -g -fpp -vec-threshold4 -openmp
FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -qopenmp -mkl=parallel
#FFLAGS_INTEL_OPT = -c -O3 -fpp -vec-threshold4 -openmp
FFLAGS_CRAY_DEV = -c -g
FFLAGS_CRAY_OPT = -c -O3
FFLAGS_GNU_DEV = -c -fopenmp -fbacktrace -fcheck=bounds -fcheck=array-temps -fcheck=pointer -g -O0
FFLAGS_GNU_OPT = -c -fopenmp -O3
FFLAGS_PGI_DEV = -c -Mcache_align -Mbounds -Mchkptr -Mstandard -g -O0
FFLAGS_PGI_OPT = -c -Mcache_align -Mstandard -O3
FFLAGS_IBM_DEV = -c -qsmp=omp -g -qkeepparm -qcheck -qsigtrap
FFLAGS_IBM_OPT = -c -qsmp=omp -O3
FFLAGS = $(FFLAGS_$(TOOLKIT)_$(BUILD_TYPE)) $(DF)$(NO_GPU) $(DF)$(NO_AMD) $(DF)$(NO_PHI) $(DF)$(NO_BLAS) $(DF)-D$(EXA_OS) $(DF)$(NO_GNU) -fPIC

#ExaTensor includes:
MY_INC1 = -I../UTILITY -I../GFC -I../DDSS -I../TALSH -I../DSVP
MY_INC2 = -I../UTILITY/OBJ -I../GFC/OBJ -I../DDSS/OBJ -I../TALSH/OBJ -I../DSVP/OBJ
MY_INC = $(MY_INC1) $(MY_INC2)
MY_LIB = ../DSVP/libdsvp.a ../TALSH/libtalsh.a ../DDSS/libddss.a ../GFC/libgfc.a ../UTILITY/libutility.a

#THREADS:
LTHREAD_GNU   = -lgomp
LTHREAD_PGI   = -mp -lpthread
LTHREAD_INTEL = -liomp5
LTHREAD_CRAY  = -L.
LTHREAD_IBM   = -lxlsmp
LTHREAD = $(LTHREAD_$(TOOLKIT))

#LINKING:
LFLAGS = $(LTHREAD) $(MPI_LINK) $(LA_LINK) $(CUDA_LINK) $(LIB)

OBJS = ./OBJ/subspaces.o ./OBJ/hardware.o ./OBJ/tensor_recursive.o ./OBJ/virta.o

$(NAME): $(OBJS) ./OBJ/test_intra.o $(MY_LIB)
	ar cr lib$(NAME).a $(OBJS) $(MY_LIB)
	$(FCOMP) ./OBJ/test_intra.o lib$(NAME).a $(MY_LIB) $(LFLAGS) -o test_$(NAME).x

./OBJ/subspaces.o: subspaces.F90 ../UTILITY/OBJ/dil_basic.o ../UTILITY/OBJ/stsubs.o ../GFC/OBJ/gfc_base.o ../GFC/OBJ/gfc_list.o ../GFC/OBJ/gfc_vec_tree.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(MY_INC) $(FFLAGS) subspaces.F90 -o ./OBJ/subspaces.o

./OBJ/hardware.o: hardware.F90 ../UTILITY/OBJ/dil_basic.o ../UTILITY/OBJ/stsubs.o ../UTILITY/OBJ/parse_prim.o ../GFC/OBJ/gfc_base.o ../GFC/OBJ/gfc_vec_tree.o ./OBJ/subspaces.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(MY_INC) $(FFLAGS) hardware.F90 -o ./OBJ/hardware.o

./OBJ/tensor_recursive.o: tensor_recursive.F90 ../TALSH/OBJ/tensor_algebra.o ../UTILITY/OBJ/stsubs.o ../UTILITY/OBJ/timers.o ../UTILITY/OBJ/combinatoric.o ../GFC/OBJ/gfc_base.o ../GFC/OBJ/gfc_list.o ../GFC/OBJ/gfc_vector.o ../GFC/OBJ/gfc_vec_tree.o ../GFC/OBJ/gfc_dictionary.o ../GFC/OBJ/multords.o ./OBJ/subspaces.o ../DDSS/OBJ/pack_prim.o ../DDSS/OBJ/distributed.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(MY_INC) $(FFLAGS) tensor_recursive.F90 -o ./OBJ/tensor_recursive.o

./OBJ/virta.o: virta.F90 ../UTILITY/OBJ/dil_basic.o ../GFC/OBJ/gfc_base.o ../GFC/OBJ/gfc_dictionary.o ../TALSH/OBJ/talshf.o ../DDSS/OBJ/service_mpi.o ../DDSS/OBJ/pack_prim.o ../DDSS/OBJ/distributed.o ../DSVP/OBJ/dsvp_base.o ./OBJ/hardware.o ./OBJ/subspaces.o ./OBJ/tensor_recursive.o
	mkdir -p ./OBJ
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(MY_INC) $(FFLAGS) virta.F90 -o ./OBJ/virta.o

./OBJ/test_intra.o: test_intra.F90 $(OBJS) ../TALSH/OBJ/tensor_algebra.o ../UTILITY/OBJ/parse_prim.o
	$(FCOMP) $(INC) $(MPI_INC) $(CUDA_INC) $(MY_INC) $(FFLAGS) test_intra.F90 -o ./OBJ/test_intra.o


.PHONY: clean
clean:
	rm -f *.x *.a ./OBJ/* *.mod *.modmic *.ptx *.log
