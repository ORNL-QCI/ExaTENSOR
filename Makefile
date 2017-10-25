NAME = ExaTensor

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
export GPU_SM_ARCH ?= 30
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

$(NAME):
	$(MAKE) -C ./UTILITY
	$(MAKE) -C ./GFC
	$(MAKE) -C ./DDSS
	$(MAKE) -C ./TALSH
	$(MAKE) -C ./DSVP
	$(MAKE) -C ./INTRAVIRT
	$(MAKE) -C ./INTERVIRT
	$(MAKE) -C ./TN
	$(MAKE) -C ./QFORCE
	rm -f ./include/*
	rm -f ./lib/*
	rm -f ./bin/*
ifeq ($(TOOLKIT),CRAY)
	cp ./INTERVIRT/OBJ/EXATENSOR.mod ./
	cp ./INTERVIRT/OBJ/EXATENSOR.mod ./include/
	cp ./INTRAVIRT/OBJ/TENSOR_RECURSIVE.mod ./include/
	cp ./INTRAVIRT/OBJ/SUBSPACES.mod ./include/
	cp ./TALSH/OBJ/TALSH.mod ./include/
	cp ./TALSH/OBJ/TENSOR_ALGEBRA.mod ./include/
	cp ./TALSH/OBJ/DIL_BASIC.mod ./include/
else
	cp ./INTERVIRT/exatensor.mod ./
	cp ./INTERVIRT/exatensor.mod ./include/
	cp ./INTRAVIRT/tensor_recursive.mod ./include/
	cp ./INTRAVIRT/subspaces.mod ./include/
	cp ./TALSH/talsh.mod ./include/
	cp ./TALSH/tensor_algebra.mod ./include/
	cp ./TALSH/dil_basic.mod ./include/
endif
	cp ./[A-Z]*/*.h ./include/
	cp ./[A-Z]*/*.hpp ./include/
	cp ./TN/*.cpp ./include/
	cp ./[A-Z]*/*.a ./lib/
	cp ./[A-Z]*/*.x ./bin/
	cp ./INTERVIRT/libExaTensor.a ./
	cp ./QFORCE/Qforce.x ./
	echo "Finished successfully!"

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.mod ./*/*.x ./*/*.a ./*/*.mod ./*/OBJ/* ./bin/* ./lib/* ./include/*
