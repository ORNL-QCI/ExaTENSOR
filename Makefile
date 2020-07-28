NAME = ExaTensor

#ADJUST THE FOLLOWING ENVIRONMENT VARIABLES ACCORDINGLY (choices are given)
#until you see "YOU ARE DONE!". The comments will guide you through (read them).
#Alternatively you can export all relevant environment variables such that this
#Makefile will pick their values, so you will not need to update anything here.
#However you will still need to read the meaning of those variables below.

#Cray cross-compiling wrappers (only for Cray): [WRAP|NOWRAP]:
export WRAP ?= NOWRAP
#Compiler: [GNU|INTEL|CRAY|IBM|PGI]:
export TOOLKIT ?= GNU
#Optimization: [DEV|OPT|PRF]:
export BUILD_TYPE ?= OPT
#MPI library base: [MPICH|OPENMPI]:
export MPILIB ?= MPICH
#BLAS: [ATLAS|MKL|OPENBLAS|ACML|LIBSCI|ESSL|NONE]:
export BLASLIB ?= NONE
#NVIDIA GPU via CUDA: [CUDA|NOCUDA]:
export GPU_CUDA ?= NOCUDA
#NVIDIA GPU architecture (two digits, >=35):
export GPU_SM_ARCH ?= 35
#Operating system: [LINUX|NO_LINUX]:
export EXA_OS ?= LINUX
#Only for Linux DEV builds with GNU: [YES|NO]:
export LINUX_GNU_ASAN ?= NO


#ADJUST EXTRAS (optional):

#LAPACK: [YES|NO]:
export WITH_LAPACK ?= NO
#Fast GPU tensor transpose (cuTT library): [YES|NO]:
export WITH_CUTT ?= NO
#In-place GPU tensor contraction (cuTensor library): [YES|NO]:
export WITH_CUTENSOR ?= NO


#ADJUST PARTIAL BUILD OPTIONS:

#Disable actual build completely (debug): [YES|NO]:
export EXA_NO_BUILD ?= NO
#Only enable TAL-SH build ($EXA_NO_BUILD must be NO): [YES|NO]:
export EXA_TALSH_ONLY ?= NO
#The build is part of the ExaTN library build: [YES|NO]:
export EXATN_SERVICE ?= NO


#SET YOUR LOCAL PATHS (for direct builds without Cray compiler wrappers):

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

#BLAS library (whichever you have chosen above):
# Set this path if you have chosen ATLAS (any default Linux BLAS):
export PATH_BLAS_ATLAS ?= /usr/lib/x86_64-linux-gnu
# Set this path to Intel root directory if you have chosen Intel MKL:
export PATH_INTEL ?= /opt/intel
#  Only reset these if Intel MKL libraries are spread in system directories:
 export PATH_BLAS_MKL ?= $(PATH_INTEL)/mkl/lib/intel64
 export PATH_BLAS_MKL_DEP ?= $(PATH_INTEL)/compilers_and_libraries/linux/lib/intel64_lin
 export PATH_BLAS_MKL_INC ?= $(PATH_INTEL)/mkl/include/intel64/lp64
# Set this path if you have chosen OpenBLAS:
export PATH_BLAS_OPENBLAS ?= /usr/lib/x86_64-linux-gnu
# Set this path if you have chosen ACML:
export PATH_BLAS_ACML ?= /opt/acml/5.3.1/gfortran64_fma4_mp/lib
# Set this path if you have chosen Cray LibSci:
export PATH_BLAS_LIBSCI ?= /opt/cray/pe/libsci/19.06.1/CRAY/8.5/x86_64/lib
# Set this path if you have chosen ESSL (also set PATH_IBM_XL_CPP, PATH_IBM_XL_FOR, PATH_IBM_XL_SMP below):
export PATH_BLAS_ESSL ?= /sw/summit/essl/6.1.0-2/essl/6.1/lib64

#IBM XL (only set these if you use IBM XL compiler and/or ESSL library):
export PATH_IBM_XL_CPP ?= /sw/summit/xl/16.1.1-3/xlC/16.1.1/lib
export PATH_IBM_XL_FOR ?= /sw/summit/xl/16.1.1-3/xlf/16.1.1/lib
export PATH_IBM_XL_SMP ?= /sw/summit/xl/16.1.1-3/xlsmp/5.1.1/lib

#LAPACK (only set these if you have chosen WITH_LAPACK=YES above):
export PATH_LAPACK_LIB ?= /usr/lib/x86_64-linux-gnu
export LAPACK_LIBS ?= -llapack

#CUDA (set these only if you build with CUDA):
export PATH_CUDA ?= /usr/local/cuda
# Only reset these if CUDA files are spread in system directories:
 export PATH_CUDA_INC ?= $(PATH_CUDA)/include
 export PATH_CUDA_LIB ?= $(PATH_CUDA)/lib64
 export PATH_CUDA_BIN ?= $(PATH_CUDA)/bin
# Reset your CUDA Host compiler if needed:
 export CUDA_HOST_COMPILER ?= /usr/bin/g++
# cuTT path (only if you use cuTT library):
export PATH_CUTT ?= /home/dima/src/cutt
# cuTensor path (only if you use cuTensor library):
export PATH_CUTENSOR ?= /home/dima/src/cutensor

#YOU ARE DONE! MAKE IT!


$(NAME):
ifeq ($(EXA_NO_BUILD),NO)
ifeq ($(EXA_TALSH_ONLY),NO)
	$(MAKE) -C ./UTILITY
	$(MAKE) -C ./GFC
	$(MAKE) -C ./DDSS
endif
	$(MAKE) -C ./TALSH
ifeq ($(EXA_TALSH_ONLY),NO)
	$(MAKE) -C ./DSVP
	$(MAKE) -C ./INTRAVIRT
	$(MAKE) -C ./INTERVIRT
#ifeq ($(EXA_OS),LINUX)
#	$(MAKE) -C ./TN
#endif
	$(MAKE) -C ./QFORCE
endif
#Gather headers, modules and libraries:
	rm -f ./include/*
	rm -f ./lib/*
	rm -f ./bin/*
ifeq ($(EXA_OS),LINUX)
ifeq ($(TOOLKIT),CRAY)
	cp -u ./[A-Z]*/OBJ/*.mod ./include/
else
	cp -u ./[A-Z]*/*.mod ./include/
endif
	cp -u ./TALSH/*.h ./include/
	cp -u ./TALSH/*.hpp ./include/
else
ifeq ($(TOOLKIT),CRAY)
	cp ./[A-Z]*/OBJ/*.mod ./include/
else
	cp ./[A-Z]*/*.mod ./include/
endif
	cp ./TALSH/*.h ./include/
	cp ./TALSH/*.hpp ./include/
endif
	cp ./[A-Z]*/*.a ./lib/
	cp ./[A-Z]*/*.x ./bin/
ifeq ($(EXA_TALSH_ONLY),NO)
	cp ./QFORCE/Qforce.x ./
#Create static and shared libraries:
	ar x ./lib/libintervirt.a
	ar x ./lib/libintravirt.a
	ar x ./lib/libdsvp.a
	ar x ./lib/libddss.a
	ar x ./lib/libtalsh.a
	ar x ./lib/libgfc.a
	ar x ./lib/libutility.a
	mv ./*.o ./lib/
	ar cr libexatensor.a ./lib/*.o
ifeq ($(EXA_OS),LINUX)
ifeq ($(WRAP),WRAP)
	CC -shared -o libexatensor.so ./lib/*.o
else
ifeq ($(TOOLKIT),IBM)
	$(PATH_$(MPILIB)_BIN)/mpicxx -qmkshrobj -o libexatensor.so ./lib/*.o
else
	$(PATH_$(MPILIB)_BIN)/mpicxx -shared -o libexatensor.so ./lib/*.o
endif
endif
	cp -u ./libexatensor.so ./lib/
	cp -u ./libexatensor.a ./lib/
	cp -u ./TALSH/libtalsh.so ./
	cp -u ./TALSH/libtalsh.so ./lib/
else
	cp ./libexatensor.a ./lib/
endif
	rm -f ./lib/*.o
else
ifeq ($(EXA_OS),LINUX)
	cp -u ./TALSH/libtalsh.so ./
	cp -u ./TALSH/libtalsh.so ./lib/
endif
endif
endif
	echo "Finished successfully!"

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.mod ./*/*.x ./*/*.a ./*/*.so ./*/*.mod ./*/OBJ/* ./bin/* ./lib/* ./include/*
