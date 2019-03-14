NAME = ExaTensor

#ADJUST THE FOLLOWING ENVIRONMENT VARIABLES ACCORDINGLY (choices are given)
#until you see "YOU ARE DONE!". The comments will guide you through (read them).
#Alternatively, you can export all relevant environment variables such that this
#Makefile will pick their values, so you will not need to update anything here.

#Cray cross-compiling wrappers (only for Cray): [WRAP|NOWRAP]:
export WRAP ?= NOWRAP
#Compiler: [GNU|PGI|INTEL|CRAY|IBM]:
export TOOLKIT ?= GNU
#Optimization: [DEV|OPT|PRF]:
export BUILD_TYPE ?= OPT
#MPI Library: [MPICH|OPENMPI]:
export MPILIB ?= MPICH
#BLAS: [ATLAS|MKL|ACML|ESSL|NONE]:
export BLASLIB ?= ATLAS
#Nvidia GPU via CUDA: [CUDA|NOCUDA]:
export GPU_CUDA ?= NOCUDA
#Nvidia GPU architecture (two digits, >=35):
export GPU_SM_ARCH ?= 35
#Operating system: [LINUX|NO_LINUX]:
export EXA_OS ?= LINUX


#ADJUST EXTRAS (optional):

#Fast GPU tensor transpose (cuTT library): [YES|NO]:
export WITH_CUTT ?= NO
#In-place GPU tensor contractions (cuTensor library): [YES|NO]:
export WITH_CUTENSOR ?= NO

#Disable actual build (debug): [YES|NO]:
export EXA_NO_BUILD ?= NO


#SET YOUR LOCAL PATHS (for unwrapped non-Cray builds):

#MPI library (whichever you have, set one):
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
# Set this path if you have chosen ATLAS (default Linux BLAS):
export PATH_BLAS_ATLAS ?= /usr/lib
# Set this path if you have chosen Intel MKL:
export PATH_INTEL ?= /opt/intel
#  Only reset these if Intel MKL libraries are spread in system directories:
export PATH_BLAS_MKL ?= $(PATH_INTEL)/mkl/lib/intel64
export PATH_BLAS_MKL_DEP ?= $(PATH_INTEL)/compilers_and_libraries/linux/lib/intel64_lin
export PATH_BLAS_MKL_INC ?= $(PATH_INTEL)/mkl/include/intel64/lp64
# Set this path if you have chosen ACML:
export PATH_BLAS_ACML ?= /opt/acml/5.3.1/gfortran64_fma4_mp/lib
# Set this path if you have chosen ESSL (also set PATH_IBM_XL_CPP, PATH_IBM_XL_FOR, PATH_IBM_XL_SMP below):
export PATH_BLAS_ESSL ?= /sw/summit/essl/6.1.0-2/essl/6.1/lib64

#IBM XL (only set these if you use IBM XL and/or ESSL):
export PATH_IBM_XL_CPP ?= /sw/summit/xl/16.1.1-1/xlC/16.1.1/lib
export PATH_IBM_XL_FOR ?= /sw/summit/xl/16.1.1-1/xlf/16.1.1/lib
export PATH_IBM_XL_SMP ?= /sw/summit/xl/16.1.1-1/xlsmp/5.1.1/lib

#CUDA (only if you build with CUDA):
export PATH_CUDA ?= /usr/local/cuda
# Only reset these if CUDA files are spread in system directories:
 export PATH_CUDA_INC ?= $(PATH_CUDA)/include
 export PATH_CUDA_LIB ?= $(PATH_CUDA)/lib64
 export PATH_CUDA_BIN ?= $(PATH_CUDA)/bin
 export CUDA_HOST_COMPILER ?= /usr/bin/g++
# cuTT path (only if you use cuTT library):
export PATH_CUTT ?= /home/dima/src/cutt
# cuTensor path (only if you use cuTensor library):
export PATH_CUTENSOR ?= /home/dima/src/cutensor

#YOU ARE DONE! MAKE IT!


$(NAME):
ifeq ($(EXA_NO_BUILD),NO)
	$(MAKE) -C ./UTILITY
	$(MAKE) -C ./GFC
	$(MAKE) -C ./DDSS
	$(MAKE) -C ./TALSH
	$(MAKE) -C ./DSVP
	$(MAKE) -C ./INTRAVIRT
	$(MAKE) -C ./INTERVIRT
ifeq ($(EXA_OS),LINUX)
	$(MAKE) -C ./TN
endif
	$(MAKE) -C ./QFORCE
#Gather headers, modules and libraries:
	rm -f ./include/*
	rm -f ./lib/*
	rm -f ./bin/*
ifeq ($(EXA_OS),LINUX)
	cp -u ./[A-Z]*/*.mod ./include/
	cp -u ./TALSH/*.h ./include/
	cp -u ./TALSH/*.hpp ./include/
else
	cp ./[A-Z]*/*.mod ./include/
	cp ./TALSH/*.h ./include/
	cp ./TALSH/*.hpp ./include/
endif
	cp ./[A-Z]*/*.a ./lib/
	cp ./[A-Z]*/*.x ./bin/
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
ifeq ($(TOOLKIT),GNU)
	g++ -shared -o libexatensor.so ./lib/*.o
else
	ld -shared -o libexatensor.so ./lib/*.o -lstdc++
endif
	cp -u ./TALSH/libtalsh.so ./
	cp -u ./TALSH/libtalsh.so ./lib/
	cp -u ./libexatensor.so ./lib/
	cp -u ./libexatensor.a ./lib/
else
	cp ./libexatensor.a ./lib/
endif
	rm -f ./lib/*.o
endif
	echo "Finished successfully!"

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.mod ./*/*.x ./*/*.a ./*/*.so ./*/*.mod ./*/OBJ/* ./bin/* ./lib/* ./include/*
