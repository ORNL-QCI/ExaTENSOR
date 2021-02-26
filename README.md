ExaTENSOR is a basic numerical tensor algebra library for
distributed HPC systems equipped with multicore CPU and NVIDIA or AMD GPU.
The hierarchical task-based parallel runtime of ExaTENSOR
is based on the virtual tensor algebra processor architecture,
i.e. a software processor specialized to numerical tensor algebra
workloads on heterogeneous HPC systems (multicore/KNL, NVIDIA or AMD GPU).

(Details can be found here: https://doi.org/10.1002/qua.25926)

AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov

Copyright (C) 2014-2021 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2021 Oak Ridge National Laboratory (UT-Battelle)

LICENSE: GNU Lesser General Public License v.3

___
** This is NOT a production version yet as the work on fixing performance
   issues and addressing numerous problems in existing MPI-3 implementations
   is still in progress. Stay tuned. **
___

PROGRAMMING MODEL: OOP FORTRAN 2003+, C/C++11, OpenMP 3+, MPI 3+, CUDA, ROCM.

SUPPORTED COMPILERS: GNU 8.x.y, GNU 11+, Intel 18+, IBM XL 16+, Cray CCE.

LIBRARY DEPENDENCIES:
 * MPI-3 (OpenMPI/4.1+ or MPICH/3.2.1+);
 * CPU BLAS (optional: ATLAS, OpenBLAS, MKL, ACML, ESSL, LibSci);
 * cuBLAS (optional);
 * cuTT (optional);
 * cuTensor (optional);
 * ROCM/HIP (optional).

DESCRIPTION: The ExaTENSOR framework further elaborates on the idea of
             domain-specific virtual processing, that is, it introduces an
             intermediate virtual processing layer which can interpret
             domain-specific code expressing numerical TENSOR ALGEBRA workloads.
             ExaTENSOR is not tied to any specific scientific domain and
             can be used for arbitrary numerical tensor computations. So far
             the ExaTENSOR framework has been aiming at adaptive hierarchical
             tensor algebra (includes dense and block-sparse tensor algebra
             as special cases). The distinguishing key elements of ExaTENSOR:
 * Fully formalized design of domain-specific virtualization for
   heterogeneous compute nodes and full HPC systems;
 * Hierarchical virtualization of the full HPC system;
 * Hierarchical data/work decomposition and task scheduling;
 * Support of new hardware accelerators via the device-unified API layer;
 * Support of multiple kinds of domain-specific virtual processors
   with their own roles cooperating with each other;
 * Adaptive hierarchical tensor representation.

BUILD: Choose the right options for your platform in the Makefile header and make it.
       Supported compilers: GNU 8.x, GNU 11+, Intel 18+, IBM XL 16+, Cray CCE. Other
       compilers are still buggy in their support of modern Fortran 2003+, unfortunately.

       Configuration of the Makefile header via environment variables (you may
       export these environment variables in advance and Makefile will pick them up):
       a) WRAP: Change to WRAP if you use Cray compiler wrappers (ftn, cc, CC).
       b) TOOLKIT: Choose you compiler (prefer GNU).
       c) BUILD_TYPE: Choose OPT for optimized binary, DEV for debugging.
       d) MPILIB: Choose MPICH if you use MPICH or its derivative (e.g. Cray-MPICH)
          or choose OPENMPI if you use OpenMPI or its derivative (e.g. Spectrum-MPI).
          For either, you will also need to set the path to your chosen MPI library below.
       e) BLASLIB: Choose an optimized BLAS implementation if you have one. If you set it,
          you will also need to provide the path to your BLAS library below.
       f) GPU_CUDA: Set to CUDA if you have CUDA and want to use it.
       g) GPU_SM_ARCH: Set the two digits of your GPU compute capability if you use CUDA.
       h) EXA_OS: Change to NO_LINUX if you have not yet had a pleasure to work in Linux.
       i) WITH_CUTT: Choose YES if you have cuTT library and want to use it. cuTT fork:
          https://github.com/DmitryLyakh/cutt.
          You will also need to provide the path to the cuTT library below.

       If you do not use Cray compiler wrappers and associated modules, you need to set these as well:
       j) PATH_MPICH or PATH_OPENMPI: Set the corresponding path to your chosen MPI library root directory.
       k) PATH_BLAS_XXX: For you chose optimized BLAS implementation XXX, provide the path
          or multiple paths to your BLAS library files following the examples in the Makefile.
          If you chose IBM ESSL, you are also expected to provide paths to the IBM XL
          C/C++/Fortran runtime libraries via PATH_IBM_XL_XXX, XXX={CPP,FOR,SMP}.
       l) PATH_CUDA: If you chose to use CUDA, provide the path to it.
       m) PATH_CUTT: If you chose to use cuTT, provide the path to it.

RUN: An example script run.sh shows how to run ExaTENSOR for a test case (Qforce.x)
     implemented in QFORCE/main.F90. In general, at least 4 MPI processes are needed
     in order to run ExaTENSOR. The test configuration utilizes >= 4 MPI processes,
     each process running at least 5 OpenMP threads. This MPI test can still run on a
     single node via oversubscription of CPU cores. Environment variables to set:

     * QF_PATH: Path to the ExaTENSOR root directory.
     * QF_NUM_PROCS: Number of MPI processes, must be greater or equal to 4.
     * QF_PROCS_PER_NODE: Number of MPI processes per node (prefer binding to socket).
     * QF_CORES_PER_PROCESS: Number of CPU cores per MPI process. In case you have less
       CPU cores than the number of MPI processes, specify the minimum of 1.
     * QF_MEM_PER_PROCESS: Host RAM memory limit (MB) per MPI process.
     * QF_NVMEM_PER_PROCESS: Non-volatile memory limit (MB) per MPI process (if your node has one).
     * QF_HOST_BUFFER_SIZE: Size of the pinned Host RAM pool (MB): Set it to QF_MEM_PER_PROCESS.
     * QF_GPUS_PER_PROCESS: Number of exclusively owned NVIDIA GPUs per MPI process.
     * QF_NUM_THREADS: Initial number of threads per MPI process (8 or more).

     At the bottom of run.sh, pick or specify your MPI execution command for Qforce.x,
     taking into account the number of MPI processes per node, oversubscription, binding, etc.
     The test should normally take less than 5 minutes to run.

LINKING WITH OTHER APPS: See link.txt after building ExaTENSOR.

KNOWN ISSUES:
     (a) GNU 9.1 gfortran introduced a regression resulting in an Internal Compiler Error
         due to inability to generate finalization wrappers for tens_signature_t and
         tens_shape_t derived types in INTRAVIRT/tensor_recursive.F90. It is fixed in GNU 11+.
