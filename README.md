ExaTENSOR: Domain-specific virtual processor specialized to
           numerical tensor algebra workloads on heterogeneous
           HPC systems (multicore + NVIDIA GPU as of now).

AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov

Copyright (C) 2014-2018 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2018 Oak Ridge National Laboratory (UT-Battelle)

LICENSE: GNU Lesser General Public License v.3

PROGRAMMING MODEL: OOP FORTRAN 2003+, C/C++, OpenMP 3+, MPI 3+, CUDA.

DESCRIPTION: The project is still UNDER ACTIVE DEVELOPMENT although
             multiple components are already in production. The ExaTENSOR
             framework further elaborates on the idea of domain-specific
             virtual processing, that is, it introduces an intermediate
             virtual processing layer which can interpret domain-specific
             code that expresses numerical TENSOR ALGEBRA workloads.
             ExaTENSOR is not tied to any specific scientific domain and
             can be used for arbitrary numerical tensor computations. So far
             the ExaTENSOR framework has been aiming at adaptive hierarchical
             tensor algebra (includes dense and block-sparse tensor algebra).

             The distinguishing key design elements:
             a) Fully formalized domain-specific virtualization of heterogeneous
                compute nodes and HPC systems;
             b) Support of new accelerators via the device-unified API layer;
             c) Recursive (hierarchical) virtualization of the HPC system;
             d) Support of multiple kinds of domain-specific virtual processors
                cooperating with each other;
             e) Recursive (hierarchical) data decomposition and task scheduling;
             f) Adaptive hierarchical tensor representation.

BUILD: Choose the options right for your system in the Makefile header and make.
       Supported compilers: gcc 8+, Intel 18+, IBM XL 16.1.1 beta5+. Other compilers
       are still buggy in their support of modern Fortran 2003+, unfortunately.

RUN: An example script run.sh shows how to run ExaTENSOR for a test case (Qforce.x).
     At least 3 MPI processes are required to run ExaTENSOR. The test configuration
     utilizes 8 MPI processes, each process running at least 5 OpenMP threads.
