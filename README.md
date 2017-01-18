ExaTENSOR: Domain-specific virtual processing framework for
           numeric tensor algebra workloads on heterogeneous
           HPC systems.

AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov

Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

LICENSE: GNU Lesser General Public License v.3

PROGRAMMING MODEL: OOP FORTRAN 2003+, C/C++, OpenMP 3+, MPI 3+, CUDA.

DESCRIPTION: The project is still UNDER ACTIVE DEVELOPMENT although
             multiple components are ready for production. The ExaTENSOR
             framework further elaborates on the idea of domain-specific
             processing, that is, it introduces an intermediate virtual
             processing layer that can interpret domain-specific codes
             which express numeric TENSOR ALGEBRA workloads in our case.
             It is important to emphasize that ExaTENSOR is not tied to
             a specific science domain. Instead it can be used in any
             science domain which relies on tensor computations. So far,
             the ExaTENSOR framework aims at adaptive block-sparse tensor
             algebra workloads (includes dense and block-sparse tensor algebra).

             The distinguishing key design elements:
             a) Fully formalized domain-specific virtualization of heterogeneous
                compute nodes;
             b) Support of new accelerators via the device-unified API layer;
             c) Recursive virtualization of the HPC system;
             d) Support of multiple kinds of domain-specific virtual processors
                cooperating with each other;
             e) Recursive data decomposition and task scheduling;
             f) Adaptive data representation.
