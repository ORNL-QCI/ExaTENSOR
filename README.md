ExaTENSOR: Domain-specific virtual processing framework for
           numeric tensor algebra workloads on heterogeneous
           HPC systems.

AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov

Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

LICENSE: GNU Lesser General Public License v.3

PROGRAMMING MODEL: OOP FORTRAN 2003+, C/C++, OpenMP 3+, MPI 3+, CUDA.

DESCRIPTION: The project is still UNDER ACTIVE DEVELOPMENT although
             multiple components are already in production. The ExaTENSOR
             framework further elaborates on the idea of domain-specific
             virtual processing, that is, it introduces an intermediate
             virtual processing layer which can interpret a domain-specific
             code that expresses a numeric TENSOR ALGEBRA workload in our case.
             ExaTENSOR is not tied to any specific scientific domain and
             can be used for arbitrary numeric tensor computations. So far
             the ExaTENSOR framework has been aiming at adaptive hierarchical
             tensor algebra (includes dense and block-sparse tensor algebra).

             The distinguishing key design elements:
             a) Fully formalized domain-specific virtualization of heterogeneous
                compute nodes and HPC systems;
             b) Support of new accelerators via the device-unified API layer;
             c) Recursive virtualization of the HPC system;
             d) Support of multiple kinds of domain-specific virtual processors
                cooperating with each other;
             e) Recursive (hierarchical) data decomposition and task scheduling;
             f) Adaptive hierarchical tensor representation.
