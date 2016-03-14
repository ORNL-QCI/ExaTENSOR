!PROJECT: Q-FORCE: Massively Parallel Quantum Many-Body Methodology on Heterogeneous HPC systems.
!BASE: ExaTensor: Massively Parallel Tensor Algebra Virtual Processor for Heterogeneous HPC systems.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/09/25

!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

!This file is part of ExaTensor.

!ExaTensor is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published
!by the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!ExaTensor is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU Lesser General Public License for more details.

!You should have received a copy of the GNU Lesser General Public License
!along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.

!COMPILATION:
! - Fortran 2003 at least (some 2008 as well).
! - MPI 3.0 at least.
! - OpenMP 3.0 at least (OpenMP 4.0 if using Intel MIC).
! - CUDA 5.0 at least.
!MPI launch notes:
! - MPI processes launched on the same node MUST have consecutive numbers!
! - Environment variable <QF_PROCS_PER_NODE> MUST be defined when using accelerators (number of processes per node)!
!PREPROCESSOR:
! - NO_AMD - do not use AMD GPU;
! - NO_PHI - do not use Intel Xeon Phi (MIC);
! - NO_GPU - do not use NVidia GPU (CUDA);
! - NO_BLAS - BLAS/LAPACK calls will be replaced by in-house routines (D.I.L.);
! - USE_OMP_MOD - use OpenMP module;
! - USE_MPI_MOD - use MPI module instead of mpif.h;
! - USE_MKL - use Intel MKL library for BLAS/LAPACK;
!OUTPUT DEVICE:
! - jo (@service.mod) - default output device handle;
!NUMERATION OF DEVICES ON A NODE:
! - device_id = 0: Host (SMP CPU node, NUMA CPU node, Self-Hosted Intel Xeon Phi, etc.);
! - device_id = [1:MAX_GPUS_PER_NODE]: NVidia GPUs (GPU#=device_id-1);
! - device_id = [MAX_GPUS_PER_NODE+1:MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE]: Intel MICs (MIC#=device_id-1-MAX_GPUS_PER_NODE);
! - device_id = [MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+1:MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE]: AMD GPUs;
       program main
        use exatensor
        implicit none
        integer:: ierr
        ierr=0
        call exa_tensor(ierr)
        if(ierr.ne.0) then
         write(*,*) 'ExaTensor terminated with an error: ',ierr
        else
         write(*,*) 'ExaTensor finished successfully!'
        endif
       end program main
