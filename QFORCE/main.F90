!PROJECT Q-FORCE: Massively Parallel Quantum Many-Body Methodology on Heterogeneous HPC systems.
!BASE: ExaTensor: Massively Parallel Tensor Algebra Virtual Processor for Heterogeneous HPC systems.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/07/10

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

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
! - Fortran 2003 at least (some minor 2008 as well).
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
        use service_mpi
        implicit none
        integer(INTL), parameter:: TEST_SPACE_DIM=100_INTL
        type(spher_symmetry_t):: basis_symmetry(1:TEST_SPACE_DIM)
        type(subspace_basis_t):: basis
        class(h_space_t), pointer:: hspace
        type(tens_rcrsv_t):: dtens,ltens,rtens
        integer:: ierr,i,my_rank,space_id
        integer(INTL):: l

!Application initializes MPI:
        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
!Application runs ExaTENSOR within MPI_COMM_WORLD:
        ierr=exatns_start(MPI_COMM_WORLD)
        if(ierr.ne.EXA_SUCCESS) write(*,*) 'Process ',my_rank,' terminated with error ',ierr
!Create a basis for a vector space:
        if(my_rank.eq.0) then
         call basis%subspace_basis_ctor(TEST_SPACE_DIM,ierr)
         if(ierr.ne.0) call quit(ierr,'subspace_basis_t.subspace_basis_ctor() failed!')
         do l=1_INTL,TEST_SPACE_DIM !set basis functions
          call basis_symmetry(l)%spher_symmetry_ctor(int((l-1)/5,INTD),0,ierr)
          if(ierr.ne.0) call quit(ierr,'spher_symmetry_t.spher_symmetry_ctor() failed!')
          call basis%set_basis_func(l,BASIS_ABSTRACT,ierr,symm=basis_symmetry(l))
          if(ierr.ne.0) call quit(ierr,'subspace_basis_t.set_basis_func() failed!')
         enddo
         call basis%finalize(ierr)
         if(ierr.ne.0) call quit(ierr,'subspace_basis_t.finalize() failed!')
!Register a vector space:
         ierr=exatns_space_register('AO_space',basis,space_id)
         if(ierr.ne.0) call quit(ierr,'exatns_space_register() failed!')
!Create tensors:
         ierr=exatns_tensor_create(dtens,R8,'dtens',(/(space_id,i=1,4)/))
         ierr=exatns_tensor_create(ltens,R8,'ltens',(/(space_id,i=1,4)/))
         ierr=exatns_tensor_create(rtens,R8,'rtens',(/(space_id,i=1,4)/))
!Contract tensors:
         ierr=exatns_tensor_contract(dtens,ltens,rtens,'D(a,b,c,d)+=L(d,i,b,j)*R(j,c,i,a)')
!Destroy tensors:
         ierr=exatns_tensor_destroy(rtens)
         ierr=exatns_tensor_destroy(ltens)
         ierr=exatns_tensor_destroy(dtens)
        endif
!Master process stop ExaTENSOR:
        if(my_rank.eq.0.and.ierr.eq.0) ierr=exatns_stop()
!Application finalizes MPI:
        call MPI_Finalize(ierr)
       end program main
