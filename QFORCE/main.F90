!PROJECT Q-FORCE: Massively Parallel Quantum Many-Body Methodology on Heterogeneous HPC systems.
!BASE: ExaTensor: Massively Parallel Tensor Algebra Virtual Processor for Heterogeneous HPC systems.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2018/09/21

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
! - USE_MPI_MOD - use MPI module instead of mpif.h;
! - USE_MKL - use Intel MKL library for BLAS/LAPACK;
!OUTPUT DEVICE:
! - jo (@service.mod) - default output device handle;
!NUMERATION OF DEVICES ON A NODE:
! - device_id = 0: Host (SMP CPU node, NUMA CPU node, Self-Hosted Intel Xeon Phi, etc.);
! - device_id = [1:MAX_GPUS_PER_NODE]: NVidia GPUs (GPU#=device_id-1);
! - device_id = [MAX_GPUS_PER_NODE+1:MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE]: Intel MICs (MIC#=device_id-1-MAX_GPUS_PER_NODE);
! - device_id = [MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+1:MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE]: AMD GPUs;

       module qforce_test
        use exatensor
        use stsubs, only: crash
        implicit none
        private

 !Tensor initialization functor:
        type, extends(tens_method_uni_t), public:: tens_init_test_t
         real(8), private:: init_val=0d0
         contains
          procedure, public:: tens_init_test_ctor=>TensInitTestCtor
          procedure, public:: apply=>TensInitTestApply
        end type tens_init_test_t

 !Tensor printing functor:
        type, extends(tens_method_uni_t), public:: tens_print_test_t
         integer(INTD), private:: dev_id=6
         contains
          procedure, public:: tens_print_test_ctor=>TensPrintTestCtor
          procedure, public:: apply=>TensPrintTestApply
        end type tens_print_test_t

       contains

        subroutine TensInitTestCtor(this,val)
         implicit none
         class(tens_init_test_t), intent(out):: this
         real(8), intent(in), optional:: val

         if(present(val)) this%init_val=val
        end subroutine TensInitTestCtor

        function TensInitTestApply(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(tens_init_test_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         class(tens_body_t), pointer:: body

         body=>tensor%get_body(ierr)
         if(ierr.eq.EXA_SUCCESS.and.associated(body)) then
          
         else
          ierr=-1
         endif
        end function TensInitTestApply

        subroutine TensPrintTestCtor(this,dev_out)
         implicit none
         class(tens_print_test_t), intent(out):: this
         integer(INTD), intent(in), optional:: dev_out

         if(present(dev_out)) this%dev_id=dev_out
        end subroutine TensPrintTestCtor

        function TensPrintTestApply(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(tens_print_test_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         class(tens_body_t), pointer:: body

         body=>tensor%get_body(ierr)
         if(ierr.eq.EXA_SUCCESS.and.associated(body)) then
          
         else
          ierr=-1
         endif
        end function TensPrintTestApply

       end module qforce_test


       program main
        use exatensor
        use service_mpi
        use stsubs, only: wait_delay
        use qforce_test
        implicit none
        integer(INTL), parameter:: TEST_SPACE_DIM=50_INTL
        type(spher_symmetry_t):: basis_symmetry(1:TEST_SPACE_DIM)
        type(subspace_basis_t):: basis
        class(h_space_t), pointer:: ao_space
        type(tens_rcrsv_t):: etens,dtens,ltens,rtens
        type(tens_init_test_t):: init1369
        type(tens_print_test_t):: tens_printer
        integer(INT_MPI):: mpi_th_provided
        integer(INTD):: ierr,i,my_rank,comm_size,ao_space_id,my_role,hsp(1:MAX_TENSOR_RANK)
        integer(INTL):: l,ao_space_root,ssp(1:MAX_TENSOR_RANK)
        real(8):: tms,tmf

!Application initializes MPI:
        call MPI_Init_Thread(MPI_THREAD_MULTIPLE,mpi_th_provided,ierr)
        if(mpi_th_provided.eq.MPI_THREAD_MULTIPLE) then
         call MPI_Comm_size(MPI_COMM_WORLD,comm_size,ierr)
         call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
!Application creates a basis for a vector space:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Creating a basis for a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
         endif
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
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
!Application registers a vector space:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering the hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
         endif
         ierr=exatns_space_register('AO_space',basis,ao_space_id,ao_space)
         if(ierr.ne.0) call quit(ierr,'exatns_space_register() failed!')
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
!Application registers user-defined methods for tensor initialization and printing:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering a user-defined tensor initialization method ... ")',ADVANCE='NO'); flush(6)
         endif
         call init1369%tens_init_test_ctor(13.69d0)
         ierr=exatns_method_register('SetTo13.69',init1369)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering a user-defined tensor printing method ... ")',ADVANCE='NO'); flush(6)
         endif
         call tens_printer%tens_print_test_ctor()
         ierr=exatns_method_register('PrintTensor',tens_printer)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
!Application runs ExaTENSOR within MPI_COMM_WORLD:
         ierr=exatns_start(MPI_COMM_WORLD)
         if(ierr.eq.EXA_SUCCESS) then
          ierr=exatns_process_role(my_role)
          if(my_role.eq.EXA_DRIVER) then
!Driver drives tensor workload:
 !Create tensors:
           ao_space_root=ao_space%get_root_id(ierr); if(ierr.ne.0) call quit(ierr,'h_space_t%get_root_id() failed!')
  !etens (scalar):
           write(6,'("Creating scalar etens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(etens,'etens',EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           ierr=exatns_dump_cache()
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           write(6,'("Tensor cache dumped")')
  !dtens:
           write(6,'("Creating tensor dtens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(dtens,'dtens',(/(ao_space_id,i=1,2)/),(/(ao_space_root,i=1,2)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Creating tensor ltens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(ltens,'ltens',(/(ao_space_id,i=1,2)/),(/(ao_space_root,i=1,2)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !rtens:
           write(6,'("Creating tensor rtens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(rtens,'rtens',(/(ao_space_id,i=1,2)/),(/(ao_space_root,i=1,2)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           ierr=exatns_dump_cache()
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           write(6,'("Tensor cache dumped")')
 !Initialize input tensors:
  !ltens:
           write(6,'("Initializing tensor ltens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(ltens,'SetTo13.69')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !rtens:
           write(6,'("Initializing tensor rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(rtens,'SetTo13.69')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           ierr=exatns_dump_cache()
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           write(6,'("Tensor cache dumped")')
#if 0
 !Contract tensors:
           write(6,'("Contracting etens+=ltens*rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(etens,ltens,rtens,'E()+=L(a,b)*R(a,b)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
#endif
 !Contract tensors:
           write(6,'("Contracting dtens+=ltens*rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(dtens,ltens,rtens,'D(a,b)+=L(c,a)*R(c,b)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Contract tensors:
           write(6,'("Contracting etens+=dtens*dtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(etens,dtens,dtens,'E()+=D(a,b)*D(a,b)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           ierr=exatns_dump_cache()
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           write(6,'("Tensor cache dumped")')
 !Destroy tensors:
  !rtens:
           write(6,'("Destroying tensor rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(rtens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Destroying tensor ltens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(ltens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !dtens:
           write(6,'("Destroying tensor dtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(dtens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           ierr=exatns_dump_cache()
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           write(6,'("Tensor cache dumped")')
  !etens:
           write(6,'("Destroying scalar etens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(etens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           ierr=exatns_dump_cache()
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           write(6,'("Tensor cache dumped")')
 !Stop ExaTENSOR runtime:
           ierr=exatns_stop()
          endif
         else
          write(6,*) 'Process ',my_rank,' terminated with error ',ierr
         endif
        else
         write(6,*) 'Your MPI library does not support MPI_THREAD_MULTIPLE! Change it!'
        endif
!Application finalizes MPI:
        call MPI_Finalize(ierr)
       end program main
