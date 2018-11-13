!PROJECT Q-FORCE: Massively Parallel Quantum Many-Body Methodology on Heterogeneous HPC systems.
!BASE: ExaTensor: Massively Parallel Tensor Algebra Virtual Processor for Heterogeneous HPC systems.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2018/11/11

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
! - C++11 at least.
! - MPI 3.0 at least.
! - OpenMP 3.0 at least (OpenMP 4.0+ if using Intel MIC).
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
        use service_mpi
        use stsubs, only: wait_delay
        use, intrinsic:: ISO_C_BINDING
        implicit none
        private

!Client (application) provides user-defined function objects:

 !Tensor initialization functor:
        type, extends(tens_method_uni_t), public:: tens_init_test_t
         complex(8), private:: init_val=(0d0,0d0)
         contains
          procedure, private:: TensInitTestCtorReal
          procedure, private:: TensInitTestCtorComplex
          generic, public:: tens_init_test_ctor=>TensInitTestCtorReal,TensInitTestCtorComplex
          procedure, public:: apply=>TensInitTestApply
        end type tens_init_test_t

 !Tensor transformation functor:
        type, extends(tens_method_uni_t), public:: tens_trans_test_t
         real(8), private:: div_val=1d0
         contains
          procedure, public:: tens_trans_test_ctor=>TensTransTestCtor
          procedure, public:: apply=>TensTransTestApply
        end type tens_trans_test_t

        public test_exatensor
        public benchmark_exatensor

       contains

!Client (application) provides implementation for user-defined function objects:

 !tens_init_test_t methods:
        subroutine TensInitTestCtorReal(this,val)
         implicit none
         class(tens_init_test_t), intent(out):: this
         real(8), intent(in):: val

         this%init_val=cmplx(val,0d0)
        end subroutine TensInitTestCtorReal

        subroutine TensInitTestCtorComplex(this,val)
         implicit none
         class(tens_init_test_t), intent(out):: this
         complex(8), intent(in):: val

         this%init_val=val
        end subroutine TensInitTestCtorComplex

        function TensInitTestApply(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(tens_init_test_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         integer(INTL):: vol
         type(tens_dense_t):: tens
         real(4), pointer:: arr_r4(:)
         real(8), pointer:: arr_r8(:)
         complex(4), pointer:: arr_c4(:)
         complex(8), pointer:: arr_c8(:)

         tens=tensor%get_dense_adapter(ierr)
         if(ierr.eq.EXA_SUCCESS) then
          vol=tens_dense_volume(tens)
          if(vol.gt.0) then
           select case(tens%data_kind)
           case(EXA_DATA_KIND_R4)
            call c_f_pointer(tens%body_ptr,arr_r4,(/vol/))
            arr_r4(:)=real(this%init_val,4)
           case(EXA_DATA_KIND_R8)
            call c_f_pointer(tens%body_ptr,arr_r8,(/vol/))
            arr_r8(:)=real(this%init_val,8)
           case(EXA_DATA_KIND_C4)
            call c_f_pointer(tens%body_ptr,arr_c4,(/vol/))
            arr_c4(:)=cmplx(real(this%init_val,4),real(imag(this%init_val),4),4)
           case(EXA_DATA_KIND_C8)
            call c_f_pointer(tens%body_ptr,arr_c8,(/vol/))
            arr_c8(:)=this%init_val
           case default
            ierr=2
           end select
          else
           ierr=1
          endif
         endif
        end function TensInitTestApply

 !tens_trans_test_t methods:
        subroutine TensTransTestCtor(this,val)
         implicit none
         class(tens_trans_test_t), intent(out):: this
         real(8), intent(in), optional:: val

         if(present(val)) this%div_val=val
        end subroutine TensTransTestCtor

        function TensTransTestApply(this,tensor,scalar) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(tens_trans_test_t), intent(in):: this
         class(tens_rcrsv_t), intent(inout):: tensor
         complex(8), intent(inout), optional:: scalar
         type(tens_dense_t):: tens

         tens=tensor%get_dense_adapter(ierr)
         !`Finish
        end function TensTransTestApply

        subroutine test_exatensor()
         implicit none
         integer(INTL), parameter:: TEST_SPACE_DIM=166_INTL !number of orbitals, must not exceed 166 here!
         integer(INTD), parameter:: color(1:166)=(/1,2,3,4,5,6,7,8,9,10,10,10,11,11,11,12,12,12,&
                                   &13,13,13,14,14,14,14,14,14,15,16,17,18,19,20,21,22,23,24,24,24,25,25,25,&
                                   &26,26,26,27,27,27,28,28,28,28,28,28,29,30,31,32,33,34,35,36,37,38,39,40,&
                                   &41,41,41,42,42,42,43,43,43,44,44,44,45,45,45,46,46,46,47,47,47,48,48,48,&
                                   &49,49,49,49,49,49,50,51,52,53,54,55,56,57,58,59,60,61,62,62,62,63,63,63,&
                                   &64,64,64,65,65,65,66,66,66,67,67,67,68,68,68,69,69,69,70,70,70,70,70,70,&
                                   &71,72,73,74,75,75,75,76,77,78,79,80,80,80,81,82,83,84,85,85,85,86,87,88,&
                                   &89,90,90,90/) !166 total
         type(color_symmetry_t):: basis_symmetry(1:TEST_SPACE_DIM)
         type(subspace_basis_t):: basis
         class(h_space_t), pointer:: ao_space
         type(tens_rcrsv_t):: etens,dtens,ltens,rtens
         type(tens_init_test_t):: init1369
         type(tens_trans_test_t):: div761
         type(tens_printer_t):: tens_printer
         integer(INTD):: ierr,i,my_rank,comm_size,ao_space_id,my_role,hsp(1:MAX_TENSOR_RANK)
         integer(INTL):: l,ao_space_root,ssp(1:MAX_TENSOR_RANK)
         complex(8):: etens_value
         real(8):: tms,tmf

         call MPI_Comm_size(MPI_COMM_WORLD,comm_size,ierr)
         call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
!Application creates a basis for a vector space:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Creating a basis for a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
         endif
         call basis%subspace_basis_ctor(TEST_SPACE_DIM,ierr)
         if(ierr.ne.0) call quit(ierr,'subspace_basis_t.subspace_basis_ctor() failed!')
         do l=1_INTL,TEST_SPACE_DIM !set basis functions
          call basis_symmetry(l)%color_symmetry_ctor(color(l),ierr)
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
!Print the registered space by levels (debug):
         if(my_rank.eq.comm_size-1) then
          !i=-1; do; i=i+1; call ao_space%print_level(i,num_subspaces=l); if(l.le.0) exit; enddo !debug
         endif
!Application registers user-defined methods for tensor initialization and printing:
 !Tensor initialization method:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering a user-defined tensor initialization method ... ")',ADVANCE='NO'); flush(6)
         endif
         call init1369%tens_init_test_ctor(13.69d0)
         ierr=exatns_method_register('SetTo13.69',init1369)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
 !Tensor transformation method:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering a user-defined tensor transformation method ... ")',ADVANCE='NO'); flush(6)
         endif
         call div761%tens_trans_test_ctor(7.61d2)
         ierr=exatns_method_register('DivideBy761',div761)
         if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_method_register() failed!')
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
 !Tensor printing method (defaults to screen):
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering a user-defined tensor printing method ... ")',ADVANCE='NO'); flush(6)
         endif
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
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
  !dtens:
           write(6,'("Creating tensor dtens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(dtens,'dtens',(/(ao_space_id,i=1,2)/),(/(ao_space_root,i=1,2)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Creating tensor ltens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(ltens,'ltens',(/(ao_space_id,i=1,2)/),(/(ao_space_root,i=1,2)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !rtens:
           write(6,'("Creating tensor rtens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(rtens,'rtens',(/(ao_space_id,i=1,2)/),(/(ao_space_root,i=1,2)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
 !Initialize tensors:
  !dtens:
           write(6,'("Initializing tensor dtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(dtens,(0d0,0d0))
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Initializing tensor ltens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(ltens,'SetTo13.69')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !rtens:
           write(6,'("Initializing tensor rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(rtens,'SetTo13.69')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
 !Transform input tensors:
  !ltens:
           write(6,'("Transforming tensor ltens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_transform(ltens,'DivideBy761')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_transform() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !rtens:
           write(6,'("Transforming tensor rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_transform(rtens,'DivideBy761')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_transform() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
 !Contract tensors:
           write(6,'("Contracting dtens+=ltens*rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(dtens,ltens,rtens,'D(a,b)+=L(c,a)*R(c,b)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Contract tensors:
           write(6,'("Contracting dtens+=ltens*rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(dtens,ltens,rtens,'D(a,b)+=L(c,a)*R(c,b)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
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
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
 !Print scalar etens:
           write(6,'("Printing scalar etens ... ")'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_traverse(etens,'PrintTensor')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_traverse() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Retrieve scalar etens directly:
           write(6,'("Retrieving directly scalar etens ... ")'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_get_scalar(etens,etens_value)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_get_scalar() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: Value = (",D21.14,1x,D21.14,"):",F16.4," sec")') etens_value,tmf-tms; flush(6)
 !Destroy tensors:
  !rtens:
           write(6,'("Destroying tensor rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(rtens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Destroying tensor ltens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(ltens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !dtens:
           write(6,'("Destroying tensor dtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(dtens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
  !etens:
           write(6,'("Destroying scalar etens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(etens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
 !Create some tensors again:
  !etens (scalar):
           write(6,'("Creating scalar etens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(etens,'etens',EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Destroy tensors again:
  !etens:
           write(6,'("Destroying scalar etens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(etens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           !ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Dump cache (debug):
           !ierr=exatns_dump_cache()
           !if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_dump_cache() failed!')
           !write(6,'("Tensor cache dumped")')
 !Stop ExaTENSOR runtime:
           ierr=exatns_stop()
          endif
         else
          write(6,*) 'Process ',my_rank,' terminated with error ',ierr
         endif
         return
        end subroutine test_exatensor

        subroutine benchmark_exatensor()
         implicit none
         integer(INTL), parameter:: TEST_SPACE_DIM=64
         integer(INTD), parameter:: BRANCHING_FACTOR=4
         type(color_symmetry_t):: basis_symmetry(1:TEST_SPACE_DIM)
         type(subspace_basis_t):: basis
         class(h_space_t), pointer:: ao_space
         type(tens_rcrsv_t):: etens,dtens,ltens,rtens
         type(tens_printer_t):: tens_printer
         integer(INTD):: ierr,i,my_rank,comm_size,ao_space_id,my_role,hsp(1:MAX_TENSOR_RANK)
         integer(INTL):: l,ao_space_root,ssp(1:MAX_TENSOR_RANK)
         complex(8):: etens_value
         real(8):: tms,tmf

         call MPI_Comm_size(MPI_COMM_WORLD,comm_size,ierr)
         call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
!Application creates a basis for a vector space:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Creating a basis for a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
         endif
         call basis%subspace_basis_ctor(TEST_SPACE_DIM,ierr)
         if(ierr.ne.0) call quit(ierr,'subspace_basis_t.subspace_basis_ctor() failed!')
         do l=1_INTL,TEST_SPACE_DIM !set basis functions
          call basis%set_basis_func(l,BASIS_ABSTRACT,ierr)
          if(ierr.ne.0) call quit(ierr,'subspace_basis_t.set_basis_func() failed!')
         enddo
         call basis%finalize(ierr)
         if(ierr.ne.0) call quit(ierr,'subspace_basis_t.finalize() failed!')
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
!Application registers a vector space:
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering the hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
         endif
         ierr=exatns_space_register('AOspace',basis,ao_space_id,ao_space,branch_factor=BRANCHING_FACTOR)
         if(ierr.ne.0) call quit(ierr,'exatns_space_register() failed!')
         if(my_rank.eq.comm_size-1) then; write(6,'("Ok")'); flush(6); endif
!Print the registered space by levels (debug):
         if(my_rank.eq.comm_size-1) then
          !i=-1; do; i=i+1; call ao_space%print_level(i,num_subspaces=l); if(l.le.0) exit; enddo !debug
         endif
!Application registers user-defined methods for tensor initialization and printing:
 !Tensor printing method (defaults to screen):
         if(my_rank.eq.comm_size-1) then
          write(6,'("Registering a user-defined tensor printing method ... ")',ADVANCE='NO'); flush(6)
         endif
         ierr=exatns_method_register('PrintTens',tens_printer)
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
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !dtens:
           write(6,'("Creating tensor dtens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(dtens,'dtens',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Creating tensor ltens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(ltens,'ltens',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !rtens:
           write(6,'("Creating tensor rtens over a hierarchical vector space ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_create(rtens,'rtens',(/(ao_space_id,i=1,4)/),(/(ao_space_root,i=1,4)/),EXA_DATA_KIND_R8)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_create() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Initialize tensors:
  !dtens:
           write(6,'("Initializing tensor dtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(dtens,(0d0,0d0))
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Initializing tensor ltens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(ltens,(1d-3,0d0))
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !rtens:
           write(6,'("Initializing tensor rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_init(rtens,(1d-2,0d0))
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_init() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Contract tensors:
           write(6,'("Contracting dtens+=ltens*rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(dtens,ltens,rtens,'D(a,b,c,d)+=L(d,i,c,j)*R(j,b,i,a)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Contract tensors again:
           write(6,'("Contracting dtens+=ltens*rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(dtens,ltens,rtens,'D(a,b,c,d)+=L(d,i,c,j)*R(j,b,i,a)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Contract tensors:
           write(6,'("Contracting etens+=dtens*dtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_contract(etens,dtens,dtens,'E()+=D(a,b,c,d)*D(a,b,c,d)')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_contract() failed!')
           ierr=exatns_sync(); if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_sync() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Print scalar etens:
           write(6,'("Printing scalar etens ... ")'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_traverse(etens,'PrintTens')
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_traverse() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Retrieve scalar etens directly:
           write(6,'("Retrieving directly scalar etens ... ")'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_get_scalar(etens,etens_value)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_get_scalar() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: Value = (",D21.14,1x,D21.14,"):",F16.4," sec")') etens_value,tmf-tms; flush(6)
 !Destroy tensors:
  !rtens:
           write(6,'("Destroying tensor rtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(rtens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !ltens:
           write(6,'("Destroying tensor ltens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(ltens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !dtens:
           write(6,'("Destroying tensor dtens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(dtens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
  !etens:
           write(6,'("Destroying scalar etens ... ")',ADVANCE='NO'); flush(6)
           tms=MPI_Wtime()
           ierr=exatns_tensor_destroy(etens)
           if(ierr.ne.EXA_SUCCESS) call quit(ierr,'exatns_tensor_destroy() failed!')
           tmf=MPI_Wtime()
           write(6,'("Ok: ",F16.4," sec")') tmf-tms; flush(6)
 !Stop ExaTENSOR runtime:
           ierr=exatns_stop()
          endif
         else
          write(6,*) 'Process ',my_rank,' terminated with error ',ierr
         endif
         return
        end subroutine benchmark_exatensor

       end module qforce_test


       program main
        use qforce_test
        use service_mpi
        implicit none
        integer(INT_MPI):: mpi_th_provided,ierr

!Application initializes MPI:
        call MPI_Init_Thread(MPI_THREAD_MULTIPLE,mpi_th_provided,ierr)
        if(mpi_th_provided.eq.MPI_THREAD_MULTIPLE) then
         call test_exatensor()
         !call benchmark_exatensor()
        else
         write(6,*) 'Your MPI library does not support MPI_THREAD_MULTIPLE! Change it!'
        endif
!Application finalizes MPI:
        call MPI_Finalize(ierr)
       end program main
