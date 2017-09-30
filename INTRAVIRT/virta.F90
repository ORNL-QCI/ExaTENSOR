!ExaTENSOR: Parallel Virtual Processing for Scale-Adaptive Hierarchical Tensor Algebra:
!Derives from abstract domain-specific virtual processing (DSVP).
!This module provides basic infrastructure for ExaTENSOR, tensor algebra virtual processor (TAVP).
!Different specializations (roles) of tensor algebra virtual processors derive from this module.
!Note that, in general, different specializations (roles) of the domain-specific virtual processor
!may have differing instruction sets (non-overlapping, overlapping, or identical). The instruction
!codes provided in this module are common for all specializations of the tensor algebra virtual processors.
!However, different specializations always have different microcodes, even for the same instruction codes.

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/09/30

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

       module virta !VIRtual Tensor Algebra
        use tensor_algebra     !basic constants
        use talsh              !on-node heterogeneous numeric tensor algebra
        use pack_prim          !object packing/unpacking primitives
        use distributed        !distributed one-sided communication layer
        use service_mpi        !basic MPI service
        use hardware           !hardware abstraction and aggregation
        use subspaces          !hierarchical vector space representation
        use tensor_recursive   !recursive (hierarchical) tensors
        use dsvp_base          !abstract domain-specific virtual processor (DSVP)
        use gfc_base           !GFC base
        use gfc_dictionary     !GFC dictionary
        implicit none
        public
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Runtime errors (ExaTENSOR aliases of basic DSVP errors):
        integer(INTD), parameter, public:: EXA_SUCCESS=DSVP_SUCCESS                         !success
        integer(INTD), parameter, public:: EXA_ERROR=DSVP_ERROR                             !generic error
        integer(INTD), parameter, public:: EXA_ERR_INVALID_ARGS=DSVP_ERR_INVALID_ARGS       !invalid arguments passed
        integer(INTD), parameter, public:: EXA_ERR_INVALID_REQ=DSVP_ERR_INVALID_REQ         !invalid request
        integer(INTD), parameter, public:: EXA_ERR_MEM_ALLOC_FAIL=DSVP_ERR_MEM_ALLOC_FAIL   !memory allocation failed
        integer(INTD), parameter, public:: EXA_ERR_MEM_FREE_FAIL=DSVP_ERR_MEM_FREE_FAIL     !memory deallocation failed
        integer(INTD), parameter, public:: EXA_ERR_BROKEN_OBJ=DSVP_ERR_BROKEN_OBJ           !broken object
        integer(INTD), parameter, public:: EXA_ERR_UNABLE_COMPLETE=DSVP_ERR_UNABLE_COMPLETE !unable to complete
        integer(INTD), parameter, public:: EXA_ERR_RSC_EXCEEDED=DSVP_ERR_RSC_EXCEEDED       !resource exceeded
 !Tensor-algebra virtual processor (TAVP) kinds (roles):
        integer(INTD), parameter, public:: EXA_NO_ROLE=DSVP_NO_KIND !undefined role
        integer(INTD), parameter, public:: EXA_DRIVER=0             !domain-specific driver process (not a TAVP)
        integer(INTD), parameter, public:: EXA_MANAGER=1            !manager (logical) process (TAVP)
        integer(INTD), parameter, public:: EXA_WORKER=2             !worker (numerical) process (TAVP)
        integer(INTD), parameter, public:: EXA_HELPER=3             !helper (auxiliary) process (TAVP)
 !TAVP hierarchy configuration:
        integer(INTD), public:: EXA_MAX_WORK_GROUP_SIZE=64 !maximal size of a work group (max number of workers per manager)
        integer(INTD), public:: EXA_MANAGER_BRANCH_FACT=2  !branching factor for the managing hierarchy
 !TAVP identification:
        integer(INTD), parameter, public:: TAVP_ANY_ID=-1              !any TAVP
 !TAVP instruction error codes:
        integer(INTD), parameter, public:: TAVP_ERR_GEN_FAILURE=-1     !unspecified generic failure
        integer(INTD), parameter, public:: TAVP_ERR_BTC_BAD=-2         !bad instruction bytecode
        integer(INTD), parameter, public:: TAVP_ERR_CHE_FAILURE=-3     !argument cache failure
        integer(INTD), parameter, public:: TAVP_ERR_ARG_UNDEFINED=-4   !instruction argument is undefined in the argument cache
        integer(INTD), parameter, public:: TAVP_ERR_ARG_DEFINED=-5     !instrcution argument is already defined in the argument cache
        integer(INTD), parameter, public:: TAVP_ERR_RSC_UNAVAILABLE=-6 !instruction is unable to obtain needed resources
        integer(INTD), parameter, public:: TAVP_ERR_RSC_FAILURE=-7     !instruction is unable to release resources cleanly
        integer(INTD), parameter, public:: TAVP_ERR_COM_FAILURE=-8     !instruction communication failed
        integer(INTD), parameter, public:: TAVP_ERR_EXC_FAILURE=-9     !instruction computation failed
 !TAVP ISA max size (the same for all TAVP specializations):
        integer(INTD), parameter, public:: TAVP_ISA_SIZE=256 !max number of TAVP instructions [0:TAVP_ISA_SIZE-1]
 !Tensor algebra virtual processor (TAVP):
  !TAVP instruction code (opcode), must be non-negative (consult TAProL spec), limited by TAVP_ISA_SIZE:
   !NOOP:
        integer(INTD), parameter, public:: TAVP_INSTR_NOOP=DS_INSTR_NOOP !no operation (empty instruction): Negative opcode
   !General control [0-15]:
        integer(INTD), parameter, public:: TAVP_INSTR_STOP=0      !stop TAVP (finishes current instructions and shutdowns TAVP)
        integer(INTD), parameter, public:: TAVP_INSTR_PAUSE=1     !pause TAVP execution (finishes active instructions and pauses TAVP)
        integer(INTD), parameter, public:: TAVP_INSTR_RESUME=2    !resume TAVP execution (resumes TAVP execution pipeline after a pause)
   !Auxiliary definitions [16-63]:
        integer(INTD), parameter, public:: TAVP_INSTR_SPACE_CREATE=16  !create a vector space
        integer(INTD), parameter, public:: TAVP_INSTR_SPACE_DESTROY=17 !destroy a vector space
   !Tensor operations [64-255]:
        integer(INTD), parameter, public:: TAVP_INSTR_CREATE=64   !tensor creation
        integer(INTD), parameter, public:: TAVP_INSTR_DESTROY=65  !tensor destruction
        integer(INTD), parameter, public:: TAVP_INSTR_LOAD=66     !tensor loading from persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_SAVE=67     !tensor saving in persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_COMM=68     !tensor communication (get/put/accumulate)
        integer(INTD), parameter, public:: TAVP_INSTR_INIT=69     !tensor initialization (assignement of a value)
        integer(INTD), parameter, public:: TAVP_INSTR_NORM1=70    !tensor 1-norm
        integer(INTD), parameter, public:: TAVP_INSTR_NORM2=71    !tensor 2-norm
        integer(INTD), parameter, public:: TAVP_INSTR_MIN=72      !tensor min element
        integer(INTD), parameter, public:: TAVP_INSTR_MAX=73      !tensor max element
        integer(INTD), parameter, public:: TAVP_INSTR_FOLD=74     !tensor dimension folding
        integer(INTD), parameter, public:: TAVP_INSTR_UNFOLD=75   !tensor dimension unfolding
        integer(INTD), parameter, public:: TAVP_INSTR_SLICE=76    !tensor slicing (taking a slice of a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_INSERT=77   !tensor insertion (inserting a tensor slice in a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_PERMUTE=78  !tensor dimension permutation (in-place)
        integer(INTD), parameter, public:: TAVP_INSTR_COPY=79     !tensor copy (copying a tensor into another tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_SCALE=80    !tensor scaling (multiplication by a number)
        integer(INTD), parameter, public:: TAVP_INSTR_ADD=81      !tensor addition (with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TRACE=82    !tensor trace (tracing over some/all tensor indices)
        integer(INTD), parameter, public:: TAVP_INSTR_PRODUCT=83  !tensor direct product (with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_CONTRACT=84 !tensor contraction (also includes tensor product, tensor addition, and tensor scaling)
!TYPES:
 !Tensor status:
        type, public:: tens_status_t
         logical, public:: created=.FALSE. !TRUE if the tensor has been created (allocated physical memory), FALSE otherwise
         logical, public:: defined=.FALSE. !TRUE if the tensor value has been defined, FALSE otherwise
         logical, public:: replica=.FALSE. !TRUE if the tensor is a replica (another defined instance of this same tensor exists), FALSE otherwise
         logical, public:: updated=.FALSE. !TRUE if the tensor value is currently being updated, FALSE otherwise
         integer(INTD), public:: is_used=0 !number of read-only references which are currently using the tensor
        end type tens_status_t
        type(tens_status_t), parameter:: tens_status_none=tens_status_t() !tensor status null
        public tens_status_none
 !Tensor instruction control fields:
#if 0
  !Tensor addition/copy/slice/insert control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_add_t
         type(permutation_t), private:: prmn                    !tensor dimension permutation
         complex(8), private:: alpha=(1d0,0d0)                  !alpha prefactor
         contains
          procedure, private:: CtrlTensAddCtor                  !ctor
          generic, public:: ctrl_tens_add_ctor=>CtrlTensAddCtor
          procedure, public:: pack=>CtrlTensAddPack             !packs the instruction control field into a plain byte packet
          procedure, public:: unpack=>CtrlTensAddUnpack         !unpacks the instruction control field from a plain byte packet
          final:: ctrl_tens_add_dtor                            !dtor
        end type ctrl_tens_add_t
#endif
  !Tensor contraction control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_contr_t
         type(contr_ptrn_ext_t), private:: contr_ptrn           !extended tensor contraction pattern
         complex(8), private:: alpha=(1d0,0d0)                  !alpha prefactor
         contains
          procedure, private:: CtrlTensContrCtor                !ctor
          generic, public:: ctrl_tens_contr_ctor=>CtrlTensContrCtor
          procedure, public:: pack=>CtrlTensContrPack           !packs the instruction control field into a plain byte packet
          procedure, public:: unpack=>CtrlTensContrUnpack       !unpacks the instruction control field from a plain byte packet
          final:: ctrl_tens_contr_dtor                          !dtor
        end type ctrl_tens_contr_t
 !Tensor argument cache entry:
        type, abstract, public:: tens_cache_entry_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL() !either owning or non-owning pointer to a tensor (ctor/dtor of extended types will decide)
         type(tens_status_t), private:: tens_status             !current status of the tensor
         contains
          procedure, public:: set_tensor=>TensCacheEntrySetTensor     !sets the pointer to a tensor
          procedure, public:: is_set=>TensCacheEntryIsSet             !returns TRUE if the tensor cache entry is set (constructed)
          procedure, public:: get_tensor=>TensCacheEntryGetTensor     !returns a non-owning pointer to the tensor
          procedure, public:: get_status=>TensCacheEntryGetStatus     !returns the tensor status object
          procedure, public:: mark_created=>TensCacheEntryMarkCreated !marks the tensor status as created (allocated memory)
          procedure, public:: mark_defined=>TensCacheEntryMarkDefined !marks the tensor status as defined to some value
          procedure, public:: mark_in_use=>TensCacheEntryMarkInUse    !marks the tensor status as in-use (read-only) and increases the reference count
          procedure, public:: mark_no_use=>TensCacheEntryMarkNoUse    !decreases the read-only usage reference count
          procedure, public:: mark_updated=>TensCacheEntryMarkUpdated !marks the tensor status as in-update (currently being updated)
          procedure, public:: mark_empty=>TensCacheEntryMarkEmpty     !marks the tensor status as empty (destroyed)
          procedure, public:: nullify_tensor=>TensCacheEntryNullifyTensor !either deallocates or simply dissociates tensor (depends on the extended dtor)
        end type tens_cache_entry_t
 !Tensor argument cache:
        type, public:: tens_cache_t
         type(dictionary_t), private:: map                        !cache dictionary: <tens_descr_t-->tens_cache_entry_t>
         contains
          procedure, public:: lookup=>TensCacheLookup             !looks up a given tensor in the cache
          procedure, public:: store=>TensCacheStore               !stores a given tensor in the cache
          procedure, public:: evict=>TensCacheEvict               !evicts a given tensor from the cache
          procedure, public:: erase=>TensCacheErase               !erases everything from the cache (regardless of pending MPI communications!)
          final:: tens_cache_dtor                                 !dtor
        end type tens_cache_t
!INTERFACES:
        abstract interface
 !tens_cache_entry_t:
  !Polymorphic allocator:
         function tens_cache_entry_alloc_i(tens_cache_entry) result(ierr)
          import:: tens_cache_entry_t,INTD
          integer(INTD):: ierr
          class(tens_cache_entry_t), allocatable, intent(out):: tens_cache_entry
         end function tens_cache_entry_alloc_i
        end interface
!DATA:
 !MPI process specialization (TAVP role, set by exatns_start):
        integer(INT_MPI), public:: process_role=EXA_NO_ROLE   !MPI process role (see above)
        integer(INT_MPI), public:: role_comm=MPI_COMM_NULL    !role-specific MPI communicator
        integer(INT_MPI), public:: role_size=0                !size of the role-specific MPI communicator
        integer(INT_MPI), public:: role_rank=-1               !process rank within the role-specific MPI communicator
        integer(INT_MPI), public:: driver_comm=MPI_COMM_NULL  !MPI intracommunicator of the driver process subspace
        integer(INT_MPI), public:: manager_comm=MPI_COMM_NULL !MPI intracommunicator of the manager process subspace
        integer(INT_MPI), public:: worker_comm=MPI_COMM_NULL  !MPI intracommunicator of the worker process subspace
        integer(INT_MPI), public:: helper_comm=MPI_COMM_NULL  !MPI intracommunicator of the helper process subspace
        integer(INT_MPI), public:: drv_mng_comm=MPI_COMM_NULL !MPI intercommunicator for the driver and managers
        integer(INT_MPI), public:: mng_wrk_comm=MPI_COMM_NULL !MPI intercommunicator for managers and workers
        integer(INT_MPI), public:: driver_mpi_process=0       !MPI rank of the Driver process (normally process 0 from MPI_COMM_WORLD)
        integer(INT_MPI), allocatable, public:: managers(:)   !ordered list of MPI processes serving as managers (ranks in MPI_COMM_WORLD)
!VISIBILITY:
 !non-member:
        public role_barrier
#if 0
 !ctrl_tens_add_t:
        private CtrlTensAddCtor
        private CtrlTensAddPack
        private CtrlTensAddUnpack
        public ctrl_tens_add_dtor
#endif
 !ctrl_tens_contr_t:
        private CtrlTensContrCtor
        private CtrlTensContrPack
        private CtrlTensContrUnpack
        public ctrl_tens_contr_dtor
 !tens_cache_entry_t:
        private TensCacheEntrySetTensor
        private TensCacheEntryIsSet
        private TensCacheEntryGetTensor
        private TensCacheEntryGetStatus
        private TensCacheEntryMarkCreated
        private TensCacheEntryMarkDefined
        private TensCacheEntryMarkInUse
        private TensCacheEntryMarkNoUse
        private TensCacheEntryMarkUpdated
        private TensCacheEntryMarkEmpty
        private TensCacheEntryNullifyTensor
 !tens_cache_t:
        private TensCacheLookup
        private TensCacheStore
        private TensCacheEvict
        private TensCacheErase
        public tens_cache_dtor

       contains
!IMPLEMENTATION:
![non-member]========================
        subroutine role_barrier(ierr)
!MPI barrier for the role communicator.
         implicit none
         integer(INTD), intent(out), optional:: ierr
         integer(INT_MPI):: errc

         call MPI_Barrier(role_comm,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine role_barrier
![ctrl_tens_contr_t]============================================
        subroutine CtrlTensContrCtor(this,contr_ptrn,ierr,alpha)
!CTOR.
         implicit none
         class(ctrl_tens_contr_t), intent(out):: this    !out: tensor contraction control field
         type(contr_ptrn_ext_t), intent(in):: contr_ptrn !in: extended tensor contraction pattern
         integer(INTD), intent(out), optional:: ierr     !out: error code
         complex(8), intent(in), optional:: alpha        !in: scalar prefactor
         integer(INTD):: errc

         if(contr_ptrn%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           this%contr_ptrn=contr_ptrn
           if(present(alpha)) this%alpha=alpha
          else
           errc=-1
          endif
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensContrCtor
!-----------------------------------------------------
        subroutine CtrlTensContrPack(this,packet,ierr)
!Packer.
         implicit none
         class(ctrl_tens_contr_t), intent(in):: this !in: tensor contraction control field
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%contr_ptrn%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           call this%contr_ptrn%pack(packet,errc)
           if(errc.eq.TEREC_SUCCESS) then
            call pack_builtin(packet,this%alpha,errc); if(errc.ne.0) errc=-1
           else
            errc=-2
           endif
          else
           errc=-3
          endif
         else
          errc=-4
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensContrPack
!-------------------------------------------------------
        subroutine CtrlTensContrUnpack(this,packet,ierr)
!Unpacker.
         implicit none
         class(ctrl_tens_contr_t), intent(inout):: this !out: tensor contraction control field
         class(obj_pack_t), intent(inout):: packet      !inout: packet
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         call this%contr_ptrn%unpack(packet,errc)
         if(errc.eq.TEREC_SUCCESS) then
          call unpack_builtin(packet,this%alpha,errc); if(errc.ne.0) errc=-1
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensContrUnpack
!--------------------------------------------
        subroutine ctrl_tens_contr_dtor(this)
         implicit none
         type(ctrl_tens_contr_t):: this

         return
        end subroutine ctrl_tens_contr_dtor
![tens_cache_entry_t]=======================================
        subroutine TensCacheEntrySetTensor(this,tensor,ierr)
!Sets a pointer to the tensor.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this   !inout: empty tensor cache entry
         class(tens_rcrsv_t), intent(in), pointer:: tensor !in: associated pointer to a tensor
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=0
         if((.not.associated(this%tensor)).and.associated(tensor)) then
          this%tensor=>tensor
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntrySetTensor
!----------------------------------------------------------
        function TensCacheEntryIsSet(this,ierr) result(ans)
         implicit none
         logical:: ans                                !out: answer
         class(tens_cache_entry_t), intent(in):: this !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=0; ans=associated(this%tensor)
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEntryIsSet
!-------------------------------------------------------------------
        function TensCacheEntryGetTensor(this,ierr) result(tensor_p)
         implicit none
         class(tens_rcrsv_t), pointer:: tensor_p      !out: pointer to the tensor
         class(tens_cache_entry_t), intent(in):: this !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=0; tensor_p=>this%tensor
         if(.not.associated(tensor_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEntryGetTensor
!--------------------------------------------------------------------
        function TensCacheEntryGetStatus(this,ierr) result(tens_stat)
         implicit none
         type(tens_status_t):: tens_stat              !out: pointer to the tensor status
         class(tens_cache_entry_t), intent(in):: this !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         tens_stat=this%tens_status
         if(.not.this%is_set(errc)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEntryGetStatus
!------------------------------------------------------
        subroutine TensCacheEntryMarkCreated(this,ierr)
!Marks the tensor as CREATED, but not yet DEFINED.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           if(.not.this%tens_status%created) then
            this%tens_status%created=.TRUE.
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryMarkCreated
!------------------------------------------------------
        subroutine TensCacheEntryMarkDefined(this,ierr)
!Marks the tensor as DEFINED to VALUE, as a result of an UPDATE.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           if(this%tens_status%created) then
            if(this%tens_status%updated) then
             if(.not.(this%tens_status%defined.or.this%tens_status%is_used.ne.0)) then
              this%tens_status%updated=.FALSE.
              this%tens_status%defined=.TRUE.
             else
              errc=-5
             endif
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryMarkDefined
!----------------------------------------------------
        subroutine TensCacheEntryMarkInUse(this,ierr)
!Associates a reference to the tensor.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           if(this%tens_status%created) then
            if(this%tens_status%defined) then
             this%tens_status%is_used=this%tens_status%is_used+1
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryMarkInUse
!----------------------------------------------------
        subroutine TensCacheEntryMarkNoUse(this,ierr)
!Dissociates a reference from the tensor.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           if(this%tens_status%created) then
            if(this%tens_status%defined) then
             if(this%tens_status%is_used.gt.0) then
              this%tens_status%is_used=this%tens_status%is_used-1
             else
              errc=-5
             endif
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryMarkNoUse
!------------------------------------------------------
        subroutine TensCacheEntryMarkUpdated(this,ierr)
!Marks the tensor status as BEING UPDATED --> NOT DEFINED.
!The tensor must not be IN USE at this time.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           if(this%tens_status%created) then
            if(this%tens_status%is_used.eq.0) then
             this%tens_status%defined=.FALSE.
             this%tens_status%updated=.TRUE.
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryMarkUpdated
!----------------------------------------------------
        subroutine TensCacheEntryMarkEmpty(this,ierr)
!Marks the tensor status as EMPTY (destroyed).
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           if(this%tens_status%created) then
            if((this%tens_status%is_used.eq.0).and.(.not.this%tens_status%updated)) then
             this%tens_status%defined=.FALSE.
             this%tens_status%created=.FALSE.
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryMarkEmpty
!-----------------------------------------------------------
        subroutine TensCacheEntryNullifyTensor(this,dealloc)
!Either deallocates or simply dissociates tensor (depends on the extended dtor).
         implicit none
         class(tens_cache_entry_t), intent(inout):: this
         logical, intent(in):: dealloc

         if(associated(this%tensor).and.dealloc) deallocate(this%tensor)
         this%tensor=>NULL()
         return
        end subroutine TensCacheEntryNullifyTensor
![tens_cache_t]========================================================
        function TensCacheLookup(this,tensor,ierr) result(tens_entry_p)
!Looks up a given tensor in the tensor cache. If found, returns a pointer
!to the corresponding tensor cache entry. If not found, returns NULL.
         implicit none
         class(tens_cache_entry_t), pointer:: tens_entry_p !out: pointer to the found tensor cache entry or NULL
         class(tens_cache_t), intent(in):: this            !in: tensor cache
         class(tens_rcrsv_t), intent(in):: tensor          !in: tensor to look up (via its descriptor as the key)
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         type(dictionary_iter_t):: dit
         type(tens_descr_t), target:: tens_descr

         tens_entry_p=>NULL()
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>NULL()
          tens_descr=tensor%get_descriptor(errc)
          if(errc.eq.TEREC_SUCCESS) then
           res=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_tens_descriptors,tens_descr,value_out=uptr)
           if(res.eq.GFC_FOUND) then
            select type(uptr); class is(tens_cache_entry_t); tens_entry_p=>uptr; end select
            if(.not.associated(tens_entry_p)) errc=-5 !trap
           else
            if(res.ne.GFC_NOT_FOUND) errc=-4
           endif
          else
           errc=-3
          endif
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheLookup
!-----------------------------------------------------------------------------------------------------
        function TensCacheStore(this,tensor,tens_cache_entry_alloc_f,ierr,tens_entry_p) result(stored)
!Given a tensor, checks whether it is present in the tensor cache. If yes, returns
!a pointer to the corresponding tensor cache entry. If no, allocates a new extended
!tensor cache entry, stores it in the tensor cache and returns a pointer to it.
!The allocation is done via a user-provided non-member allocator of an extension of
!tens_cache_entry_t supplied with a proper dtor. Note that the newly created extended
!tensor cache entry may need further construction. Here only the <tensor> component of
!an extended tens_cache_entry_t is set.
         implicit none
         logical:: stored                                                !out: TRUE on successful new store, FALSE otherwise
         class(tens_cache_t), intent(inout):: this                       !inout: tensor cache
         class(tens_rcrsv_t), pointer, intent(in):: tensor               !in: pointer to a tensor to be stored
         procedure(tens_cache_entry_alloc_i):: tens_cache_entry_alloc_f  !in: non-member allocator of an extended(tens_cache_entry_t) class with a proper dtor
         integer(INTD), intent(out), optional:: ierr                     !out: error code
         class(tens_cache_entry_t), pointer, intent(out), optional:: tens_entry_p !out: tensor cache entry (newly created or existing)
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         type(dictionary_iter_t):: dit
         type(tens_descr_t), target:: tens_descr
         class(tens_cache_entry_t), allocatable, target:: tce
         class(tens_cache_entry_t), pointer:: tcep

         stored=.FALSE.
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          if(associated(tensor)) then
           tens_descr=tensor%get_descriptor(errc)
           if(errc.eq.TEREC_SUCCESS) then
            errc=tens_cache_entry_alloc_f(tce) !allocates an empty instance of an extended(tens_cache_entry_t) with a proper dtor
            if(errc.eq.0) then
             uptr=>NULL()
             res=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_tens_descriptors,tens_descr,tce,GFC_BY_VAL,value_out=uptr)
             tcep=>NULL(); select type(uptr); class is(tens_cache_entry_t); tcep=>uptr; end select
             if(associated(tcep)) then
              if(res.eq.GFC_NOT_FOUND) then
               stored=.TRUE.
               call tcep%set_tensor(tensor,errc)
              else
               if(res.ne.GFC_FOUND) errc=-7
              endif
              if(errc.eq.0.and.present(tens_entry_p)) tens_entry_p=>tcep
              tcep=>NULL()
             else
              errc=-6
             endif
             deallocate(tce) !a clone has been stored in this%map
            else
             errc=-5
            endif
           else
            errc=-4
           endif
          else
           errc=-3
          endif
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(errc.ne.0.and.present(tens_entry_p)) tens_entry_p=>NULL()
         if(present(ierr)) ierr=errc
         return
        end function TensCacheStore
!----------------------------------------------------------------
        function TensCacheEvict(this,tensor,ierr) result(evicted)
!Evicts a specific tensor cache entry from the tensor cache.
!If the corresponding tensor cache entry is not found,
!no error is risen, but <evicted>=FALSE.
         implicit none
         logical:: evicted                           !out: TRUE if the tensor cache entry has been found and evicted, FALSE otherwise
         class(tens_cache_t), intent(inout):: this   !inout: tensor cache
         class(tens_rcrsv_t), intent(in):: tensor    !in: tensor to find and evict from the cache
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,res
         type(dictionary_iter_t):: dit
         type(tens_descr_t), target:: tens_descr

         evicted=.FALSE.
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          tens_descr=tensor%get_descriptor(errc)
          if(errc.eq.TEREC_SUCCESS) then
           res=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_tens_descriptors,tens_descr)
           if(res.eq.GFC_FOUND) then; evicted=.TRUE.; else; if(res.ne.GFC_NOT_FOUND) errc=-4; endif
          else
           errc=-3
          endif
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEvict
!-------------------------------------------
        subroutine TensCacheErase(this,ierr)
!Erases the tensor cache completely and unconditionally.
         implicit none
         class(tens_cache_t), intent(inout):: this   !inout: tensor cache
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,res
         type(dictionary_iter_t):: dit

         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          call dit%delete_all(errc); if(errc.ne.GFC_SUCCESS) errc=-3
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheErase
!---------------------------------------
        subroutine tens_cache_dtor(this)
         implicit none
         type(tens_cache_t):: this

         call this%erase()
         return
        end subroutine tens_cache_dtor

       end module virta
