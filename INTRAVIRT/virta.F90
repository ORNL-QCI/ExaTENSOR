!ExaTENSOR: Parallel Virtual Processing for Scale-Adaptive Hierarchical Tensor Algebra:
!Derives from abstract domain-specific virtual processing (DSVP).
!This module provides basic infrastructure for ExaTENSOR, tensor algebra virtual processor (TAVP).
!Different specializations (roles) of tensor algebra virtual processors derive from this module.
!Note that, in general, different specializations (roles) of the domain-specific virtual processor
!may have differing instruction sets (non-overlapping, overlapping, or identical). The instruction
!codes provided in this module are common for all specializations of the tensor algebra virtual processors.
!However, different specializations always have different microcodes, even for the same instruction codes.

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/11/17

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
        use gfc_vector         !GFC vector
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
        integer(INTD), public:: EXA_MAX_WORK_GROUP_SIZE=2 !maximal size of a work group (max number of workers per manager)
        integer(INTD), public:: EXA_MANAGER_BRANCH_FACT=2 !branching factor for the managing hierarchy
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
 !Tensor algebra virtual processor (TAVP) ISA:
  !General:
        integer(INTD), parameter, public:: TAVP_ISA_SIZE=256       !max number of TAVP instruction codes [0:TAVP_ISA_SIZE-1]
        integer(INTD), parameter, public:: TAVP_ISA_CTRL_FIRST=0   !first TAVP_INSTR_CTRL_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_CTRL_LAST=15   !last TAVP_INSTR_CTRL_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_SPACE_FIRST=16 !first TAVP_INSTR_SPACE_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_SPACE_LAST=63  !last TAVP_INSTR_SPACE_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_TENS_FIRST=64  !first TAVP_INSTR_TENS_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_TENS_LAST=255  !last TAVP_INSTR_TENS_XXX code
  !TAVP instruction code (opcode), must be non-negative (consult TAProL spec), limited by TAVP_ISA_SIZE:
   !NOOP:
        integer(INTD), parameter, public:: TAVP_INSTR_NOOP=DS_INSTR_NOOP !no operation (empty instruction): Negative opcode
   !General control [0-15]:
        integer(INTD), parameter, public:: TAVP_INSTR_CTRL_STOP=0       !stop TAVP (finishes current instructions and shutdowns TAVP)
        integer(INTD), parameter, public:: TAVP_INSTR_CTRL_RESUME=1     !resume TAVP execution (resumes TAVP execution pipeline after a pause, ignored if no pause has been posted)
        integer(INTD), parameter, public:: TAVP_INSTR_CTRL_PAUSE=2      !pause TAVP execution (finishes active instructions and pauses TAVP)
   !Auxiliary definitions [16-63]:
        integer(INTD), parameter, public:: TAVP_INSTR_SPACE_CREATE=16   !create a (hierarchical) vector space
        integer(INTD), parameter, public:: TAVP_INSTR_SPACE_DESTROY=17  !destroy a (hierarchical) vector space
   !Tensor operations (decomposable) [64-255]:
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_CREATE=64    !tensor creation
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_DESTROY=65   !tensor destruction
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_LOAD=66      !tensor loading from persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_SAVE=67      !tensor saving in persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_COMM=68      !tensor communication (get/put/accumulate)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_INIT=69      !tensor initialization (assignement of a value)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_NORM1=70     !tensor 1-norm
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_NORM2=71     !tensor 2-norm
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_MIN=72       !tensor min element
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_MAX=73       !tensor max element
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_FOLD=74      !tensor dimension folding
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_UNFOLD=75    !tensor dimension unfolding
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_SLICE=76     !tensor slicing (taking a slice of a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_INSERT=77    !tensor insertion (inserting a tensor slice in a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_PERMUTE=78   !tensor dimension permutation (in-place)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_COPY=79      !tensor copy (copying a tensor into another tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_SCALE=80     !tensor scaling (multiplication by a number)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_ADD=81       !tensor addition (with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_TRACE=82     !tensor trace (tracing over some/all tensor indices)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_PRODUCT=83   !tensor direct product (with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_CONTRACT=84  !tensor contraction (also includes tensor product, tensor addition, and tensor scaling)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_REPLICATE=85 !tensor replication
!TYPES:
 !Tensor status:
        type, public:: tens_status_t
         logical, public:: created=.FALSE. !TRUE if the tensor has been created (allocated physical memory), FALSE otherwise
         logical, public:: defined=.FALSE. !TRUE if the tensor value has been defined, FALSE otherwise
         logical, public:: replica=.FALSE. !TRUE if the tensor is a replica (another defined instance of this same tensor exists), FALSE otherwise
         logical, public:: updated=.FALSE. !TRUE if the tensor value is currently being updated, FALSE otherwise
         integer(INTD), public:: is_used=0 !number of read-only references which are currently using the tensor
        end type tens_status_t
        type(tens_status_t), parameter:: TENS_STATUS_NONE=tens_status_t() !tensor status null
        public TENS_STATUS_NONE
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
         integer(INTD), private:: ref_count=0                   !reference count: Number of existing tensor operands associated with the tensor cache entry
         contains
          procedure, public:: set_tensor=>TensCacheEntrySetTensor         !sets the pointer to a tensor
          procedure, public:: is_set=>TensCacheEntryIsSet                 !returns TRUE if the tensor cache entry is set (constructed)
          procedure, public:: get_tensor=>TensCacheEntryGetTensor         !returns a non-owning pointer to the tensor
          procedure, public:: get_status=>TensCacheEntryGetStatus         !returns the tensor status object
          procedure, public:: mark_created=>TensCacheEntryMarkCreated     !marks the tensor status as created (allocated memory)
          procedure, public:: mark_defined=>TensCacheEntryMarkDefined     !marks the tensor status as defined to some value
          procedure, public:: mark_in_use=>TensCacheEntryMarkInUse        !marks the tensor status as in-use (read-only) and increments the reference count
          procedure, public:: mark_no_use=>TensCacheEntryMarkNoUse        !decrements the read-only usage reference count
          procedure, public:: mark_updated=>TensCacheEntryMarkUpdated     !marks the tensor status as in-update (currently being updated)
          procedure, public:: mark_empty=>TensCacheEntryMarkEmpty         !marks the tensor status as empty (destroyed)
          procedure, public:: incr_ref_count=>TensCacheEntryIncrRefCount  !increments the reference count
          procedure, public:: decr_ref_count=>TensCacheEntryDecrRefCount  !decrements the reference count
          procedure, public:: get_ref_count=>TensCacheEntryGetRefCount    !returns the current reference count: Number of existing tensor operands associated with the tensor cache entry
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
 !External data register:
        type, public:: data_register_t
         type(dictionary_t), private:: ext_data                           !string --> talsh_tens_data_t
         contains
          procedure, public:: register_data=>DataRegisterRegisterData     !registers external local data
          procedure, public:: unregister_data=>DataRegisterUnregisterData !unregisters previously registered data
          procedure, public:: retrieve_data=>DataRegisterRetrieveData     !retrieves previously registered data
          final:: data_register_dtor                                      !dtor
        end type data_register_t
 !External method register:
        type, public:: method_register_t
         type(dictionary_t), private:: ext_methods                              !string --> talsh_tens_definer_t
         contains
          procedure, public:: register_method=>MethodRegisterRegisterMethod     !registers an external tensor defining procedure
          procedure, public:: unregister_method=>MethodRegisterUnregisterMethod !unregisters a previously registered tensor defining procedure
          procedure, public:: retrieve_method=>MethodRegisterRetrieveMethod     !retrieves a previously registered tensor defining procedure
          final:: method_register_dtor                                          !dtor
        end type method_register_t
!INTERFACES:
        abstract interface
 !tens_cache_entry_t:
  !Polymorphic allocator:
         function tens_cache_entry_alloc_i(tens_cache_entry) result(ierr)
          import:: tens_cache_entry_t,INTD
          integer(INTD):: ierr
          class(tens_cache_entry_t), allocatable, intent(out):: tens_cache_entry
         end function tens_cache_entry_alloc_i
 !External method (user-defined tensor operation):
         function exatns_method_i(tens_args,scal_args) result(ierr)
          import:: INTD,tens_rcrsv_t
          implicit none
          integer(INTD):: ierr                                !out: error code
          class(tens_rcrsv_t), intent(inout):: tens_args(0:)  !inout: tensor arguments
          complex(8), intent(inout), optional:: scal_args(0:) !inout: scalar arguments
         end function exatns_method_i
        end interface
        public tens_cache_entry_alloc_i
        public exatns_method_i
!DATA:
 !MPI process specialization (TAVP role, set by exatns_start):
        integer(INT_MPI), public:: process_role=EXA_NO_ROLE   !MPI process role (see above)
        integer(INT_MPI), public:: role_comm=MPI_COMM_NULL    !role-specific MPI communicator
        integer(INT_MPI), public:: role_size=0                !size of the role-specific MPI communicator
        integer(INT_MPI), public:: role_rank=-1               !process rank within the role-specific MPI communicator
        integer(INT_MPI), public:: driver_gl_rank=-1          !driver process rank in the global MPI communicator
        integer(INT_MPI), public:: top_manager_gl_rank=-1     !root manager process rank in the global MPI communicator
        integer(INT_MPI), public:: drv_mng_comm=MPI_COMM_NULL !MPI intercommunicator for the driver and managers
        integer(INT_MPI), public:: mng_wrk_comm=MPI_COMM_NULL !MPI intercommunicator for managers and workers
 !External universal tensor dimension strength assessing function:
        procedure(tens_rcrsv_dim_strength_i), pointer, public:: dim_strength_assess=>NULL() !assesses the strength of tensor dimensions
        real(8), public:: dim_strength_thresh=0d0         !tensor dimension strength threshold above which the dimension will split
 !External data register:
        type(data_register_t), public:: data_register     !string --> talsh_tens_data_t
 !External method register:
        type(method_register_t), public:: method_register !string --> talsh_tens_definer_t
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
        private TensCacheEntryIncrRefCount
        private TensCacheEntryDecrRefCount
        private TensCacheEntryGetRefCount
        private TensCacheEntryNullifyTensor
 !tens_cache_t:
        private TensCacheLookup
        private TensCacheStore
        private TensCacheEvict
        private TensCacheErase
        public tens_cache_dtor
 !data_register_t:
        private DataRegisterRegisterData
        private DataRegisterUnregisterData
        private DataRegisterRetrieveData
        public data_register_dtor
 !method_register_t:
        private MethodRegisterRegisterMethod
        private MethodRegisterUnregisterMethod
        private MethodRegisterRetrieveMethod
        public method_register_dtor

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
!--------------------------------------------------
        subroutine TensCacheEntryIncrRefCount(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry

!$OMP ATOMIC UPDATE
         this%ref_count=this%ref_count+1
         return
        end subroutine TensCacheEntryIncrRefCount
!--------------------------------------------------
        subroutine TensCacheEntryDecrRefCount(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry

!$OMP ATOMIC UPDATE
         this%ref_count=this%ref_count-1
         return
        end subroutine TensCacheEntryDecrRefCount
!----------------------------------------------------------------
        function TensCacheEntryGetRefCount(this) result(refcount)
         implicit none
         integer(INTD):: refcount                     !out: reference count
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry

!$OMP ATOMIC READ
         refcount=this%ref_count
         return
        end function TensCacheEntryGetRefCount
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
         tens_descr=tensor%get_descriptor(errc,skip_body=.TRUE.,skip_location=.TRUE.)
         if(errc.eq.TEREC_SUCCESS) then
!!!$OMP CRITICAL (TAVP_CACHE)
          errc=dit%init(this%map)
          if(errc.eq.GFC_SUCCESS) then
           uptr=>NULL()
           res=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_tens_descriptors,tens_descr,value_out=uptr)
           if(res.eq.GFC_FOUND) then
            select type(uptr); class is(tens_cache_entry_t); tens_entry_p=>uptr; end select
            if(.not.associated(tens_entry_p)) errc=-5 !trap
           else
            if(res.ne.GFC_NOT_FOUND) errc=-4
           endif
           res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
          else
           errc=-2
          endif
!!!$OMP END CRITICAL (TAVP_CACHE)
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
         if(associated(tensor)) then
          tens_descr=tensor%get_descriptor(errc,skip_body=.TRUE.,skip_location=.TRUE.)
          if(errc.eq.TEREC_SUCCESS) then
!!!$OMP CRITICAL (TAVP_CACHE)
           errc=dit%init(this%map)
           if(errc.eq.GFC_SUCCESS) then
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
            res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-4
           else
            errc=-3
           endif
!!!$OMP END CRITICAL (TAVP_CACHE)
          else
           errc=-2
          endif
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
         tens_descr=tensor%get_descriptor(errc,skip_body=.TRUE.,skip_location=.TRUE.)
         if(errc.eq.TEREC_SUCCESS) then
!!!$OMP CRITICAL (TAVP_CACHE)
          errc=dit%init(this%map)
          if(errc.eq.GFC_SUCCESS) then
           res=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_tens_descriptors,tens_descr)
           if(res.eq.GFC_FOUND) then; evicted=.TRUE.; else; if(res.ne.GFC_NOT_FOUND) errc=-4; endif
           res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
          else
           errc=-2
          endif
!!!$OMP END CRITICAL (TAVP_CACHE)
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

!!!$OMP CRITICAL (TAVP_CACHE)
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%delete_all(); if(errc.ne.GFC_SUCCESS) errc=-3
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
!!!$OMP END CRITICAL (TAVP_CACHE)
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
![data_register_t]---------------------------------------------------------
        subroutine DataRegisterRegisterData(this,data_name,extrn_data,ierr)
         implicit none
         class(data_register_t), intent(inout):: this     !inout: data register
         character(*), intent(in):: data_name             !in: data name
         type(talsh_tens_data_t), intent(in):: extrn_data !in: external data
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_strings,data_name,extrn_data)
          if(errc.eq.GFC_NOT_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DataRegisterRegisterData
!-----------------------------------------------------------------
        subroutine DataRegisterUnregisterData(this,data_name,ierr)
         implicit none
         class(data_register_t), intent(inout):: this !inout: data register
         character(*), intent(in):: data_name         !in: data name
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_strings,data_name)
          if(errc.eq.GFC_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DataRegisterUnregisterData
!--------------------------------------------------------------------------------
        function DataRegisterRetrieveData(this,data_name,ierr) result(extrn_data)
         implicit none
         type(talsh_tens_data_t), pointer:: extrn_data !out: pointer to the data
         class(data_register_t), intent(in):: this     !inout: data register
         character(*), intent(in):: data_name          !in: data name
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit
         class(*), pointer:: uptr

         extrn_data=>NULL()
         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>NULL()
          errc=dit%search(GFC_DICT_JUST_FIND,cmp_strings,data_name,value_out=uptr)
          if(errc.eq.GFC_FOUND) then
           if(associated(uptr)) then
            select type(uptr); type is(talsh_tens_data_t); extrn_data=>uptr; end select
            if(.not.associated(extrn_data)) errc=-6
           else
            errc=-5
           endif
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function DataRegisterRetrieveData
!------------------------------------------
        subroutine data_register_dtor(this)
         implicit none
         type(data_register_t):: this  !inout: data register
         integer(INTD):: errc
         type(dictionary_iter_t):: dit

         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%delete_all()
          errc=dit%release()
         endif
         return
        end subroutine data_register_dtor
![method_register_t]---------------------------------------------------------------
        subroutine MethodRegisterRegisterMethod(this,method_name,extrn_method,ierr)
         implicit none
         class(method_register_t), intent(inout):: this         !inout: method register
         character(*), intent(in):: method_name                 !in: method name
         class(talsh_tens_definer_t), intent(in):: extrn_method !in: external tensor defining method
         integer(INTD), intent(out), optional:: ierr            !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_strings,method_name,extrn_method)
          if(errc.eq.GFC_NOT_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine MethodRegisterRegisterMethod
!-----------------------------------------------------------------------
        subroutine MethodRegisterUnregisterMethod(this,method_name,ierr)
         implicit none
         class(method_register_t), intent(inout):: this !inout: method register
         character(*), intent(in):: method_name         !in: method name
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_strings,method_name)
          if(errc.eq.GFC_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine MethodRegisterUnregisterMethod
!----------------------------------------------------------------------------------------
        function MethodRegisterRetrieveMethod(this,method_name,ierr) result(extrn_method)
         implicit none
         class(talsh_tens_definer_t), pointer:: extrn_method !out: pointer to the method
         class(method_register_t), intent(in):: this         !inout: method register
         character(*), intent(in):: method_name              !in: method name
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit
         class(*), pointer:: uptr

         extrn_method=>NULL()
         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>NULL()
          errc=dit%search(GFC_DICT_JUST_FIND,cmp_strings,method_name,value_out=uptr)
          if(errc.eq.GFC_FOUND) then
           if(associated(uptr)) then
            select type(uptr); class is(talsh_tens_definer_t); extrn_method=>uptr; end select
            if(.not.associated(extrn_method)) errc=-6
           else
            errc=-5
           endif
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function MethodRegisterRetrieveMethod
!--------------------------------------------
        subroutine method_register_dtor(this)
         implicit none
         type(method_register_t):: this  !inout: method register
         integer(INTD):: errc
         type(dictionary_iter_t):: dit

         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%delete_all()
          errc=dit%release()
         endif
         return
        end subroutine method_register_dtor

       end module virta
