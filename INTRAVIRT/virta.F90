!ExaTENSOR: Parallel Virtual Processing for Scale-Adaptive Hierarchical Tensor Algebra:
!Derives from abstract domain-specific virtual processing (DSVP).
!This module provides basic infrastructure for ExaTENSOR, tensor algebra virtual processor (TAVP).
!Different specializations (roles) of tensor algebra virtual processors derive from this module.
!Note that, in general, different specializations (roles) of the domain-specific virtual processor
!may have differing instruction sets (non-overlapping, overlapping, or identical). The instruction
!codes provided in this module are common for all specializations of the tensor algebra virtual processors.
!However, different specializations always have different microcodes, even for the same instruction codes.

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/08/01

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
        use tensor_algebra                       !basic constants
        use talsh                                !on-node heterogeneous numeric tensor algebra
        use pack_prim                            !object packing/unpacking primitives
        use distributed                          !distributed one-sided communication layer
        use service_mpi                          !basic MPI service
        use hardware                             !hardware abstraction and aggregation
        use subspaces                            !hierarchical vector space representation
        use tensor_recursive                     !recursive (hierarchical) tensors
        use dsvp_base                            !abstract domain-specific virtual processor (DSVP)
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
 !Tensor-algebra virtual processor (TAVP) kinds (roles):
        integer(INTD), parameter, public:: EXA_NO_ROLE=DSVP_NO_KIND !undefined role
        integer(INTD), parameter, public:: EXA_DRIVER=0             !domain-specific driver process (not a TAVP)
        integer(INTD), parameter, public:: EXA_MANAGER=1            !manager (logic) process (TAVP)
        integer(INTD), parameter, public:: EXA_WORKER=2             !worker (numeric) process (TAVP)
        integer(INTD), parameter, public:: EXA_HELPER=3             !helper (auxiliary) process (TAVP)
        integer(INTD), public:: EXA_MAX_WORK_GROUP_SIZE=64 !maximal size of a work group (max number of workers per manager)
 !TAVP identification:
        integer(INTD), parameter, public:: TAVP_ANY_ID=-1              !any TAVP
 !TAVP instruction error codes:
        integer(INTD), parameter, public:: TAVP_ERR_GEN_FAILURE=-1     !unspecified generic failure
        integer(INTD), parameter, public:: TAVP_ERR_BTC_BAD=-2         !bad instruction bytecode
        integer(INTD), parameter, public:: TAVP_ERR_CHE_FAILURE=-3     !argument cache failure
        integer(INTD), parameter, public:: TAVP_ERR_ARG_UNDEFINED=-4   !instruction argument is undefined in the argument cache
        integer(INTD), parameter, public:: TAVP_ERR_ARG_DEFINED=-5     !instrcution argument is already defined in the argument cache
        integer(INTD), parameter, public:: TAVP_ERR_RSC_UNAVAILABLE=-6 !instruction is unable to obtain needed resources
        integer(INTD), parameter, public:: TAVP_ERR_COM_FAILURE=-7     !instruction communication failed
        integer(INTD), parameter, public:: TAVP_ERR_EXC_FAILURE=-8     !instruction execution (computation) failed
 !TAVP ISA size (the same for all TAVP specializations):
        integer(INTD), parameter, public:: TAVP_ISA_SIZE=256 !max number of TAVP instructions [0:TAVP_ISA_SIZE-1]
 !Tensor algebra virtual processor (TAVP):
  !TAVP instruction code (opcode), must be non-negative (consult TAProL spec), limited by TAVP_ISA_SIZE:
   !NOOP:
        integer(INTD), parameter, public:: TAVP_INSTR_NOOP=DS_INSTR_NOOP !no operation (empty instruction)
   !General control [0-15]:
        integer(INTD), parameter, public:: TAVP_INSTR_STOP=0      !stop TAVP (finishes current instructions and shutdowns TAVP)
        integer(INTD), parameter, public:: TAVP_INSTR_PAUSE=1     !pause TAVP execution (finishes active instructions and pauses TAVP)
        integer(INTD), parameter, public:: TAVP_INSTR_RESUME=2    !resume TAVP execution (resumes TAVP execution pipeline)
   !Tensor operations [16-255]:
        integer(INTD), parameter, public:: TAVP_INSTR_CREATE=16   !tensor creation
        integer(INTD), parameter, public:: TAVP_INSTR_DESTROY=17  !tensor destruction
        integer(INTD), parameter, public:: TAVP_INSTR_LOAD=18     !tensor loading from persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_SAVE=19     !tensor saving in persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_COMM=20     !tensor communication (get/put/accumulate)
        integer(INTD), parameter, public:: TAVP_INSTR_INIT=21     !tensor initialization (assignement of a value)
        integer(INTD), parameter, public:: TAVP_INSTR_NORM1=22    !tensor 1-norm
        integer(INTD), parameter, public:: TAVP_INSTR_NORM2=23    !tensor 2-norm
        integer(INTD), parameter, public:: TAVP_INSTR_MIN=24      !tensor min element
        integer(INTD), parameter, public:: TAVP_INSTR_MAX=25      !tensor max element
        integer(INTD), parameter, public:: TAVP_INSTR_FOLD=26     !tensor dimension folding
        integer(INTD), parameter, public:: TAVP_INSTR_UNFOLD=27   !tensor dimension unfolding
        integer(INTD), parameter, public:: TAVP_INSTR_SLICE=28    !tensor slicing (taking a slice of a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_INSERT=29   !tensor insertion (inserting a tensor slice in a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_PERMUTE=30  !tensor dimension permutation (in-place)
        integer(INTD), parameter, public:: TAVP_INSTR_COPY=31     !tensor copy (copying a tensor into another tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_SCALE=32    !tensor scaling (multiplication by a number)
        integer(INTD), parameter, public:: TAVP_INSTR_ADD=33      !tensor addition (with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TRACE=34    !tensor trace (tracing over some/all tensor indices)
        integer(INTD), parameter, public:: TAVP_INSTR_PRODUCT=35  !tensor direct product (with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_CONTRACT=36 !tensot contraction (also includes tensor product, tensor addition, and tensor scaling)
!TYPES:
 !Tensor status:
        type, public:: tens_status_t
         logical, public:: created=.FALSE. !TRUE if the tensor has been created (allocated physical memory), FALSE otherwise
         logical, public:: defined=.FALSE. !TRUE if the tensor value has been defined, FALSE otherwise
         logical, public:: replica=.FALSE. !TRUE if the tensor is a replica (another defined instance of this same tensor exists), FALSE otherwise
         logical, public:: is_used=.FALSE. !TRUE if the tensor is currently participating in a computation, FALSE otherwise
         logical, public:: updated=.FALSE. !TRUE if the tensor is currently being updated in a computation, FALSE otherwise
        end type tens_status_t
 !Tensor resource (local resource):
        type, extends(ds_resrc_t), public:: tens_resrc_t
         type(C_PTR), private:: base_addr=C_NULL_PTR   !local buffer address for tensor body storage
         integer(C_SIZE_T), private:: bytes=0_C_SIZE_T !size of the tensor body storage buffer in bytes
         logical, private:: pinned=.FALSE.             !whether or not the buffer is pinned
         integer(C_INT), private:: dev_id=DEV_NULL     !flat device id
         integer(C_INT), private:: ref_count=0         !reference count
         contains
          procedure, public:: is_empty=>TensResrcIsEmpty               !returns TRUE of the tensor resource is empty (unallocated)
          procedure, public:: allocate_buffer=>TensResrcAllocateBuffer !allocates a local buffer for tensor body storage
          procedure, public:: free_buffer=>TensResrcFreeBuffer         !frees the local buffer
          procedure, public:: get_mem_ptr=>TensResrcGetMemPtr          !returns a C pointer to the local memory buffer
          procedure, public:: get_mem_size=>TensResrcGetMemSize        !returns the size of the memory buffer in bytes
          procedure, private:: incr_ref_count=>TensResrcIncrRefCount   !increments the reference count
          procedure, private:: decr_ref_count=>TensResrcDecrRefCount   !decrements the reference count
        end type tens_resrc_t
 !Tensor operand:
        type, extends(ds_oprnd_t), public:: tens_oprnd_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL()   !pointer to a persistent recursive tensor
         class(tens_resrc_t), pointer, private:: resource=>NULL() !pointer to a persistent local tensor resource
         type(talsh_tens_t), private:: talsh_tens                 !TAL-SH tensor object (for performing actual computations)
         contains
          procedure, private:: TensOprndCtor                    !ctor
          generic, public:: tens_oprnd_ctor=>TensOprndCtor
          procedure, public:: set_resource=>TensOprndSetResource!sets the resource component if it has not been set via constructor
          procedure, public:: get_resource=>TensOprndGetResource!returns a pointer to the tensor resource
          procedure, public:: get_tensor=>TensOprndGetTensor    !returns a pointer to the tensor
          procedure, public:: set_talsh_tens=>TensOprndSetTalshTens !sets up the TAL-SH tensor object for further processing with TAL-SH
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty state
          final:: tens_oprnd_dtor                               !dtor
        end type tens_oprnd_t
 !Tensor instruction control fields:
#if 0
  !Tensor addition/copy/slice/insert control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_add_t
         type(permutation_t), private:: prmn                    !tensor dimension permutation
         complex(8), private:: alpha=(1d0,0d0)                  !alpha prefactor
         contains
          procedure, private:: CtrlTensAddCtor                  !ctor
          generic, public:: ctrl_tens_add_ctor=>CtrlTensAddCtor,CtrlTensAddUnpack
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
!DATA:
 !MPI process specialization (TAVP role):
        integer(INT_MPI), public:: process_role=EXA_NO_ROLE !MPI process role (see above)
        integer(INT_MPI), public:: role_comm=MPI_COMM_NULL  !role-specific MPI communicator
        integer(INT_MPI), public:: role_size=0              !size of the role-specific MPI communicator
        integer(INT_MPI), public:: role_rank=-1             !process rank within the role-specific MPI communicator
        integer(INT_MPI), public:: driver_mpi_process=0     !MPI rank of the Driver process (normally process 0 from MPI_COMM_WORLD)
        integer(INT_MPI), allocatable, public:: managers(:) !ordered list of MPI processes serving as managers (ranks in MPI_COMM_WORLD)
 !TAVP address space:
        type(DistrSpace_t), public:: tavp_addr_space        !shared DDSS address space among processes of the same role
!VISIBILITY:
 !non-member:
        public role_barrier
 !tens_resrc_t:
        private TensResrcIsEmpty
        private TensResrcAllocateBuffer
        private TensResrcFreeBuffer
        private TensResrcGetMemPtr
        private TensResrcGetMemSize
        private TensResrcIncrRefCount
        private TensResrcDecrRefCount
 !tens_oprnd_t:
        private TensOprndCtor
        private TensOprndSetResource
        private TensOprndGetResource
        private TensOprndGetTensor
        private TensOprndSetTalshTens
        private TensOprndIsRemote
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndRelease
        private TensOprndDestruct
        public tens_oprnd_dtor
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
!tens_resrc_t]=====================================
        function TensResrcIsEmpty(this) result(ans)
!Returns TRUE if the tensor resource is empty.
         implicit none
         logical:: ans                          !out: answer
         class(tens_resrc_t), intent(in):: this !in: tensor resource

         ans=(this%bytes.le.0_C_SIZE_T)
         return
        end function TensResrcIsEmpty
!---------------------------------------------------------------------------
        subroutine TensResrcAllocateBuffer(this,bytes,ierr,in_buffer,dev_id)
!Allocates local memory either from a system or from a custom buffer.
         implicit none
         class(tens_resrc_t), intent(inout):: this    !inout: tensor resource
         integer(INTL), intent(in):: bytes            !in: size in bytes
         integer(INTD), intent(out), optional:: ierr  !out: error code
         logical, intent(in), optional:: in_buffer    !in: if TRUE the memory will be allocated from a custom buffer, FALSE from the system
         integer(INTD), intent(in), optional:: dev_id !in: flat device id (defaults to Host)
         integer(INTD):: errc
         integer(C_INT):: in_buf,dev
         type(C_PTR):: addr

         errc=0
         if(this%is_empty()) then
          if(bytes.gt.0_INTL) then
           in_buf=NOPE; if(present(in_buffer)) then; if(in_buffer) in_buf=YEP; endif
           dev=talsh_flat_dev_id(DEV_HOST,0); if(present(dev_id)) dev=dev_id
           errc=mem_allocate(dev,int(bytes,C_SIZE_T),in_buf,addr)
           if(errc.eq.0) then
            this%base_addr=addr
            this%bytes=bytes
            this%pinned=(in_buf.ne.NOPE)
            this%dev_id=dev
           else
            errc=-1
           endif
          else
           errc=-2
          endif
         else
          errc=-3
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensResrcAllocateBuffer
!------------------------------------------------
        subroutine TensResrcFreeBuffer(this,ierr)
!Frees the tensor resource buffer.
         implicit none
         class(tens_resrc_t), intent(inout):: this    !inout: tensor resource
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=0
         if(.not.this%is_empty()) then !free only allocated resources
          errc=mem_free(this%dev_id,this%base_addr)
          if(errc.eq.0) then
           this%base_addr=C_NULL_PTR
           this%bytes=0_C_SIZE_T
           this%pinned=.FALSE.
           this%dev_id=DEV_NULL
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensResrcFreeBuffer
!-----------------------------------------------------------------
        function TensResrcGetMemPtr(this,ierr,bytes) result(mem_p)
!Returns a C pointer to the local memory buffer used by the resource.
         implicit none
         type(C_PTR):: mem_p                          !out: C pointer
         class(tens_resrc_t), intent(in):: this       !in: tensor resource
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTL), intent(out), optional:: bytes !out: number of bytes
         integer(INTD):: errc

         errc=0; mem_p=C_NULL_PTR
         if(.not.this%is_empty()) then
          mem_p=this%base_addr
          if(present(bytes)) bytes=this%bytes
         else
          if(present(bytes)) bytes=0_INTL
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensResrcGetMemPtr
!------------------------------------------------------------
        function TensResrcGetMemSize(this,ierr) result(bytes)
!Returns the size of the memory buffer in bytes.
         implicit none
         integer(INTL):: bytes                       !out: size in bytes
         class(tens_resrc_t), intent(in):: this      !in: tensor resource
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(.not.this%is_empty()) then
          bytes=this%bytes; errc=0
         else
          bytes=0_INTL; errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensResrcGetMemSize
!---------------------------------------------
        subroutine TensResrcIncrRefCount(this)
!Increments the reference count.
         implicit none
         class(tens_resrc_t), intent(inout):: this !inout: tensor resource

         this%ref_count=this%ref_count+1
         return
        end subroutine TensResrcIncrRefCount
!---------------------------------------------
        subroutine TensResrcDecrRefCount(this)
!Decrements the reference count.
         implicit none
         class(tens_resrc_t), intent(inout):: this !inout: tensor resource

         this%ref_count=this%ref_count-1
         return
        end subroutine TensResrcDecrRefCount
![tens_oprnd_t]=================================================
        subroutine TensOprndCtor(this,tensor,ierr,tens_resource)
!Constructs a tensor operand. The <tensor> must be set.
!The associated tensor resource is optional and may still be empty.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                            !inout: tensor operand
         class(tens_rcrsv_t), target, intent(in):: tensor                     !in: tensor
         integer(INTD), intent(out), optional:: ierr                          !out: error code
         class(tens_resrc_t), target, intent(inout), optional:: tens_resource !in: local tensor resource, may still be empty
         integer(INTD):: errc

         errc=0
         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%tensor=>tensor
             if(present(tens_resource)) then
              call tens_resource%incr_ref_count()
              this%resource=>tens_resource
             endif
             call this%mark_active(errc)
             if(errc.ne.DSVP_SUCCESS) errc=-1
            else
             errc=-2
            endif
           else
            errc=-3
           endif
          else
           errc=-4
          endif
         else
          errc=-5
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndCtor
!---------------------------------------------------------------
        subroutine TensOprndSetResource(this,tens_resource,ierr)
!Sets the resource component if it has not been set via constructor.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                  !inout: active tensor operand
         class(tens_resrc_t), target, intent(inout):: tens_resource !inout: tensor resource (may be empty)
         integer(INTD), intent(out), optional:: ierr                !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.0) then
           if(.not.associated(this%resource)) then
            call tens_resource%incr_ref_count()
            this%resource=>tens_resource
           else
            errc=-1
           endif
          else
           errc=-2
          endif
         else
          errc=-3
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetResource
!------------------------------------------------------------------
        function TensOprndGetResource(this,ierr) result(resource_p)
!Returns a pointer to the tensor resource.
         implicit none
         type(tens_resrc_t), pointer:: resource_p    !out: pointer to the tensor resource
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; resource_p=>this%resource
         if(.not.associated(resource_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetResource
!------------------------------------------------------------
        function TensOprndGetTensor(this,ierr) result(tens_p)
!Returns a pointer to the tensor.
         implicit none
         type(tens_rcrsv_t), pointer:: tens_p        !out: pointer to the tensor
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; tens_p=>this%tensor
         if(.not.associated(tens_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetTensor
!--------------------------------------------------
        subroutine TensOprndSetTalshTens(this,ierr)
!Sets up the TAL-SH tensor object for further processing with TAL-SH.
!The local tensor body location is imported from the tensor resource.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,nd,unres,lay,data_kind,dims(1:MAX_TENSOR_RANK)
         integer(INTL):: dim_ext(1:MAX_TENSOR_RANK)
         logical:: shpd,layd,locd
         type(C_PTR):: mem_p
         class(tens_rcrsv_t), pointer:: tens_p
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(tens_header_t), pointer:: header_p

         if(this%is_active(errc)) then
          if(errc.eq.0) then
           if(talsh_tensor_is_empty(this%talsh_tens)) then
            tens_p=>this%tensor
            if(tens_p%is_set(errc,num_dims=nd,shaped=shpd,unresolved=unres,layed=layd,located=locd)) then
             if((errc.eq.0).and.(unres.eq.0).and.shpd.and.layd) then
              body_p=>tens_p%get_body(errc)
              if(errc.eq.TEREC_SUCCESS) then
               layout_p=>body_p%get_layout(errc)
               if(errc.eq.TEREC_SUCCESS) then
                lay=layout_p%get_layout_kind(errc)
                if(errc.eq.TEREC_SUCCESS.and.lay.eq.TEREC_LAY_FDIMS) then
                 data_kind=layout_p%get_data_type(errc)
                 if(errc.eq.TEREC_SUCCESS) then
                  if(.not.this%resource%is_empty()) then
                   mem_p=this%resource%get_mem_ptr(errc)
                   if(errc.eq.0) then
                    header_p=>tens_p%get_header(errc)
                    if(errc.eq.TEREC_SUCCESS) then
                     call header_p%get_dims(dim_ext,nd,errc)
                     if(errc.eq.TEREC_SUCCESS) then
                      dims(1:nd)=dim_ext(1:nd)
                      errc=talsh_tensor_construct(this%talsh_tens,data_kind,dims(1:nd),ext_mem=mem_p)
                      if(errc.ne.TALSH_SUCCESS) errc=-14
                     else
                      errc=-13
                     endif
                    else
                     errc=-12
                    endif
                    header_p=>NULL()
                   else
                    errc=-11
                   endif
                  else
                   errc=-10
                  endif
                 else
                  errc=-9
                 endif
                else
                 errc=-8
                endif
               else
                errc=-7
               endif
               layout_p=>NULL()
              else
               errc=-6
              endif
              body_p=>NULL()
             else
              errc=-5
             endif
            else
             errc=-4
            endif
            tens_p=>NULL()
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
        end subroutine TensOprndSetTalshTens
!--------------------------------------------------------
        function TensOprndIsRemote(this,ierr) result(res)
!Returns TRUE if the tensor operand is remote, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INT_MPI):: host_proc_rank,mpi_comm,my_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p

         errc=0
         if(this%is_active()) then
          body_p=>this%tensor%get_body(errc)
          if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
           layout_p=>body_p%get_layout(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
            descr_p=>layout_p%get_data_descr(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
             if(descr_p%is_set(errc,host_proc_rank,mpi_comm)) then
              if(errc.eq.0) then
               call MPI_Comm_Rank(mpi_comm,my_rank,errc)
               if(errc.eq.0) then
                res=.not.(host_proc_rank.eq.my_rank)
               else
                errc=-1
               endif
              else
               errc=-2
              endif
             else
              errc=-3
             endif
            else
             errc=-4
            endif
           else
            errc=-5
           endif
          else
           errc=-6
          endif
         else
          errc=-7
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndIsRemote
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for the remote tensor operand.
!If the resources have already been allocated, does nothing.
!If the resource component is not set, an error will be returned.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand (with an associated resource component)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: buf_size
         integer(INT_MPI):: host_proc_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p

         errc=0
         if(this%is_active()) then
          if(associated(this%resource)) then
           body_p=>this%tensor%get_body(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
            layout_p=>body_p%get_layout(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
             if(layout_p%is_set(errc)) then
              if(errc.eq.TEREC_SUCCESS) then
               buf_size=layout_p%get_body_size(errc)
               if(errc.eq.TEREC_SUCCESS.and.buf_size.gt.0_INTL) then
                if(this%resource%is_empty()) then
                 call this%resource%allocate_buffer(buf_size,errc); if(errc.ne.0) errc=-1
                endif
               else
                errc=-2
               endif
              else
               errc=-3
              endif
             else
              errc=-4
             endif
            else
             errc=-5
            endif
           else
            errc=-6
           endif
          else
           errc=-7
          endif
         else
          errc=-8
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndAcquireRsc
!----------------------------------------------
        subroutine TensOprndPrefetch(this,ierr)
!Starts prefetching the (remote) tensor operand using the local tensor resource.
!If the resource component has not been set, an error will be returned.
!If the local resource has not been allocated, it will be allocated here.
!If the tensor operand has been delivered before, does nothing.
!If there is a pending communication on the tensor operand, returns an error.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand (with an associated resource component)
         integer(INTD), intent(out), optional:: ierr !out: error code, may return TRY_LATER
         integer(INTD):: errc
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         type(C_PTR):: cptr

         if(.not.this%is_present(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%resource)) then
            body_p=>this%tensor%get_body(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
             layout_p=>body_p%get_layout(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
              descr_p=>layout_p%get_data_descr(errc)
              if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
               if(descr_p%is_set(errc)) then
                if(errc.eq.0) then
                 if(this%resource%is_empty()) call this%acquire_rsc(errc)
                 if(errc.eq.0) then
                  if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
                   cptr=this%resource%base_addr
                   call descr_p%get_data(cptr,errc,MPI_ASYNC_REQ)
                   if(errc.ne.0.and.errc.ne.TRY_LATER) errc=-1
                   if(errc.eq.0) then
                    call this%set_comm_stat(DS_OPRND_FETCHING,errc); if(errc.ne.DSVP_SUCCESS) errc=-2
                   endif
                  else
                   errc=-3
                  endif
                 else
                  errc=-4
                 endif
                else
                 errc=-5
                endif
               else
                errc=-6
               endif
              else
               errc=-7
              endif
             else
              errc=-8
             endif
            else
             errc=-9
            endif
           else
            errc=-10
           endif
          else
           errc=-11
          endif
         else
          errc=-12
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrefetch
!--------------------------------------------
        subroutine TensOprndUpload(this,ierr)
!Starts uploading the (remote) tensor operand from the local tensor resource.
!The tensor operand must be marked as delivered (present), even if it is local.
!If there is a pending communication on the tensor operand, returns an error.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         type(C_PTR):: cptr

         if(this%is_present(errc)) then !assumes that the local resource is allocated
          if(errc.eq.DSVP_SUCCESS) then
           body_p=>this%tensor%get_body(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
            layout_p=>body_p%get_layout(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
             descr_p=>layout_p%get_data_descr(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
              if(descr_p%is_set(errc)) then
               if(errc.eq.0) then
                if(.not.this%resource%is_empty()) then !trap
                 if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
                  cptr=this%resource%base_addr
                  call descr_p%acc_data(cptr,errc,MPI_ASYNC_REQ)
                  if(errc.ne.0.and.errc.ne.TRY_LATER) errc=-1
                  if(errc.eq.0) then
                   call this%set_comm_stat(DS_OPRND_UPLOADING,errc); if(errc.ne.DSVP_SUCCESS) errc=-2
                  endif
                 else
                  errc=-3
                 endif
                else
                 errc=-4
                endif
               else
                errc=-5
               endif
              else
               errc=-6
              endif
             else
              errc=-7
             endif
            else
             errc=-8
            endif
           else
            errc=-9
           endif
          else
           errc=-10
          endif
         else
          errc=-11
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUpload
!---------------------------------------------------------
        function TensOprndSync(this,ierr,wait) result(res)
!Synchronizes the pending prefetch/upload, either TEST or WAIT.
!A successful synchronization on prefetch will mark the tensor operand
!as delivered (present). A successful synchronization on upload will
!not change the status of the tensor operand (which is present).
         implicit none
         logical:: res                               !out: TRUE on communication completion, FALSE otherwise
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc,sts
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         logical:: tw

         errc=0; res=.FALSE.
         tw=.FALSE.; if(present(wait)) tw=wait
         if(this%is_active()) then
          body_p=>this%tensor%get_body(errc)
          if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
           layout_p=>body_p%get_layout(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
            descr_p=>layout_p%get_data_descr(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
             if(descr_p%is_set(errc)) then
              if(errc.eq.0) then
               sts=this%get_comm_stat()
               if(sts.ne.DS_OPRND_NO_COMM) then
                if(tw) then
                 call descr_p%wait_data(errc); if(errc.eq.0) then; res=.TRUE.; else; errc=-1; endif
                else
                 res=descr_p%test_data(errc); if(errc.ne.0) errc=-2
                endif
                if(sts.eq.DS_OPRND_FETCHING.and.res) then
                 call this%mark_delivered(errc); if(errc.ne.0) errc=-3
                endif
               else
                errc=-4
               endif
              else
               errc=-5
              endif
             else
              errc=-6
             endif
            else
             errc=-7
            endif
           else
            errc=-8
           endif
          else
           errc=-9
          endif
         else
          errc=-10
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!---------------------------------------------
        subroutine TensOprndRelease(this,ierr)
!Releases local tensor resources occupied by the tensor operand,
!unless there are other active tensor operands sharing the same resource.
!In the latter case, nothing will be done and no error raised.
!Also, if the resource component is not set, nothing will be done either.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: delivered

         errc=0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%resource)) then
            if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
             delivered=this%sync(errc,wait=.TRUE.)
             if(.not.delivered.or.errc.ne.DSVP_SUCCESS) errc=-1
            endif
            if(this%resource%ref_count.eq.1) then !only one (last) tensor operand is associated with this resource
             call this%resource%free_buffer(errc); if(errc.ne.0) errc=-2
            endif
           endif
          else
           errc=-3
          endif
         else
          if(errc.ne.DSVP_SUCCESS) errc=-4
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndRelease
!----------------------------------------------
        subroutine TensOprndDestruct(this,ierr)
!Destructs the tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         logical:: delivered

         errc=0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if(.not.delivered.or.errc.ne.DSVP_SUCCESS) errc=-1
           endif
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
           if(.not.talsh_tensor_is_empty(this%talsh_tens)) then
            ier=talsh_tensor_destruct(this%talsh_tens)
            if(ier.ne.TALSH_SUCCESS.and.errc.eq.0) errc=-3
           endif
           if(associated(this%resource)) then
            call this%resource%decr_ref_count()
            this%resource=>NULL()
           endif
           this%tensor=>NULL()
          else
           errc=-4
          endif
         else
          if(errc.ne.DSVP_SUCCESS) errc=-5
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndDestruct
!---------------------------------------
        subroutine tens_oprnd_dtor(this)
         implicit none
         type(tens_oprnd_t):: this
         integer(INTD):: errc

         call this%destruct(errc)
         if(errc.ne.0) call quit(errc,'#FATAL(virta:tens_oprnd_dtor): Destructor failed!')
         return
        end subroutine tens_oprnd_dtor
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

       end module virta
