!ExaTENSOR: Parallel Virtual Processing for Scale-Adaptive Hierarchical Tensor Algebra:
!Derives from abstract domain-specific virtual processing (DSVP).
!This module provides basic infrastructure for ExaTENSOR, tensor algebra virtual processor (TAVP).
!The logical and numerical tensor algebra virtual processors (L-TAVP, N-TAVP) derive from this module.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/05/31

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
        use dil_basic                            !basic constants
        use talsh                                !on-node numeric tensor algebra
        use pack_prim                            !object packing primitives
        use distributed                          !distributed communication layer
#ifndef NO_LINUX
        use service_mpi, only: get_memory_status,MPI_COMM_NULL
#else
        use service_mpi, only: MPI_COMM_NULL
#endif
        use hardware                             !hardware abstraction
        use subspaces                            !hierarchical vector space representation
        use tensor_recursive                     !recursive tensors
        use dsvp_base                            !abstract domain-specific virtual processor (DSVP)
        implicit none
        public
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Errors (ExaTENSOR aliases of DSVP errors):
        integer(INTD), parameter, public:: EXA_SUCCESS=DSVP_SUCCESS                       !success
        integer(INTD), parameter, public:: EXA_ERROR=DSVP_ERROR                           !generic error
        integer(INTD), parameter, public:: EXA_ERR_INVALID_ARGS=DSVP_ERR_INVALID_ARGS     !invalid arguments passed to a procedure
        integer(INTD), parameter, public:: EXA_ERR_INVALID_REQ=DSVP_ERR_INVALID_REQ       !invalid request
        integer(INTD), parameter, public:: EXA_ERR_MEM_ALLOC_FAIL=DSVP_ERR_MEM_ALLOC_FAIL !memory allocation failed
        integer(INTD), parameter, public:: EXA_ERR_MEM_FREE_FAIL=DSVP_ERR_MEM_FREE_FAIL   !memory deallocation failed
        integer(INTD), parameter, public:: EXA_ERR_BROKEN_OBJ=DSVP_ERR_BROKEN_OBJ         !broken object
 !Tensor-algebra virtual processor kinds (roles):
        integer(INTD), parameter, public:: EXA_NO_ROLE=0   !undefined role
        integer(INTD), parameter, public:: EXA_MANAGER=1   !manager process (global root is a manager as well)
        integer(INTD), parameter, public:: EXA_WORKER=2    !worker process (aka C-process)
        integer(INTD), parameter, public:: EXA_HELPER=3    !helper process
        integer(INTD), public:: EXA_MAX_WORK_GROUP_SIZE=64 !maximal size of a work group (max number of workers per manager)
 !Elementary tensor instruction (ETI) granularity classification:
        real(8), public:: EXA_FLOPS_MEDIUM=1d10 !minimal number of Flops to consider the operation as medium-cost
        real(8), public:: EXA_FLOPS_HEAVY=1d12  !minimal number of Flops to consider the operation as heavy-cost
        real(8), public:: EXA_COST_TO_SIZE=1d2  !minimal cost (Flops) to size (Words) ratio to consider the operation compute intensive
 !Tensor algebra virtual processor (TAVP):
  !Tensor instruction code (opcode), must be non-negative (consult TAProL spec):
        integer(INTD), parameter, public:: TENSOR_INSTR_NOOP=DS_INSTR_NOOP !no operation (empty instruction)
        integer(INTD), parameter, public:: TENSOR_INSTR_COMM=0      !tensor communication (get/put/accumulate)
        integer(INTD), parameter, public:: TENSOR_INSTR_INIT=1      !tensor initialization (assignement to a value)
        integer(INTD), parameter, public:: TENSOR_INSTR_NORM1=2     !tensor 1-norm
        integer(INTD), parameter, public:: TENSOR_INSTR_NORM2=3     !tensor 2-norm
        integer(INTD), parameter, public:: TENSOR_INSTR_MIN=4       !tensor min element
        integer(INTD), parameter, public:: TENSOR_INSTR_MAX=5       !tensor max element
        integer(INTD), parameter, public:: TENSOR_INSTR_FOLD=6      !tensor dimension folding
        integer(INTD), parameter, public:: TENSOR_INSTR_UNFOLD=7    !tensor dimension unfolding
        integer(INTD), parameter, public:: TENSOR_INSTR_SCALE=8     !tensor scaling (multiplication by a number)
        integer(INTD), parameter, public:: TENSOR_INSTR_SLICE=9     !tensor slicing (taking a slice of a tensor)
        integer(INTD), parameter, public:: TENSOR_INSTR_INSERT=10   !tensor insertion (inserting a tensor slice in a tensor)
        integer(INTD), parameter, public:: TENSOR_INSTR_PERMUTE=11  !tensor dimension permutation (in-place)
        integer(INTD), parameter, public:: TENSOR_INSTR_COPY=12     !tensor copy (copying a tensor into another tensor with an optional index permutation)
        integer(INTD), parameter, public:: TENSOR_INSTR_ADD=13      !tensor addition
        integer(INTD), parameter, public:: TENSOR_INSTR_TRACE=14    !tensor trace (tracing over some/all tensor indices)
        integer(INTD), parameter, public:: TENSOR_INSTR_PRODUCT=15  !tensor direct (Cartesian) product
        integer(INTD), parameter, public:: TENSOR_INSTR_CONTRACT=16 !tensot contraction (also includes tensor product, tensor addition, and tensor scaling)
        integer(INTD), parameter, public:: TENSOR_ISA_SIZE=17       !total number of non-negative tensor instruction codes
!TYPES:
 !Tensor operand:
        type, extends(ds_oprnd_t), public:: tens_oprnd_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL() !pointer to a recursive tensor
         contains
          procedure, private:: TensOprndCtor                    !ctor
          generic, public:: tens_oprnd_ctor=>TensOprndCtor,TensOprndUnpack
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: pack=>TensOprndPack               !packs the specification of the tensor operand into a plain byte packet
          procedure, public:: unpack=>TensOprndUnpack           !unpacks the specification of the tensor operand from a plain byte packet
          procedure, public:: destruct=>TensOprndDestruct       !dtor
        end type tens_oprnd_t
 !Tensor instruction control fields:
  !Tensor copy/addition control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_add_t
         type(permutation_t), private:: prmn                    !tensor dimension permutation
         complex(8), private:: alpha                            !alpha prefactor
         contains
          procedure, private:: CtrlTensAddCtor                  !ctor
          generic, public:: ctrl_tens_add_ctor=>CtrlTensAddCtor,CtrlTensAddUnpack
          procedure, public:: pack=>CtrlTensAddPack             !pack the instruction control field in a plain byte packet
          procedure, public:: unpack=>CtrlTensAddUnpack         !unpacks the instruction control field from a plain byte packet
          final:: ctrl_tens_add_dtor                            !dtor
        end type ctrl_tens_add_t
  !Tensor contraction control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_contr_t
         type(contr_ptrn_ext_t), private:: contr_ptrn           !extended tensor contraction pattern
         complex(8), private:: alpha                            !alpha prefactor
         contains
          procedure, private:: CtrlTensContrCtor                !ctor
          generic, public:: ctrl_tens_contr_ctor=>CtrlTensContrCtor,CtrlTensContrUnpack
          procedure, public:: pack=>CtrlTensContrPack           !pack the instruction control field in a plain byte packet
          procedure, public:: unpack=>CtrlTensContrUnpack       !unpacks the instruction control field from a plain byte packet
          final:: ctrl_tens_contr_dtor                          !dtor
        end type ctrl_tens_contr_t
 !Tensor instruction:
        type, extends(ds_instr_t), public:: tens_instr_t
        contains
         procedure, private:: TensInstrCtor                     !ctor
         generic, public:: tens_instr_ctor=>TensInstrCtor
         procedure, public:: decode=>TensInstrDecode            !decoding procedure: Unpacks the raw byte packet and constructs a TAVP instruction
         procedure, public:: encode=>TensInstrEncode            !encoding procedure: Packs the TAVP instruction into a raw byte packet
         final:: tens_instr_dtor                                !dtor
        end type tens_instr_t
!DATA`Remove:
 !Current role of the tensor alegbra virtual processor:
        integer(INTD), protected:: my_role=EXA_NO_ROLE         !role of this virtual processor (set at run-time)
        integer(INTD), protected:: my_group=-1                 !computing group the virtual processor belongs to (set at run-time): [0..MAX]
        integer(INTD), protected:: my_group_size=0             !size of the computing group the virtual processor belongs to (set at runtime): [1..EXA_MAX_WORK_GROUP_SIZE]
        integer(INTD), protected:: my_group_index=-1           !virtual processor ID within its computing group (set at run-time): [0..my_group_size-1]
        integer(INTD), protected:: my_group_comm=MPI_COMM_NULL !group MPI communicator (if any)
!VISIBILITY:
 !non-member:
        public tavp_establish_role !`Remove
 !tens_oprnd_t:
        private TensOprndCtor
        private TensOprndIsRemote
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndRelease
        private TensOprndPack
        private TensOprndUnpack
        private TensOprndDestruct
 !ctrl_tens_add_t:
        private CtrlTensAddCtor
        private CtrlTensAddPack
        private CtrlTensAddUnpack
        public ctrl_tens_add_dtor
 !ctrl_tens_contr_t:
        private CtrlTensContrCtor
        private CtrlTensContrPack
        private CtrlTensContrUnpack
        public ctrl_tens_contr_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrDecode
        private TensInstrEncode
        public tens_instr_dtor

       contains
!IMPLEMENTATION:
![Non-member]====================================
        subroutine tavp_establish_role(role,ierr)
!Specializes a tensor algebra virtual processor (TAVP).
         implicit none
         integer(INTD), intent(in):: role  !in: specific role of the tensor algebra virtual processor
         integer(INTD), intent(out):: ierr !out: error code

         ierr=EXA_SUCCESS
         select case(role)
         case(EXA_MANAGER,EXA_WORKER,EXA_HELPER)
          my_role=role
         case(EXA_NO_ROLE)
          ierr=EXA_ERROR
         case default
          ierr=EXA_ERR_INVALID_ARGS
         end select
         return
        end subroutine tavp_establish_role
![tens_oprnd_t]===================================
        subroutine TensOprndCtor(this,tensor,ierr)
!Constructs a tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this         !inout: tensor operand
         class(tens_rcrsv_t), pointer, intent(in):: tensor !in: tensor
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=0
         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(tensor)) then
            if(tensor%is_set(errc)) then
             if(errc.eq.TEREC_SUCCESS) then
              this%tensor=>tensor
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
         else
          errc=-6
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndCtor
!-------------------------------------------------
        subroutine TensOprndPack(this,packet,ierr)
!Packs the tensor operand into a packet.
         implicit none
         class(tens_oprnd_t), intent(in):: this      !in: tensor operand
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%tensor)) then
            call this%tensor%pack(packet,errc)
            if(errc.ne.TEREC_SUCCESS) errc=-1
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
        end subroutine TensOprndPack
!---------------------------------------------------
        subroutine TensOprndUnpack(this,packet,ierr)
!Unpacks the tensor operand from a packet (ctor).
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         class(obj_pack_t), intent(inout):: packet   !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUnpack
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
           delivered=this%sync(errc,wait=.TRUE.)
           if(.not.delivered.or.errc.ne.DSVP_SUCCESS) errc=-1
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
          else
           errc=-3
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndDestruct

       end module virta
