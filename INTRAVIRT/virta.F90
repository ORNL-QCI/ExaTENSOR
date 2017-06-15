!ExaTENSOR: Parallel Virtual Processing for Scale-Adaptive Hierarchical Tensor Algebra:
!Derives from abstract domain-specific virtual processing (DSVP).
!This module provides basic infrastructure for ExaTENSOR, tensor algebra virtual processor (TAVP).
!The logical and numerical tensor algebra virtual processors (L-TAVP, N-TAVP) derive from this module.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/06/15

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
        use talsh                                !on-node heterogeneous numeric tensor algebra
        use pack_prim                            !object packing/unpacking primitives
        use distributed                          !distributed one-sided communication layer
        use service_mpi                          !basic MPI service
        use hardware                             !hardware abstraction
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
 !Errors (ExaTENSOR aliases of DSVP errors):
        integer(INTD), parameter, public:: EXA_SUCCESS=DSVP_SUCCESS                         !success
        integer(INTD), parameter, public:: EXA_ERROR=DSVP_ERROR                             !generic error
        integer(INTD), parameter, public:: EXA_ERR_INVALID_ARGS=DSVP_ERR_INVALID_ARGS       !invalid arguments passed
        integer(INTD), parameter, public:: EXA_ERR_INVALID_REQ=DSVP_ERR_INVALID_REQ         !invalid request
        integer(INTD), parameter, public:: EXA_ERR_MEM_ALLOC_FAIL=DSVP_ERR_MEM_ALLOC_FAIL   !memory allocation failed
        integer(INTD), parameter, public:: EXA_ERR_MEM_FREE_FAIL=DSVP_ERR_MEM_FREE_FAIL     !memory deallocation failed
        integer(INTD), parameter, public:: EXA_ERR_BROKEN_OBJ=DSVP_ERR_BROKEN_OBJ           !broken object
        integer(INTD), parameter, public:: EXA_ERR_UNABLE_COMPLETE=DSVP_ERR_UNABLE_COMPLETE !unable to complete
 !Tensor-algebra virtual processor kinds (roles):
        integer(INTD), parameter, public:: EXA_NO_ROLE=DSVP_NO_KIND !undefined role
        integer(INTD), parameter, public:: EXA_PARSER=0             !domain-specific parser process (not a TAVP)
        integer(INTD), parameter, public:: EXA_WORKER=1             !worker (numeric) process (TAVP)
        integer(INTD), parameter, public:: EXA_MANAGER=2            !manager (logic) process (TAVP)
        integer(INTD), parameter, public:: EXA_HELPER=3             !helper (auxiliary) process (TAVP)
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
          generic, public:: tens_oprnd_ctor=>TensOprndCtor
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty (undefined) state
          final:: tens_oprnd_dtor                               !dtor
        end type tens_oprnd_t
#if 0
 !Tensor instruction control fields:
  !Tensor copy/addition control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_add_t
         type(permutation_t), private:: prmn                    !tensor dimension permutation
         complex(8), private:: alpha                            !alpha prefactor
         contains
          procedure, private:: CtrlTensAddCtor                  !ctor
          generic, public:: ctrl_tens_add_ctor=>CtrlTensAddCtor,CtrlTensAddUnpack
          procedure, public:: pack=>CtrlTensAddPack             !packs the instruction control field into a plain byte packet
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
          procedure, public:: pack=>CtrlTensContrPack           !packs the instruction control field into a plain byte packet
          procedure, public:: unpack=>CtrlTensContrUnpack       !unpacks the instruction control field from a plain byte packet
          final:: ctrl_tens_contr_dtor                          !dtor
        end type ctrl_tens_contr_t
 !Tensor instruction:
        type, extends(ds_instr_t), public:: tens_instr_t
        contains
         procedure, private:: TensInstrCtor                     !ctor
         generic, public:: tens_instr_ctor=>TensInstrCtor
         procedure, public:: decode=>TensInstrDecode            !decoding procedure: Unpacks the raw byte packet (bytecode) and constructs a TAVP instruction
         procedure, public:: encode=>TensInstrEncode            !encoding procedure: Packs the TAVP instruction into a raw byte packet (bytecode)
         final:: tens_instr_dtor                                !dtor
        end type tens_instr_t
#endif
!DATA:
 !MPI process specialization (TAVP role):
        integer(INT_MPI), public:: process_role=EXA_NO_ROLE !TAVP role
        integer(INT_MPI), public:: group_comm=MPI_COMM_NULL !role specific MPI communicator
        integer(INT_MPI), public:: group_size=0             !size of the role specific MPI communicator
        integer(INT_MPI), public:: group_rank=-1            !process rank within the role specific MPI communicator
!VISIBILITY:
 !tens_oprnd_t:
        private TensOprndCtor
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
#endif
       contains
!IMPLEMENTATION:
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
!--------------------------------------------------------
        function TensOprndIsRemote(this,ierr) result(res)
!Returns TRUE if the tensor operand is remote, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INT_MPI):: host_proc_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p

         errc=0
         body_p=>this%tensor%get_body(errc)
         if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
          layout_p=>body_p%get_layout(errc)
          if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
           descr_p=>layout_p%get_data_descr(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
            if(descr_p%is_set(errc,host_proc_rank)) then
             res=.not.(host_proc_rank.eq.impir)
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
         if(present(ierr)) ierr=errc
         return
        end function TensOprndIsRemote
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: buf_size
         integer(INT_MPI):: host_proc_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p

         errc=0
         body_p=>this%tensor%get_body(errc)
         if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
          layout_p=>body_p%get_layout(errc)
          if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
           descr_p=>layout_p%get_data_descr(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
            if(descr_p%is_set(errc,host_proc_rank)) then
             if(host_proc_rank.ne.impir) then !remote tensor operand
              buf_size=descr_p%data_size(errc)
              if(errc.eq.0) then
               !`Allocate local buffer
              else
               errc=-1
              endif
             else
              if(VERBOSE)&
                &write(CONS_OUT,'("WARNING(virta:tens_oprnd_t:acquire_rsc): Attempt to acquire resources for a local operand!")')
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
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndAcquireRsc
!----------------------------------------------
        subroutine TensOprndPrefetch(this,ierr)
!Starts prefetching the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrefetch
!--------------------------------------------
        subroutine TensOprndUpload(this,ierr)
!Starts uploading the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUpload
!---------------------------------------------------------
        function TensOprndSync(this,ierr,wait) result(res)
!Synchronizes the pending prefetch/upload, either TEST or WAIT.
         implicit none
         logical:: res                               !out: TRUE on communication completion, FALSE otherwise
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!---------------------------------------------
        subroutine TensOprndRelease(this,ierr)
!Release local resources occupied by the tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
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
           delivered=this%sync(errc,wait=.TRUE.)
           if(.not.delivered.or.errc.ne.DSVP_SUCCESS) errc=-1
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
           this%tensor=>NULL()
          else
           errc=-3
          endif
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

       end module virta
