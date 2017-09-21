!ExaTENSOR: TAVP-Manager (TAVP-MNG) implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/09/21

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

       module tavp_manager
        use virta
        use gfc_base
        use gfc_list
        use gfc_dictionary
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output device
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
!TYPES:
 !Tensor argument cache entry (TAVP-specific):
        type, extends(tens_cache_entry_t), private:: tens_entry_mng_t
         integer(INTD), private:: owner_id=-1                         !tensor owner id (non-negative TAVP-MNG id)
         contains
          procedure, private:: TensEntryMngCtor                       !ctor
          generic, public:: tens_entry_mng_ctor=>TensEntryMngCtor
          procedure, public:: get_owner_id=>TensEntryMngGetOwnerId    !returns the owner id
          final:: tens_entry_mng_dtor
        end type tens_entry_mng_t
 !Tensor operand (encapsulated tensor data processible by a specific TAVP):
        type, extends(ds_oprnd_t), private:: tens_oprnd_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL() !non-owning pointer to a persistent recursive tensor
         integer(INTD), private:: owner_id=-1                   !tensor owner id
         contains
          procedure, private:: TensOprndCtor                    !ctor
          generic, public:: tens_oprnd_ctor=>TensOprndCtor
          procedure, public:: get_tensor=>TensOprndGetTensor    !returns a pointer to the tensor
          procedure, public:: get_owner_id=>TensOprndGetOwnerId !returns the tensor owner id
          procedure, public:: set_owner_id=>TensOprndSetOwnerId !sets the tensor owner id
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty state
          final:: tens_oprnd_dtor
        end type tens_oprnd_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), private:: tens_instr_t
         integer(INTD), private:: num_out_oprnds=0                    !number of the output tensor instruction operands
         integer(INTD), private:: out_oprnds(0:MAX_TENSOR_OPERANDS-1) !tensor instruction operands which are considered output (to sync on them for completion)
         contains
          procedure, private:: TensInstrCtor               !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: encode=>TensInstrEncode      !encoding procedure: Packs the TAVP instruction into a raw byte packet
          final:: tens_instr_dtor                          !dtor
        end type tens_instr_t
 !TAVP-MNG decoder:
        type, extends(ds_decoder_t), private:: tavp_mng_decoder_t
         integer(INTD), private:: source_comm                       !bytecode source communicator
         integer(INTD), private:: source_rank=-1                    !bytecode source process rank
         type(pack_env_t), private:: bytecode                       !incoming bytecode
         contains
          procedure, public:: configure=>TAVPMNGDecoderConfigure    !configures TAVP-MNG decoder
          procedure, public:: start=>TAVPMNGDecoderStart            !starts TAVP-MNG decoder
          procedure, public:: shutdown=>TAVPMNGDecoderShutdown      !shuts down TAVP-MNG decoder
          procedure, public:: decode=>TAVPMNGDecoderDecode          !decodes the DS bytecode into DS instructions
        end type tavp_mng_decoder_t
 !TAVP-MNG decoder configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_decoder_conf_t
         integer(INTD), public:: source_comm !MPI communicator of the source process
         integer(INTD), public:: source_rank !source process rank from which the bytecode is coming
        end type tavp_mng_decoder_conf_t
 !TAVP-MNG retirer:
        type, extends(ds_encoder_t), private:: tavp_mng_retirer_t
         integer(INTD), private:: retire_comm                       !retired bytecode destination communicator
         integer(INTD), private:: retire_rank=-1                    !retired bytecode destination process rank
         type(pack_env_t), private:: bytecode                       !outgoing bytecode
         contains
          procedure, public:: configure=>TAVPMNGRetirerConfigure    !configures TAVP-MNG retirer
          procedure, public:: start=>TAVPMNGRetirerStart            !starts TAVP-MNG retirer
          procedure, public:: shutdown=>TAVPMNGRetirerShutdown      !shuts down TAVP-MNG retirer
          procedure, public:: encode=>TAVPMNGRetirerEncode          !encodes a DS instruction into the DS bytecode
        end type tavp_mng_retirer_t
 !TAVP-MNG retirer configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_retirer_conf_t
         integer(INTD), public:: retire_comm                        !MPI communicator of the retired bytecode destination process
         integer(INTD), public:: retire_rank                        !destination process rank to which the retired bytecode is going
        end type tavp_mng_retirer_conf_t
 !TAVP-MNG locator:
        type, extends(ds_encoder_t), private:: tavp_mng_locator_t
         integer(INTD), private:: ring_comm                         !MPI communicator of the locating ring (tree level)
         integer(INTD), private:: ring_send=-1                      !MPI process rank in the locating ring to which instructions are sent
         integer(INTD), private:: ring_recv=-1                      !MPI process rank in the locating ring from which instructions are received
         type(pack_env_t), private:: bytecode                       !circulating bytecode
         contains
          procedure, public:: configure=>TAVPMNGLocatorConfigure    !configures TAVP-MNG locator
          procedure, public:: start=>TAVPMNGLocatorStart            !start TAVP-MNG locator
          procedure, public:: shutdown=>TAVPMNGLocatorShutdown      !shuts down TAVP-MNG locator
          procedure, public:: encode=>TAVPMNGLocatorEncode          !encodes a DS instruction into the DS bytecode
          procedure, public:: send=>TAVPMNGLocatorSend              !sends a packet of DS instructions to the next process in the locating ring
          procedure, public:: receive=>TAVPMNGLocatorReceive        !receives a packet of DS instructions from the previous process in the locating ring
          procedure, public:: locate=>TAVPMNGLocatorLocate          !scans DS instructions and fills in information for own tensors
        end type tavp_mng_locator_t
 !TAVP_MNG locator configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_locator_conf_t
         integer(INTD), public:: ring_comm                           !MPI communicator of the locating ring (tree level)
         integer(INTD), public:: ring_send                           !MPI process rank in the locating ring to which instructions are sent
         integer(INTD), public:: ring_recv                           !MPI process rank in the locating ring from which instructions are received
        end type tavp_mng_locator_conf_t
 !TAVP-MNG:
        type, extends(dsvp_t), public:: tavp_mng_t
         type(tens_cache_t), private:: tens_cache                    !tensor argument cache
         type(tavp_mng_decoder_t), private:: decoder                 !DSVU: decodes incoming tensor instructions from the manager
         type(tavp_mng_retirer_t), private:: retirer                 !DSVU: retires processed tensor instructions and sends them back to the manager
         type(tavp_mng_locator_t), private:: locator                 !DSVU: locates metadata for remote tensor arguments
#if 0
         type(tavp_mng_decomposer_t), private:: decomposer           !DSVU: decomposes tensors and tensor instructions into smaller pieces
         type(tavp_mng_dispatcher_t), private:: dispatcher           !DSVU: dispatches ready to be executed tensor instructions to the lower tree level
         type(tavp_mng_replicator_t), private:: replicator           !DSVU: replicates tensor blocks for communicaton avoiding
         type(tavp_mng_collector_t), private:: collector             !DSVU: collects processed bytecode from the lower tree level and updates tables
#endif
         contains
          procedure, public:: configure=>TAVPMNGConfigure            !configures the TAVP-MNG DSVP
        end type tavp_mng_t
 !TAVP-MNG configuration:
        type, extends(dsv_conf_t), public:: tavp_mng_conf_t
         character(:), allocatable, public:: description       !TAVP description
         integer(INTD), public:: tavp_id                       !TAVP id
         integer(INTD), public:: source_comm                   !MPI communicator of the bytecode source (higher-level manager)
         integer(INTD), public:: source_rank                   !MPI process rank of the bytecode source (higher-level manager)
         integer(INTD), public:: retire_comm                   !MPI communicator of the retired bytecode destination (higher-level manager)
         integer(INTD), public:: retire_rank                   !MPI process rank of the retired bytecode destination (higher-level manager)
         integer(INTD), public:: ring_comm                     !MPI communicator of the locating ring (tree level)
         integer(INTD), public:: ring_send_rank                !MPI process rank in the locating ring to which instructions are sent
         integer(INTD), public:: ring_recv_rank                !MPI process rank in the locating ring from which instructions are received
         integer(INTD), public:: dispatch_comm                 !MPI communicator of the processes to which the bytecode is further dispatched (lower-level)
         integer(INTD), allocatable, public:: dispatch_rank(:) !MPI ranks of the processes to which the bytecode is further dispatched (lower-level)
         integer(INTD), public:: collect_comm                  !MPI communicator of the processes from which the retired bytecode is collected (lower-level)
         integer(INTD), allocatable, public:: collect_rank(:)  !MPI communicator of the processes from which the retired bytecode is collected (lower-level)
        end type tavp_mng_conf_t
!VISIBILITY:
 !tens_entry_mng_t:
        private TensEntryMngCtor
        private TensEntryMngGetOwnerId
        public tens_entry_mng_dtor
        private tens_entry_mng_alloc
 !tens_oprnd_t:
        private TensOprndCtor
        private TensOprndGetTensor
        private TensOprndGetOwnerId
        private TensOprndSetOwnerId
        private TensOprndIsRemote
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndRelease
        private TensOprndDestruct
        public tens_oprnd_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrEncode
        public tens_instr_dtor
 !tavp_mng_decoder_t:
        private TAVPMNGDecoderConfigure
        private TAVPMNGDecoderStart
        private TAVPMNGDecoderShutdown
        private TAVPMNGDecoderDecode
 !tavp_mng_retirer_t:
        private TAVPMNGRetirerConfigure
        private TAVPMNGRetirerStart
        private TAVPMNGRetirerShutdown
        private TAVPMNGRetirerEncode
 !tavp_mng_locator_t:
        private TAVPMNGLocatorConfigure
        private TAVPMNGLocatorStart
        private TAVPMNGLocatorShutdown
        private TAVPMNGLocatorEncode
        private TAVPMNGLocatorSend
        private TAVPMNGLocatorReceive
        private TAVPMNGLocatorLocate
 !tavp_mng_t:
        private TAVPMNGConfigure
!IMPLEMENTATION:
       contains
![tens_entry_mng_t]========================================
        subroutine TensEntryMngCtor(this,tensor,owner,ierr)
!Constructs a <tens_entry_mng_t>. Note move semantics for <tensor>!
         implicit none
         class(tens_entry_mng_t), intent(out):: this          !out: specialized tensor cache entry
         class(tens_rcrsv_t), pointer, intent(inout):: tensor !inout: pointer to an allocated tensor (ownership transfer will occur here!)
         integer(INTD), intent(in):: owner                    !in: tensor owner id (>=0)
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(tensor)) then
          if(owner.ge.0) then
           call this%set_tensor(tensor,errc)
           if(errc.eq.0) then
            tensor=>NULL() !transfer the ownership
            this%owner_id=owner
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
        end subroutine TensEntryMngCtor
!------------------------------------------------------------
        function TensEntryMngGetOwnerId(this,ierr) result(id)
         implicit none
         integer(INTD):: id                          !out: tensor owner id
         class(tens_entry_mng_t), intent(in):: this  !in: specialized tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          id=this%owner_id
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensEntryMngGetOwnerId
!-------------------------------------------
        subroutine tens_entry_mng_dtor(this)
         implicit none
         type(tens_entry_mng_t):: this !inout: specialized tensor cache entry

         call this%nullify_tensor(.TRUE.) !deallocate the tensor component (if it is set)
         this%owner_id=-1
         return
        end subroutine tens_entry_mng_dtor
!-------------------------------------------------------------
        function tens_entry_mng_alloc(tens_entry) result(ierr)
!Non-member allocator for tens_entry_mng_t.
         implicit none
         integer(INTD):: ierr
         class(tens_cache_entry_t), allocatable, intent(out):: tens_entry

         allocate(tens_entry_mng_t::tens_entry,STAT=ierr)
         return
        end function tens_entry_mng_alloc
![tens_oprnd_t]=========================================
        subroutine TensOprndCtor(this,tensor,ierr,owner)
!Constructs a tensor operand. The <tensor> must be set.
!The owner id is optional.
         implicit none
         class(tens_oprnd_t), intent(inout):: this        !inout: undefined tensor operand (on entrance)
         class(tens_rcrsv_t), target, intent(in):: tensor !in: defined tensor
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: owner      !in: tensor owner id (no restrictions)
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%tensor=>tensor
             if(present(owner)) this%owner_id=owner
             call this%mark_active(errc)
             if(errc.ne.DSVP_SUCCESS) errc=-5
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
        end subroutine TensOprndCtor
!------------------------------------------------------------
        function TensOprndGetTensor(this,ierr) result(tens_p)
!Returns a pointer to the tensor.
         implicit none
         class(tens_rcrsv_t), pointer:: tens_p       !out: pointer to the tensor
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; tens_p=>this%tensor
         if(.not.associated(tens_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetTensor
!---------------------------------------------------------
        function TensOprndGetOwnerId(this,ierr) result(id)
!Returns the tensor owner id.
         implicit none
         integer(INTD):: id                          !out: tensor owner id (negative value means self)
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(this%tensor)) then
          id=this%owner_id
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetOwnerId
!------------------------------------------------------
        subroutine TensOprndSetOwnerId(this,owner,ierr)
!Sets the tensor owner id.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(in):: owner           !in: tensor owner id (no restrictions)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           this%owner_id=owner
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetOwnerId
!--------------------------------------------------------
        function TensOprndIsRemote(this,ierr) result(res)
!Returns TRUE if the tensor operand is remote, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,id

         res=.TRUE. !assume remote by default
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           id=this%get_owner_id(errc)
           if(errc.eq.0) then
            res=(id.ne.role_rank)
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
        end function TensOprndIsRemote
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !No local resources are currently needed
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

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(.not.this%is_present(errc)) then
            if(errc.eq.DSVP_SUCCESS) then
             if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
              call this%set_comm_stat(DS_OPRND_FETCHING,errc); if(errc.ne.DSVP_SUCCESS) errc=-6
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
        end subroutine TensOprndPrefetch
!----------------------------------------------
        subroutine TensOprndUpload(this,ierr)
!Starts uploading the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_present(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
            call this%set_comm_stat(DS_OPRND_UPLOADING,errc); if(errc.ne.DSVP_SUCCESS) errc=-4
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
        end subroutine TensOprndUpload
!---------------------------------------------------------
        function TensOprndSync(this,ierr,wait) result(res)
!Synchronizes a pending prefetch/upload.
         implicit none
         logical:: res                               !out: TRUE on communication completion, FALSE otherwise
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            call this%mark_delivered(errc); if(errc.ne.DSVP_SUCCESS) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!---------------------------------------------
        subroutine TensOprndRelease(this,ierr)
!Releases local resources acquired for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: delivered

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if((.not.delivered).or.(errc.ne.0)) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         !No local resources are currently needed
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

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if((.not.delivered).or.(errc.ne.0)) errc=-4
           endif
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
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
         if(errc.ne.0) call quit(errc,'#FATAL(TAVP-MNG:tens_oprnd_dtor): Destructor failed!')
        end subroutine tens_oprnd_dtor
![tens_instr_t]============================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of a given tensor operation
!for a specific TAVP kind.
         implicit none
         class(tens_instr_t), intent(inout):: this        !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(in):: op_code              !in: instruction code (see top of this module)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         class(*), intent(in), target, optional:: op_spec !in: formal operation specification
         integer(INTD):: errc,ier

         if(this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Construct the instruction:
           select case(op_code)
           case(TAVP_INSTR_NOOP)
           case(TAVP_INSTR_STOP)
           case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
            call construct_instr_create_destroy(errc); if(errc.ne.0) errc=-6
           case(TAVP_INSTR_CONTRACT)
            call construct_instr_contract(errc); if(errc.ne.0) errc=-5
           case default
            errc=-4 !invalid instruction opcode (or not implemented)
           end select
!Activate the instruction:
           if(errc.eq.0) then
            call this%activate(op_code,errc); if(errc.ne.0) errc=-3
           endif
           if(errc.ne.0) then
            call this%set_status(DS_INSTR_RETIRED,ier,TAVP_ERR_GEN_FAILURE)
            call tens_instr_dtor(this) !`Bad: dtor() expects type(tens_instr_t), not class
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine construct_instr_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !op_spec={tens_rcrsv_t}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_oprnd_t), pointer:: tens_oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0; tensor=>NULL()
          select type(op_spec); class is(tens_rcrsv_t); tensor=>op_spec; end select
          if(associated(tensor)) then
           if(tensor%is_set()) then
            call this%alloc_operands(1,jerr)
            if(jerr.eq.DSVP_SUCCESS) then
             allocate(tens_oprnd,STAT=jerr)
             if(jerr.eq.0) then
              call tens_oprnd%tens_oprnd_ctor(tensor,jerr) !`tensor owner is omitted here
              if(jerr.eq.0) then
               oprnd=>tens_oprnd
               call this%set_operand(0,oprnd,jerr); if(jerr.ne.DSVP_SUCCESS) jerr=-6
               oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
               this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/)
              else
               jerr=-5
              endif
              tens_oprnd=>NULL()
             else
              jerr=-4
             endif
            else
             jerr=-3
            endif
           else
            jerr=-2
           endif
           tensor=>NULL() !<tensor> pointed to an external object
          else
           jerr=-1
          endif
          return
         end subroutine construct_instr_create_destroy

         subroutine construct_instr_contract(jerr)
          !CONTRACT two tensors into another tensor:
          !op_spec={tens_contraction_t}
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj
          class(ds_oprnd_t), pointer:: oprnd
          class(ds_instr_ctrl_t), pointer:: instr_ctrl
          class(contr_ptrn_ext_t), pointer:: contr_ptrn
          class(tens_contraction_t), pointer:: tens_contr
          class(ctrl_tens_contr_t), pointer:: tens_contr_ctrl
          class(tens_rcrsv_t), pointer:: tensor
          class(tens_oprnd_t), pointer:: tens_oprnd

          jerr=0; tens_contr=>NULL()
          select type(op_spec); class is(tens_contraction_t); tens_contr=>op_spec; end select
          if(associated(tens_contr)) then
           if(tens_contr%is_set()) then
            contr_ptrn=>tens_contr%get_ext_contr_ptrn(jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             allocate(tens_contr_ctrl,STAT=jerr)
             if(jerr.eq.0) then
              call tens_contr_ctrl%ctrl_tens_contr_ctor(contr_ptrn,jerr,tens_contr%get_prefactor()) !contraction pattern is cloned by value
              if(jerr.eq.0) then
               instr_ctrl=>tens_contr_ctrl
               call this%set_control(instr_ctrl,jerr) !ownership transfer for instr_ctrl=tens_contr_ctrl
               if(jerr.eq.DSVP_SUCCESS) then
                call this%alloc_operands(3,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 do jj=0,2
                  tensor=>tens_contr%get_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) exit
                  allocate(tens_oprnd,STAT=jerr); if(jerr.ne.0) exit
                  call tens_oprnd%tens_oprnd_ctor(tensor,jerr); if(jerr.ne.0) exit
                  oprnd=>tens_oprnd; call this%set_operand(jj,oprnd,jerr); if(jerr.ne.DSVP_SUCCESS) exit !ownership transfer for oprnd=tens_oprnd
                  oprnd=>NULL(); tens_oprnd=>NULL(); tensor=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
                 enddo
                 this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/)
                else
                 jerr=-7
                endif
               else
                jerr=-6
               endif
               instr_ctrl=>NULL()
              else
               jerr=-5
              endif
              tens_contr_ctrl=>NULL() !<tens_contr_ctrl> pointer was saved in the tensor instruction and will later be deallocated
             else
              jerr=-4
             endif
            else
             jerr=-3
            endif
            contr_ptrn=>NULL()
           else
            jerr=-2
           endif
           tens_contr=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine construct_instr_contract

        end subroutine TensInstrCtor
!---------------------------------------------------------
        subroutine TensInstrEncode(this,instr_packet,ierr)
!Encodes a tensor instruction into the bytecode packet:
! 1. Instruction code;
! 2. Instruction status;
! 3. Instruction error code;
! 4. Instruction control field;
! 5. Instruction operands.
         implicit none
         class(tens_instr_t), intent(in):: this          !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,op_code,stat,err_code

!Pack the instruction attributes (op_code,status,error):
         if(.not.this%is_empty(errc)) then
          op_code=this%get_code(errc)
          if(errc.eq.DSVP_SUCCESS) then
           call pack_builtin(instr_packet,op_code,errc)
           if(errc.eq.0) then
            stat=this%get_status(errc,err_code)
            if(errc.eq.DSVP_SUCCESS) then
             call pack_builtin(instr_packet,stat,errc)
             if(errc.eq.0) then
              call pack_builtin(instr_packet,err_code,errc)
              if(errc.eq.0) then
!Pack the instruction body:
               select case(op_code)
               case(TAVP_INSTR_NOOP)
               case(TAVP_INSTR_STOP)
               case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
                call encode_instr_create_destroy(errc); if(errc.ne.0) errc=-9
               case(TAVP_INSTR_CONTRACT)
                call encode_instr_contract(errc); if(errc.ne.0) errc=-8
               case default
                errc=-7 !invalid instruction opcode (or not implemented)
               end select
              else
               errc=-6
              endif
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

        contains

         subroutine encode_instr_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !Packet format: {op_code|status|error|tensor}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0
          oprnd=>this%get_operand(0,jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(oprnd)
           class is(tens_oprnd_t)
            tensor=>oprnd%get_tensor(jerr)
            if(jerr.eq.0) then
             if(tensor%is_set()) then
              call tensor%pack(instr_packet,jerr)
              if(jerr.ne.0) jerr=-5
             else
              jerr=-4
             endif
             tensor=>NULL()
            else
             jerr=-3
            endif
           class default
            jerr=-2
           end select
           oprnd=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_create_destroy

         subroutine encode_instr_contract(jerr)
          !CONTRACT two tensors: tensor0+=tensor1*tensor2*scalar:
          !Packed format: {op_code|status|error|ctrl_tens_contr_t|tensor0,tensor1,tensor2}
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor
          class(ds_instr_ctrl_t), pointer:: tens_contr_ctrl

          jerr=0
          tens_contr_ctrl=>this%get_control(jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(tens_contr_ctrl)
           class is(ctrl_tens_contr_t)
            call tens_contr_ctrl%pack(instr_packet,jerr); if(jerr.ne.0) jerr=-8
           class default
            jerr=-7
           end select
           if(jerr.eq.0) then
            do jj=0,2
             oprnd=>this%get_operand(jj,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-6; exit; endif
             select type(oprnd)
             class is(tens_oprnd_t)
              tensor=>oprnd%get_tensor(jerr); if(jerr.ne.0) then; jerr=-5; exit; endif
              if(.not.tensor%is_set()) then; jerr=-4; exit; endif !trap
              call tensor%pack(instr_packet,jerr); if(jerr.ne.0) then; jerr=-3; exit; endif
             class default
              jerr=-2; exit
             end select
            enddo
            tensor=>NULL(); oprnd=>NULL()
           endif
           tens_contr_ctrl=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_contract

        end subroutine TensInstrEncode
!---------------------------------------
        subroutine tens_instr_dtor(this)
         implicit none
         type(tens_instr_t):: this !inout: empty or retired tensor instruction
         integer(INTD):: sts,errc

         sts=this%get_status(errc)
         if((sts.eq.DS_INSTR_EMPTY.or.sts.eq.DS_INSTR_RETIRED).and.errc.eq.DSVP_SUCCESS) then
          call this%clean(errc)
          if(errc.ne.DSVP_SUCCESS) call quit(errc,'#FATAL(TAVP-MNG:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          call quit(-1,'#FATAL(TAVP-MNG:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor
![tavp_mng_decoder_t]=====================================
        subroutine TAVPMNGDecoderConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_mng_decoder_conf_t)
          if(conf%source_rank.ge.0) then
           this%source_rank=conf%source_rank
           this%source_comm=conf%source_comm
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecoderConfigure
!------------------------------------------------
        subroutine TAVPMNGDecoderStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this !inout: TAVP-MNG decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier

         errc=0
         if(DEBUG.gt.0) write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder started as DSVU # ",i2)') impir,this%get_id() !debug
         !`Implement
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecoderStart
!---------------------------------------------------
        subroutine TAVPMNGDecoderShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this !inout: TAVP-MNG decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         if(DEBUG.gt.0) write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder stopped as DSVU # ",i2)') impir,this%get_id() !debug
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecoderShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGDecoderDecode(this,ds_instr,instr_packet,ierr)
!Decodes a tensor instruction from a plain bytecode.
!The bytecode format is specified in tens_instr_t.encode().
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this     !inout: TAVP-MNG decoder
         class(ds_instr_t), intent(inout), target:: ds_instr !out: decoded tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet     !in: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,op_code,stat,err_code
         class(dsvp_t), pointer:: dsvp
         class(tens_cache_t), pointer:: arg_cache

         if(ds_instr%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Retrieve the TAVP argument cache:
           arg_cache=>NULL(); dsvp=>this%get_dsvp() !host DSVP
           select type(dsvp); class is(tavp_mng_t); arg_cache=>dsvp%tens_cache; end select
           if(associated(arg_cache)) then
!Extract the instruction attributes (opcode,status,error):
            call unpack_builtin(instr_packet,op_code,errc)
            if(errc.eq.0) then
             call unpack_builtin(instr_packet,stat,errc)
             if(errc.eq.0) then
              call unpack_builtin(instr_packet,err_code,errc)
!Extract the instruction body:
              if(errc.eq.0) then
               select case(op_code)
               case(TAVP_INSTR_NOOP)
               case(TAVP_INSTR_STOP)
               case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
                call decode_instr_create_destroy(errc); if(errc.ne.0) errc=-10
               case(TAVP_INSTR_CONTRACT)
                call decode_instr_contract(errc); if(errc.ne.0) errc=-9
               case default
                errc=-8 !unknown instruction opcode (or not implemented)
               end select
!Activate the instruction:
               if(errc.eq.0) then
                call ds_instr%activate(op_code,errc,stat,err_code)
                if(errc.ne.DSVP_SUCCESS) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
                 errc=-7
                endif
               endif
              else
               errc=-6
              endif
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

         contains

          subroutine decode_instr_create_destroy(jerr)
           !CREATE/DESTROY a tensor
           integer(INTD), intent(out):: jerr
           class(tens_rcrsv_t), pointer:: tensor
           class(tens_oprnd_t), pointer:: tens_oprnd
           class(ds_oprnd_t), pointer:: oprnd
           class(tens_cache_entry_t), pointer:: tens_entry
           class(tens_entry_mng_t), pointer:: tens_mng_entry
           logical:: res

           jerr=0
           tensor=>NULL(); allocate(tensor,STAT=jerr)
           if(jerr.eq.0) then
            call tensor%tens_rcrsv_ctor(instr_packet,jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             tens_entry=>NULL(); tens_entry=>arg_cache%lookup(tensor,jerr)
             if(jerr.eq.0) then
              select case(op_code)
              case(TAVP_INSTR_CREATE) !CREATE a tensor
               if(.not.associated(tens_entry)) then
                res=arg_cache%store(tensor,tens_entry_mng_alloc,jerr,tens_entry_p=tens_entry)
                if(res.and.(jerr.eq.0).and.associated(tens_entry)) then
                 tens_mng_entry=>NULL()
                 select type(tens_entry); class is(tens_entry_mng_t); tens_mng_entry=>tens_entry; end select
                 if(associated(tens_mng_entry)) then
                  call ds_instr%alloc_operands(1,jerr)
                  if(jerr.eq.DSVP_SUCCESS) then
                   tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                   if(jerr.eq.0) then
                    call tens_oprnd%tens_oprnd_ctor(tensor,jerr) !tensor is stored in a tensor cache
                    if(jerr.eq.0) then
                     oprnd=>tens_oprnd; call ds_instr%set_operand(0,oprnd,jerr)
                     if(jerr.ne.DSVP_SUCCESS) then
                      call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE)
                      jerr=-16
                     endif
                    else
                     call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE)
                     jerr=-15
                    endif
                   else
                    call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE)
                    jerr=-14
                   endif
                  else
                   call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE)
                   jerr=-13
                  endif
                 else
                  call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE)
                  jerr=-12
                 endif
                else
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE)
                 jerr=-11
                endif
               else
                call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_ARG_DEFINED)
                jerr=-10
               endif
              case(TAVP_INSTR_DESTROY) !DESTROY a tensor
               if(associated(tens_entry)) then
                deallocate(tensor); tensor=>NULL() !deallocate the temporary tensor
                tens_mng_entry=>NULL()
                select type(tens_entry); class is(tens_entry_mng_t); tens_mng_entry=>tens_entry; end select
                if(associated(tens_mng_entry)) then
                 tensor=>tens_mng_entry%get_tensor() !use the same tensor from the tensor cache
                 call ds_instr%alloc_operands(1,jerr)
                 if(jerr.eq.DSVP_SUCCESS) then
                  tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                  if(jerr.eq.0) then
                   call tens_oprnd%tens_oprnd_ctor(tensor,jerr) !tensor is from tensor cache
                   if(jerr.eq.0) then
                    oprnd=>tens_oprnd; call ds_instr%set_operand(0,oprnd,jerr)
                    if(jerr.ne.DSVP_SUCCESS) then
                     call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE)
                     jerr=-9
                    endif
                   else
                    call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE)
                    jerr=-8
                   endif
                  else
                   call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE)
                   jerr=-7
                  endif
                 else
                  call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE)
                  jerr=-6
                 endif
                else
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE)
                 jerr=-5
                endif
               else
                call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_ARG_UNDEFINED)
                jerr=-4
               endif
              end select
             else
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE)
              jerr=-3
             endif
            else
             call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD)
             jerr=-2
            endif
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE)
            jerr=-1
           endif
           return
          end subroutine decode_instr_create_destroy

          subroutine decode_instr_contract(jerr)
           !CONTRACT two tensors
           integer(INTD), intent(out):: jerr
           class(ds_instr_ctrl_t), pointer:: instr_ctrl
           class(ctrl_tens_contr_t), pointer:: tens_contr_ctrl
           class(tens_rcrsv_t), pointer:: tensor
           class(tens_oprnd_t), pointer:: tens_oprnd
           class(ds_oprnd_t), pointer:: oprnd
           class(tens_cache_entry_t), pointer:: tens_entry
           class(tens_entry_mng_t), pointer:: tens_mng_entry
           integer(INTD):: jj

           jerr=0
           tens_contr_ctrl=>NULL(); allocate(tens_contr_ctrl,STAT=jerr)
           if(jerr.eq.0) then
            call tens_contr_ctrl%unpack(instr_packet,jerr)
            if(jerr.eq.0) then
             instr_ctrl=>tens_contr_ctrl; call ds_instr%set_control(instr_ctrl,jerr)
             if(jerr.eq.DSVP_SUCCESS) then
              call ds_instr%alloc_operands(3,jerr)
              if(jerr.eq.DSVP_SUCCESS) then
               do jj=0,2
                tensor=>NULL(); allocate(tensor,STAT=jerr)
                if(jerr.ne.0) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-12; exit
                endif
                call tensor%tens_rcrsv_ctor(instr_packet,jerr)
                if(jerr.ne.TEREC_SUCCESS) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-11; exit
                endif
                tens_entry=>NULL(); tens_entry=>arg_cache%lookup(tensor,jerr)
                if(jerr.ne.0) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-10; exit
                endif
                if(.not.associated(tens_entry)) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_ARG_UNDEFINED); jerr=-9; exit
                endif
                deallocate(tensor); tensor=>NULL() !deallocate the temporary tensor
                tens_mng_entry=>NULL(); select type(tens_entry); class is(tens_entry_mng_t); tens_mng_entry=>tens_entry; end select
                if(.not.associated(tens_mng_entry)) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-8; exit
                endif
                tensor=>tens_mng_entry%get_tensor() !use the same tensor from the tensor cache
                tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                if(jerr.ne.0) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-7; exit
                endif
                call tens_oprnd%tens_oprnd_ctor(tensor,jerr) !tensor is stored in tensor cache
                if(jerr.ne.0) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-6; exit
                endif
                oprnd=>tens_oprnd; call ds_instr%set_operand(jj,oprnd,jerr)
                if(jerr.ne.DSVP_SUCCESS) then
                 call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-5; exit
                endif
               enddo
              else
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE)
               jerr=-4
              endif
             else
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE)
              jerr=-3
             endif
            else
             call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD)
             jerr=-2
            endif
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE)
            jerr=-1
           endif
           return
          end subroutine decode_instr_contract

        end subroutine TAVPMNGDecoderDecode
![tavp_mng_retirer_t]=====================================
        subroutine TAVPMNGRetirerConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,i,n

         errc=0
         select type(conf)
         type is(tavp_mng_retirer_conf_t)
          if(conf%retire_rank.ge.0) then
           this%retire_rank=conf%retire_rank
           this%retire_comm=conf%retire_comm
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerConfigure
!------------------------------------------------
        subroutine TAVPMNGRetirerStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this !inout: TAVP-MNG retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier

         errc=0
         if(DEBUG.gt.0) write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer started as DSVU # ",i2)') impir,this%get_id() !debug
         !`Implement
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerStart
!---------------------------------------------------
        subroutine TAVPMNGRetirerShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this !inout: TAVP-MNG retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         if(DEBUG.gt.0) write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer stopped as DSVU # ",i2)') impir,this%get_id() !debug
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGRetirerEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this  !inout: retirer
         class(ds_instr_t), intent(in), target:: ds_instr !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet  !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         if(ds_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           call ds_instr%encode(instr_packet,errc); if(errc.ne.0) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerEncode
![tavp_mng_locator_t]=====================================
        subroutine TAVPMNGLocatorConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,i,n

         errc=0
         select type(conf)
         type is(tavp_mng_locator_conf_t)
          if(conf%ring_send.ge.0.and.conf%ring_recv.ge.0) then
           this%ring_comm=conf%ring_comm
           this%ring_send=conf%ring_send
           this%ring_recv=conf%ring_recv
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorConfigure
!------------------------------------------------
        subroutine TAVPMNGLocatorStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier

         errc=0
         if(DEBUG.gt.0) write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator started as DSVU # ",i2)') impir,this%get_id() !debug
         !`Implement
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorStart
!---------------------------------------------------
        subroutine TAVPMNGLocatorShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         if(DEBUG.gt.0) write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator stopped as DSVU # ",i2)') impir,this%get_id() !debug
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGLocatorEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this  !inout: TAVP-MNG locator DSVU
         class(ds_instr_t), intent(in), target:: ds_instr !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet  !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         if(ds_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           call ds_instr%encode(instr_packet,errc); if(errc.ne.0) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorEncode
!-----------------------------------------------
        subroutine TAVPMNGLocatorSend(this,ierr)
!Sends the current bytecode packet to the next locating process.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorSend
!--------------------------------------------------
        subroutine TAVPMNGLocatorReceive(this,ierr)
!Receives a bytecode packet from the previous locating process.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorReceive
!-------------------------------------------------
        subroutine TAVPMNGLocatorLocate(this,ierr)
!Scans DS instructions and fills in information for own tensors.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorLocate
![tavp_mng_t]======================================
        subroutine TAVPMNGConfigure(this,conf,ierr)
!Configures TAVP-MNG DSVP:
! * Configures static DSVU;
! * Allocates and configures dynamic DSVU;
! * Sets up global DSVU table in DSVP;
! * Sets up DSVP description and id;
         implicit none
         class(tavp_mng_t), intent(inout), target:: this !out: configured DSVP (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVP configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,num_units
         type(tavp_mng_decoder_conf_t):: decoder_conf
         type(tavp_mng_retirer_conf_t):: retirer_conf
         type(tavp_mng_locator_conf_t):: locator_conf
#if 0
         type(tavp_mng_decomposer_conf_t):: decomposer_conf
         type(tavp_mng_dispatcher_conf_t):: dispatcher_conf
         type(tavp_mng_replicator_conf_t):: replicator_conf
         type(tavp_mng_collector_conf_t):: collector_conf
#endif

         if(.not.this%is_configured(errc)) then
          if(errc.eq.0) then
           select type(conf)
           type is(tavp_mng_conf_t)
            if(conf%tavp_id.ge.0.and.allocated(conf%description)) then
             num_units=0 !increment by one after each unit configuration
 !Configure static DSVU:
  !Decoder:
             decoder_conf=tavp_mng_decoder_conf_t(conf%source_comm,conf%source_rank)
             call this%decoder%configure(decoder_conf,errc)
             if(errc.eq.0) then
              num_units=num_units+1
  !Retirer:
              retirer_conf=tavp_mng_retirer_conf_t(conf%retire_comm,conf%retire_rank)
              call this%retirer%configure(retirer_conf,errc)
              if(errc.eq.0) then
               num_units=num_units+1
  !Locator:
               locator_conf=tavp_mng_locator_conf_t(conf%ring_comm,conf%ring_send_rank,conf%ring_recv_rank)
               call this%locator%configure(locator_conf,errc)
               if(errc.eq.0) then
                num_units=num_units+1
 !Set up global DSVU table (references to all DSVU):
                call this%alloc_units(num_units,errc)
                if(errc.eq.DSVP_SUCCESS) then
                 call this%set_unit(this%decoder,errc)
                 if(errc.eq.DSVP_SUCCESS) then
                  call this%set_unit(this%retirer,errc)
                  if(errc.eq.DSVP_SUCCESS) then
                   call this%set_unit(this%locator,errc)
                   if(errc.eq.DSVP_SUCCESS) then
 !Set the DSVP id and description:
                    call this%set_description(int(conf%tavp_id,INTL),conf%description,errc)
                    if(errc.ne.DSVP_SUCCESS) errc=-12
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
              else
               errc=-6
              endif
             else
              errc=-5
             endif
            else
             errc=-4
            endif
           class default
            errc=-3
           end select
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(errc.ne.0) call this%destroy()
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGConfigure

       end module tavp_manager
